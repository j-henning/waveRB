clear
close all
clc

addpath('../splines');
addpath('../solvers');


%% Specify all needed input data
% Example 1 (stationarly waves, with analytical solution)
n = 3; m = 2;
f_time = {@(t) 2+ 0*t; @(t) t.^2*(n^2 + m^2)*pi^2;};
f_space_x = {@(x) sin(n*pi*x), @(x) sin(n*pi*x)};
f_space_y = {@(y) sin(m*pi*y), @(y) sin(m*pi*y)};
u_0_x = @(x) 0*x; u_0_y = @(y)  0*y;
u_1_x = @(x) 0*x; u_1_y = @(y) 0 *y;
u_analytical = @(x,y,t) t.^2 .* sin(n*pi*x).*sin(m*pi*y);
has_analytical_solution = true;


% Example 2
% f_time = {@(t) 0*t};
% f_space_x = {@(x) 0*x};
% f_space_y = {@(y) 0*y};
% u_0_x = @(x) x* (x < 0.5) + (1-x).*(x > 0.5); u_0_y = @(y)  y* (y < 0.5) + (1-y).*(y > 0.5);
% u_1_x = @(x) 0*x; u_1_y = @(y) 0 *y;
% inital_conditions = true;
% has_analytical_solution = false;
% load('2D-example2') % Reference solution

tolerance = 1e-8; % For the CG or Galerkin method
bSplineOrder_time = 3;
bSplineOrder_space = 3;

offset_time_ansatz = [0, 2];
offset_time_test = [0, 2];
offset_space_ansatz = [1, 1];
offset_space_test = [1, 1];

plotting = false;

% Specify the solution for plotting and the error calculation
x = linspace(0,1,2^8+1);
y = linspace(0,1,2^8+1);
t = linspace(0,1,2^8+1);


l2error = true; % Calculate the l2 error

space_refinement = 1:3;
for refinementLevel_space = space_refinement
    for refinementLevel_time = refinementLevel_space
        %% Calculate the matrices
        disp('Compute matrices')
        
        [Q_time, M_time, N_time, Q_space, M_space, N_space, rhs,M_space_local,Q_space_local,N_space_local] ...
            = computeMatrices2D(refinementLevel_time, bSplineOrder_time, ...
            refinementLevel_space, bSplineOrder_space, f_time, f_space_x, f_space_y, ...
            u_0_x, u_0_y, u_1_x, u_1_y);
        
        disp('done')
        for ii=1:3 % Test all three solutions
            
            %% Sove the linear equation system
            
            % tic;
            % B = kron(Q_time, M_space) ...
            %     + kron(N_time,N_space')...
            %     + kron(N_time',N_space)...
            %     + kron(M_time, Q_space);
            % B=sparse(B);
            % U = B \ reshape(rhs, [numel(rhs) 1]);
            % time = toc;
            % disp(['Solving time with mldivide: ' num2str(time)])
            
            funA=@(X)( M_space*X*Q_time'+N_space'*X*N_time'+...
                N_space*X*N_time+Q_space*X*M_time');
            funB=@(x) reshape(funA(reshape(x, [nthroot(length(x),3)* nthroot(length(x),3), nthroot(length(x),3)])), ...
                [length(x) 1]); % only works if space ref = time ref
            
            
            
            rhsfull=reshape(rhs, [numel(rhs) 1]);
            rhsfull=rhsfull/norm(rhsfull);
            params.Q_time=Q_time; params.M_time=M_time;
            params.M_space=M_space; params.Q_space=Q_space;
            params.A_space=N_space; params.D_time=N_time;
            params.M_spacelocal=M_space_local;
            params.Q_spacelocal=Q_space_local;
            params.N_spacelocal=N_space_local;
            params.d=2;
            exactflag=1;
            switch ii
                case 1
                    % Product operator preconditioner
                    
                    params.precond='optimal';
                    tt=tic;
                    [U, iterCGopt(refinementLevel_space, refinementLevel_time)]= ...
                        pcg_fun4(funA,rhsfull,0*rhsfull,params,30,tolerance,...
                        size(M_space,1),size(M_time,2),exactflag);
                    U_cg_optimal=U(:)*norm(rhs(:));
                    timeCGopt(refinementLevel_space, refinementLevel_time) = toc(tt);
                    
                    solvingErrorCGopt(refinementLevel_space, refinementLevel_time) = ...
                        norm(funB(U_cg_optimal) - reshape(rhs, [numel(rhs), 1])) / norm(reshape(rhs, [numel(rhs), 1]));
                    
                    fprintf('Time refinement: %d, Space refinement: %d, Time to solve: %f\n', ...
                        refinementLevel_time, refinementLevel_space,  timeCGopt(refinementLevel_space, refinementLevel_time))
                    
                case 2
                    % Lyap operator preconditioner
                    params.precond='lyap';
                    tt=tic;
                    [U, iterCGlyap(refinementLevel_space, refinementLevel_time)]=...
                        pcg_fun4(funA,rhsfull,0*rhsfull,params,30,tolerance,size(M_space,1),size(M_time,2),exactflag);
                    timeCGlyap(refinementLevel_space, refinementLevel_time) = toc(tt);
                    U_cg_lyap=U(:)*norm(rhs(:));
                    
                    solvingErrorCGlyap(refinementLevel_space, refinementLevel_time) = ...
                        norm(funB(U_cg_lyap) - reshape(rhs, [numel(rhs), 1])) / norm(reshape(rhs, [numel(rhs), 1]));
                    
                    fprintf('Time refinement: %d, Space refinement: %d, Time to solve: %f\n', ...
                        refinementLevel_time, refinementLevel_space,  timeCGlyap(refinementLevel_space, refinementLevel_time))
                    
                case 3
                    % Galerkin
                    [uu,ss,vv]=svds(rhs,1);
                    rhs1=uu(:,1)*sqrt(ss(1,1));
                    rhs2=vv(:,1)*sqrt(ss(1,1));
                    
                    
                    path(path,'./valeria/')
                    maxitG=100;
                    tolG=tolerance;
                    info=1;
                    tt=tic;
                    [X1,X2,restot,iterGalerkin(refinementLevel_space, refinementLevel_time)]= ...
                        Galerkin3(M_space,2*N_space,Q_space,Q_time,(N_time+N_time')/2,M_time,rhs1,rhs2,maxitG,tolG,info);
                    timeGalerkin(refinementLevel_space, refinementLevel_time) = toc(tt);
                    U=X1*X2';
                    U_galerkin=U(:);
                    solvingErrorGalerkin(refinementLevel_space, refinementLevel_time) = ...
                        norm(funB(U_galerkin) - reshape(rhs, [numel(rhs), 1])) / norm(reshape(rhs, [numel(rhs), 1]));
                    
                    fprintf('Galerkin: Time refinement: %d, Space refinement: %d, Time to solve: %f\n', ...
                        refinementLevel_time, refinementLevel_space,  timeGalerkin(refinementLevel_space, refinementLevel_time))
            end
            
            
            switch ii
                case 1
                    U = U_cg_optimal;
                    fprintf('Testing the cg_opt solution\n');
                case 2
                    U = U_cg_lyap;
                    fprintf('Testing the cg_lyap solution\n');
                case 3
                    U = U_galerkin;
                    fprintf('Testing the galerkin solution\n');
            end
            
            %% Plot the solution
            if ~plotting && ~l2error
                continue
            end
            tic;
            
            % Recalculate some needed variables
            T_time = InitUniNodes(refinementLevel_time, bSplineOrder_time);
            T_space = InitUniNodes(refinementLevel_space, bSplineOrder_space);
            
            
            nsol_time = length(T_time) - bSplineOrder_time;
            nsol_space = length(T_space) - bSplineOrder_space;
            
            
            dim_time =  nsol_time-(offset_time_ansatz(1)+offset_time_ansatz(2));
            
            u = reshape(full(U)', sqrt(numel(U) / dim_time), sqrt(numel(U) / dim_time), dim_time);
            
            % Preparing boundary for plot
            u_full = zeros(nsol_space, nsol_space, nsol_time);
            
            
            u_full(1+offset_space_ansatz(1):nsol_space-offset_space_ansatz(2),...
                1+offset_space_ansatz(1):nsol_space-offset_space_ansatz(2), ...
                1+offset_time_ansatz(1):nsol_time-offset_time_ansatz(2)) = u;
            
            
            sol = zeros(length(x), length(y), length(t));
            
            total_iterations = numel(u_full);
            iteration = 0;
            
            
            % Precalculate the splines
            splines_space = zeros(size(u_full,1), length(x));
            splines_2_space = zeros(size(u_full,1), length(x));
            
            for i=1:size(splines_space,1)
                for j=1:size(splines_space,2)
                    g = @(s) Ndiff(T_space, bSplineOrder_space,0,i,(s==0)*(s+eps) + (s==1)*(s-eps) + (s > 0 && s < 1)*s);
                    h = @(s) Ndiff(T_space, bSplineOrder_space,2,i,(s==0)*(s+eps) + (s==1)*(s-eps) + (s > 0 && s < 1)*s);
                    
                    splines_space(i,j) = g(x(j));
                    splines_2_space(i,j) = h(x(j));
                end
            end
            
            splines_time = zeros(size(u_full,3), length(t));
            splines_2_time = zeros(size(u_full,3), length(t));
            
            for i=1:size(splines_time,1)
                for j=1:size(splines_time,2)
                    g = @(s) Ndiff(T_time, bSplineOrder_time,0,i,(s==0)*(s+eps) + (s==1)*(s-eps) + (s > 0 && s < 1)*s);
                    h = @(s) Ndiff(T_time, bSplineOrder_time,2,i,(s==0)*(s+eps) + (s==1)*(s-eps) + (s > 0 && s < 1)*s);
                    
                    splines_time(i,j) = g(t(j));
                    splines_2_time(i,j) = h(t(j));
                end
            end
            
            for i=1:size(u_full,1) % every space node in x
                if mod(i, size(u_full,1)/8) == 0
                    disp(['Plotting preparation: ' num2str(i / size(u_full,1) * 100) '%'])
                end
                spline_x = splines_space(i,:);
                spline_x_2 = splines_2_space(i,:);
                space_x_indices = min(find(spline_x, 1, 'first'), find(spline_x_2, 1, 'first')):...
                    max(find(spline_x, 1, 'last'), find(spline_x_2, 1, 'last'));
                spline_x = spline_x(space_x_indices);
                spline_x_2 = spline_x_2(space_x_indices);
                
                for j=1:size(u_full, 2) % every space node in y
                    
                    spline_y = splines_space(j,:);
                    spline_y_2 = splines_2_space(j,:);
                    space_y_indices = min(find(spline_y, 1, 'first'), find(spline_y_2, 1, 'first')):...
                        max(find(spline_y, 1, 'last'), find(spline_y_2, 1, 'last'));
                    spline_y = spline_y(space_y_indices);
                    spline_y_2 = spline_y_2(space_y_indices);
                    
                    % Precalculate some plinees
                    temp1= kron(spline_y, spline_x);
                    temp2 = kron(spline_y, spline_x_2) ...
                        + kron(spline_y_2, spline_x);
                    
                    for k=1:size(u_full,3) % every time node
                        if abs(u_full(i,j,k)) < eps
                            continue;
                        end
                        
                        spline_t = splines_time(k,:);
                        spline_t_2 = splines_2_time(k,:);
                        time_indices = min(find(spline_t, 1, 'first'), find(spline_t_2, 1, 'first')):...
                            max(find(spline_t, 1, 'last'), find(spline_t_2, 1, 'last'));
                        spline_t = spline_t(time_indices);
                        spline_t_2 = spline_t_2(time_indices);
                        
                        % Build the whole spline
                        
                        if  k < 5 || k > size(u_full,3) - 5
                            spline = kron(spline_t_2, temp1) ...
                                - kron(spline_t, temp2);
                            
                            
                            spline = reshape(spline, ...
                                [length(space_x_indices), length(space_y_indices), ...
                                length(time_indices)]);
                        end
                        sol(space_x_indices, space_y_indices, time_indices) = ...
                            sol(space_x_indices, space_y_indices, time_indices) + ...
                            u_full(i,j,k) .* spline;
                    end
                end
            end
            toc
            
            
            
            
            if plotting
                zlimits = [min(sol, [], 'all'), max(sol, [], 'all')];
                [X, Y] = meshgrid(y, x);
                fh = figure();
                fh.WindowState = 'maximized';
                for i=1:2:length(t)
                    if has_analytical_solution
                        subplot(1,2,1)
                    end
                    surf(X, Y, sol(:,:,i))
                    xlabel('x')
                    ylabel('y')
                    title(['Numerical solution at t = ' num2str(t(i))])
                    zlim(zlimits)
                    
                    if has_analytical_solution
                        subplot(1,2,2)
                        surf(X,Y, u_analytical(X, Y, t(i)))
                        xlabel('x')
                        ylabel('y')
                        title(['Analytical solution at t = ' num2str(t(i))])
                        zlim(zlimits)
                    end
                    drawnow
                    pause(0.2)
                end
                
                %             if has_analytical_solution
                %                 figure
                %                 plot(t, u_analytical(1/6, .25, t), 'k'), hold on
                %                 plot(t, squeeze(sol(round(size(sol,1) / 4),  round(size(sol,2) / 6), :)))
                %                 plot(t, squeeze(sol(round(size(sol,1) / 4),  round(size(sol,2)*5 / 6), :)))
                %                 plot(t, squeeze(sol(round(size(sol,1)*3 / 4),  round(size(sol,2) / 2), :)))
                %                 plot(t, u_analytical(1/2, .25, t), 'k')
                %                 plot(t, squeeze(sol(round(size(sol,1) * 3 / 4),  round(size(sol,2) / 6), :)))
                %                 plot(t, squeeze(sol(round(size(sol,1) / 4),  round(size(sol,2) * 3 / 6), :)))
                %                 plot(t, squeeze(sol(round(size(sol,1) * 3 / 4),  round(size(sol,2) * 5 / 6), :)))
                %                 grid on
                %             end
            end
            
            
            if l2error
                if has_analytical_solution
                    [X, Y, T] = meshgrid(x, y, t);
                    sol_ana = u_analytical(X,Y,T);
                else
                    sol_ana = solRef;
                end
                
                switch ii
                    case 1
                        errorCGopt(refinementLevel_space, refinementLevel_time) = sqrt(mean( (sol-sol_ana).^2, 'all'));
                    case 2
                        errorCGlyap(refinementLevel_space, refinementLevel_time) = sqrt(mean( (sol-sol_ana).^2, 'all'));
                    case 3
                        errorGalerkin(refinementLevel_space, refinementLevel_time) = sqrt(mean( (sol-sol_ana).^2, 'all'));
                end
            end
            
        end
    end
end

% errorCGopt = diag(l2errCG_opt);
% errorCGlyap = diag(l2errCG_lyap);
% errorGalerkin = diag(l2err_galekin);
% 
% solvingErrorCGopt = diag(solvingErrorCGopt);
% solvingErrorCGlyap = diag(solvingErrorCGlyap);
% solvingErrorGalerkin = diag(solvingErrorGalerkin);
% 
% figure
% subplot(1,2,1)
% semilogy(errorCGopt, '*--', 'LineWidth', 3), hold on, grid on
% semilogy(errorCGlyap, 'd--', 'LineWidth', 3)
% semilogy(errorGalerkin, 's--', 'LineWidth', 3)
% title('L2 Error')
% legend('CG opt', 'CG lyap', 'Galerkin')
% xlabel('Space and time refinements')
% 
% subplot(1,2,2)
% semilogy(solvingErrorCGopt, '*--', 'LineWidth', 3), hold on, grid on
% semilogy(solvingErrorCGlyap, 'd--', 'LineWidth', 3)
% semilogy(solvingErrorGalerkin, 's--', 'LineWidth', 3)
% title('||Ax-b||/||b||')
% legend('CG opt', 'CG lyap', 'Galerkin')
% xlabel('Space and time refinements')


% figure
% % subplot(1,2,1)
% surf(l2err)
% xlabel('space refinement')
% ylabel('time refinement')
% title('L2 error')
%
% subplot(1,2,2)
% surf(solvingTimes)
% xlabel('space refinement')
% ylabel('time refinement')
% title('Time to solve')

% figure
% semilogy(l2err', '*--', 'LineWidth', 2)
% grid on
% xlabel('Space refinements')
% legend('Time refinement = 1', 'Time refinement = 2', 'Time refinement = 3',...
%     'Time refinement = 4','Time refinement = 5','Time refinement = 6', ...
%     'Time refinement = 7')
% title('L2 error')


% figure
% subplot(1,2,1)
% semilogy(solvingTimes_cg_optimal', '*--', 'LineWidth', 2)
% grid on
% xlabel('Space refinements')
% legend('Time refinement = 1' , 'Time refinement = 2', 'Time refinement = 3',...
%     'Time refinement = 4','Time refinement = 5','Time refinement = 6', ...
%     'Time refinement = 7','Time refinement = 8','Time refinement = 9')
% title('Time to solve')
% subplot(1,2,2)
% semilogy(iter_opt', '*--', 'LineWidth', 2)
% grid on
% xlabel('Space refinements')
% legend('Time refinement = 1' , 'Time refinement = 2', 'Time refinement = 3',...
%     'Time refinement = 4','Time refinement = 5','Time refinement = 6', ...
%     'Time refinement = 7','Time refinement = 8','Time refinement = 9')
% title('CG Iterations')

%% Write the whole data into files

% space_refinement = space_refinement';
% writetable(table(space_refinement,iter_opt, solvingTimes_cg_opt, l2errCG_opt), ...
%     'data/2Dexample1-opt-exact-1','Delimiter',' ')
% writetable(table(space_refinement,iter_lyap, solvingTimes_cg_lyap, l2errCG_lyap), ...
%     'data/2Dexample1-lyap-exact-1','Delimiter',' ')
% writetable(table(space_refinement,iter_galerkin, solvingTimes_Galerkin, l2err_galekin), ...
%     'data/2Dexample1-galerkin','Delimiter',' ')

% for i=1:100
%    surf(abs(sol(:,:,i) - sol_ana(:,:,i))), drawnow
%    pause(0.1)
% end

% l2err_galekin
%
% figure
% subplot(1,3,1)
% plot(squeeze(max(abs(sol - sol_ana), [], [1 2])))
% subplot(1,3,2)
% surf(sol(:,:,1) - sol_ana(:,:,1))
% subplot(1,3,3)
% surf(sol(:,:,end) - sol_ana(:,:,end))
%
% figure
% for i=1:100
%     surf(sol(:,:,i) ), drawnow
%     pause(0.1)
% end

refinement = space_refinement';

errorCGopt = diag(errorCGopt);
timeCGopt = diag(timeCGopt);
iterCGopt = diag(iterCGopt);

errorCGlyap = diag(errorCGlyap);
timeCGlyap = diag(timeCGlyap);
iterCGlyap = diag(iterCGlyap);

errorGalerkin = diag(errorGalerkin);
timeGalerkin = diag(timeGalerkin);
iterGalerkin = diag(iterGalerkin);

solvingErrorCGopt = diag(solvingErrorCGopt);
solvingErrorCGlyap = diag(solvingErrorCGlyap);
solvingErrorGalerkin = diag(solvingErrorGalerkin);


figure
subplot(3,2,1)
plot(iterCGopt, '*--', 'LineWidth', 3), hold on, grid on
plot(iterCGlyap, 'd--', 'LineWidth', 3)
plot(iterGalerkin, 's--', 'LineWidth', 3)
title('Number of Iterations')
legend('CG opt', 'CG lyap', 'Galerkin')
xlabel('Space and time refinements')

subplot(3,2,2)
semilogy(timeCGopt, '*--', 'LineWidth', 3), hold on, grid on
semilogy(timeCGlyap, 'd--', 'LineWidth', 3)
semilogy(timeGalerkin, 's--', 'LineWidth', 3)
title('Walltime')
legend('CG opt', 'CG lyap', 'Galerkin')
xlabel('Space and time refinements')


subplot(3,2,3)
semilogy(errorCGopt, '*--', 'LineWidth', 3), hold on, grid on
semilogy(errorCGlyap, 'd--', 'LineWidth', 3)
semilogy(errorGalerkin, 's--', 'LineWidth', 3)
title('L2 Error')
legend('CG opt', 'CG lyap', 'Galerkin')
xlabel('Space and time refinements')

subplot(3,2,4)
loglog(errorCGopt, timeCGopt, '*--', 'LineWidth', 3), hold on, grid on
loglog(errorCGlyap, timeCGlyap, 'd--', 'LineWidth', 3)
loglog(errorGalerkin, timeGalerkin,  's--', 'LineWidth', 3)
title('L2 Error vs Walltime')
legend('CG opt', 'CG lyap', 'Galerkin')
xlabel('L2 Error')

subplot(3,2,5)
semilogy(solvingErrorCGopt, '*--', 'LineWidth', 3), hold on, grid on
semilogy(solvingErrorCGlyap, 'd--', 'LineWidth', 3)
semilogy(solvingErrorGalerkin, 's--', 'LineWidth', 3)
title('||Ax-b||/||b||')
legend('CG opt', 'CG lyap', 'Galerkin')
xlabel('Space and time refinements')


