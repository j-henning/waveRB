clear
close all
clc

addpath('../splines');
addpath('../solvers');



%% Specify all needed input data
% Example 1 (very smooth)
mm = 3; nn = 2; oo = 1;
f_time = {@(t) 2+ 0.*t; @(t) (t.^2)*(mm^2 + nn^2 + oo^2)*(pi^2)};
f_space_x = {@(x) sin(mm*pi*x); @(x) sin(mm*pi*x)};
f_space_y = {@(y) sin(nn*pi*y); @(y) sin(nn*pi*y)};
f_space_z = {@(z) sin(oo*pi*z); @(z) sin(oo*pi*z)};
u_analytical = @(x,y,z,t) t.^2 .* sin(mm*pi*x).*sin(nn*pi*y).*sin(oo*pi*z);
u_0_x = @(x) 0*x; u_0_y = @(y) 0*y; u_0_z = @(z) 0*z;
u_1_x = @(x) 0*x; u_1_y = @(y) 0*y; u_1_z = @(z) 0*z;
has_analytical_solution = true;

% Example 2
% f_time = {@(t) 0*t};
% f_space_x = {@(x) 0*x};
% f_space_y = {@(y) 0*y};
% f_space_z = {@(z) 0*z};
% u_0_x = @(x) x* (x < 0.5) + (1-x).*(x > 0.5);
% u_0_y = @(y) y* (y < 0.5) + (1-y).*(y > 0.5);
% u_0_z = @(z) z* (z < 0.5) + (1-z).*(z > 0.5);
% u_1_x = @(x) 0*x;
% u_1_y = @(y) 0*y;
% u_1_z = @(z) 0*z;
% inital_conditions = true;
% has_analytical_solution = false;
% load('3D-example2.mat')

tolerance = 1e-8; % For the CG or Galerkin method
bSplineOrder_time = 3;
bSplineOrder_space = 3;

offset_time_ansatz = [0, 2];
offset_time_test = [0, 2];
offset_space_ansatz = [1, 1];
offset_space_test = [1, 1];

plotting = false;
% Specify the resolution for plotting and the error calculation
x = linspace(0,1,2^6+1); % Increase for large space refinement levels
y = x;
z = x;
t = linspace(0,1,2^6+1);

l2error = true; % Calculate the l2 error


for refinementLevel_space = 1:4
    for refinementLevel_time = refinementLevel_space
        
        
        
        %% Calculate the matrices
        [Q_time, M_time, N_time, Q_space, M_space, N_space, rhs,M_space_local,Q_space_local,N_space_local] ...
            = computeMatrices3D(refinementLevel_time, bSplineOrder_time, ...
            refinementLevel_space, bSplineOrder_space, ...
            f_time, f_space_x, f_space_y, f_space_z, ...
            u_0_x, u_0_y, u_0_z, u_1_x, u_1_y, u_1_z);
        fprintf('\n')
        fprintf('\n')
        disp(['Size of the time matrices: ' num2str(length(Q_time)) ' x ' num2str(length(Q_time))])
        disp(['Size of the space matrices: ' num2str(length(Q_space)) ' x ' num2str(length(Q_space))])
        
        %% Sove the linear equation system
        
        
        % B = kron(Q_time, M_space) ...
        %     + kron(N_time,N_space')...
        %     + kron(N_time',N_space)...
        %     + kron(M_time, Q_space);
        % tic;
        % U = B \ rhs(:);
        % time = toc;
        % disp(['Solving time: ' num2str(time)])
        
        
        
        funA=@(X)( M_space*X*Q_time'+N_space'*X*N_time'+N_space*X*N_time+Q_space*X*M_time');
        funB=@(x) reshape(funA(reshape(x, [nthroot(length(x),4)^3, nthroot(length(x),4)])), ...
            [length(x) 1]); % only works if space ref = time ref
        
        rhsfull=reshape(rhs, [numel(rhs) 1]);
        rhsfull=rhsfull/norm(rhsfull);
        params.Q_time=Q_time; params.M_time=M_time;
        params.M_space=M_space; params.Q_space=Q_space;
        params.A_space=N_space; params.D_time=N_time;
        params.M_spacelocal=M_space_local;
        params.Q_spacelocal=Q_space_local;
        params.N_spacelocal=N_space_local;
        params.d=3;
        exactflag=0;
        
        for ii=1 % Test all three solutions
            switch ii
                case 1
                    % Product operator preconditioner
                    params.precond='optimal';
                    tt=tic;
                    [U, iter_opt(refinementLevel_space, refinementLevel_time)]=pcg_fun4(funA,rhsfull,0*rhsfull,params,30,tolerance,size(M_space,1),size(M_time,2),exactflag);
                    U_cg_optimal=U(:)*norm(rhs(:));
                    solvingTimes_cg_opt(refinementLevel_space, refinementLevel_time) = toc(tt);
                    
                    solvingErrorCGopt(refinementLevel_space, refinementLevel_time) = ...
                        norm(funB(U_cg_optimal) - reshape(rhs, [numel(rhs), 1])) / norm(reshape(rhs, [numel(rhs), 1]));
                    
                    fprintf('Time refinement: %d, Space refinement: %d, Time to solve: %f\n', ...
                        refinementLevel_time, refinementLevel_space,  solvingTimes_cg_opt(refinementLevel_space, refinementLevel_time))
                    
                case 2
                    % Lyap operator preconditioner
                    params.precond='lyap';
                    tt=tic;
                    [U, iter_lyap(refinementLevel_space, refinementLevel_time)]=pcg_fun4(funA,rhsfull,0*rhsfull,params,30,tolerance,size(M_space,1),size(M_time,2),exactflag);
                    solvingTimes_cg_lyap(refinementLevel_space, refinementLevel_time) = toc(tt);
                    U_cg_lyap=U(:)*norm(rhs(:));
                    
                    solvingErrorCGlyap(refinementLevel_space, refinementLevel_time) = ...
                        norm(funB(U_cg_lyap) - reshape(rhs, [numel(rhs), 1])) / norm(reshape(rhs, [numel(rhs), 1]));
                    
                    fprintf('Time refinement: %d, Space refinement: %d, Time to solve: %f\n', ...
                        refinementLevel_time, refinementLevel_space,  solvingTimes_cg_lyap(refinementLevel_space, refinementLevel_time))
                case 3
                    % Galerkin
                    [uu,ss,vv]=svds(rhs,1);
                    rhs1=uu(:,1)*sqrt(ss(1,1));
                    rhs2=vv(:,1)*sqrt(ss(1,1));
                    path(path,'./valeria/')
                    maxitG=50;
                    tolG=tolerance;
                    info=1;
                    tt=tic;
                    [X1,X2,restot,iter_galerkin(refinementLevel_space, refinementLevel_time)]=Galerkin3(M_space,2*N_space,Q_space,Q_time,(N_time+N_time')/2,M_time,rhs1,rhs2,maxitG,tolG,info);
                    solvingTimes_Galerkin(refinementLevel_space, refinementLevel_time) = toc(tt);
                    U=X1*X2';
                    U_galerkin=U(:);
                    
                    solvingErrorGalerkin(refinementLevel_space, refinementLevel_time) = ...
                        norm(funB(U_galerkin) - reshape(rhs, [numel(rhs), 1])) / norm(reshape(rhs, [numel(rhs), 1]));
                    
                    
                    fprintf('Galerkin: Time refinement: %d, Space refinement: %d, Time to solve: %f\n', ...
                        refinementLevel_time, refinementLevel_space,  solvingTimes_Galerkin(refinementLevel_space, refinementLevel_time))
                    
            end
            
            
            %% Plot the solution (only recommended for small problems)
            if ~plotting && ~l2error
                continue;
            else
                
                
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
                
                
                
                % NOTE!
                % This code is very inefficient, due to the spline implementation
                % and could be implemented much more efficient.
                tic;
                
                % Recalculate some needed variables
                T_time = InitUniNodes(refinementLevel_time, bSplineOrder_time);
                T_space = InitUniNodes(refinementLevel_space, bSplineOrder_space);
                
                nsol_time = length(T_time) - bSplineOrder_time;
                nsol_space = length(T_space) - bSplineOrder_space;
                
                
                dim_time =  nsol_time-(offset_time_ansatz(1)+offset_time_ansatz(2));
                
                u = reshape(full(U), round((numel(U) / dim_time)^(1/3)), ...
                    round((numel(U) / dim_time)^(1/3)), ...
                    round((numel(U) / dim_time)^(1/3)), ...
                    dim_time);
                
                % u = permute(u, [2 3 1 4]);
                % Preparing boundary for plot
                
                u_full = zeros(nsol_space, nsol_space, nsol_space, nsol_time);
                
                u_full(1+offset_space_ansatz(1):nsol_space-offset_space_ansatz(2),...
                    1+offset_space_ansatz(1):nsol_space-offset_space_ansatz(2), ...
                    1+offset_space_ansatz(1):nsol_space-offset_space_ansatz(2), ...
                    1+offset_time_ansatz(1):nsol_time-offset_time_ansatz(2)) = u;
                
                
                sol = zeros(length(x), length(y), length(z), length(t));
                
                total_iterations = numel(u_full);
                iteration = 0;
                
                % Precalculate the splines
                splines_space = zeros(size(u_full,1), length(x));
                splines_2_space = zeros(size(u_full,1), length(x));
                
                for i=1:size(splines_space,1)
                    for j=1:size(splines_space,2)
                        g = @(s) Ndiff(T_space, bSplineOrder_space,0,i,(s==1)*(s-eps) + (s < 1)*s);
                        h = @(s) Ndiff(T_space, bSplineOrder_space,2,i,(s==1)*(s-eps) + (s < 1)*s);
                        
                        splines_space(i,j) = g(x(j));
                        splines_2_space(i,j) = h(x(j));
                    end
                end
                
                splines_time = zeros(size(u_full,4), length(t));
                splines_2_time = zeros(size(u_full,4), length(t));
                
                for i=1:size(splines_time,1)
                    for j=1:size(splines_time,2)
                        g = @(s) Ndiff(T_time, bSplineOrder_time,0,i,(s==1)*(s-eps) + (s < 1)*s);
                        h = @(s) Ndiff(T_time, bSplineOrder_time,2,i,(s==1)*(s-eps) + (s < 1)*s);
                        
                        splines_time(i,j) = g(t(j));
                        splines_2_time(i,j) = h(t(j));
                    end
                end
                
                
                for i=1:size(u_full,1) % every space node in x
                    if mod(i, size(u_full,1)/8) == 0
                        disp(['Getting solution: ' num2str(i / size(u_full,1) * 100) '%'])
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
                        
                        for k=1:size(u_full, 3) % every space node in z
                            spline_z = splines_space(k,:);
                            spline_z_2 = splines_2_space(k,:);
                            space_z_indices = min(find(spline_z, 1, 'first'), find(spline_z_2, 1, 'first')):...
                                max(find(spline_z, 1, 'last'), find(spline_z_2, 1, 'last'));
                            spline_z = spline_z(space_z_indices);
                            spline_z_2 = spline_z_2(space_z_indices);
                            
                            spl_space = kron(spline_z, kron(spline_y, spline_x));
                            spl_space_mixed = kron(spline_z, kron(spline_y, spline_x_2))...
                                + kron(spline_z, kron(spline_y_2, spline_x)) ...
                                + kron(spline_z_2, kron(spline_y, spline_x));
                            
                            for l=1:size(u_full,4) % every time node
                                
                                if abs(u_full(i,j,k,l)) < eps
                                    continue;
                                end
                                
                                spline_t = splines_time(l,:);
                                spline_t_2 = splines_2_time(l,:);
                                time_indices = min(find(spline_t, 1, 'first'), find(spline_t_2, 1, 'first')):...
                                    max(find(spline_t, 1, 'last'), find(spline_t_2, 1, 'last'));
                                spline_t = spline_t(time_indices);
                                spline_t_2 = spline_t_2(time_indices);
                                
                                % Build the whole spline
                                %                     spline = kron(spline_t_2, kron(spline_z, kron(spline_y, spline_x))) ...
                                %                         - kron(spline_t, kron(spline_z, kron(spline_y, spline_x_2))) ...
                                %                         - kron(spline_t, kron(spline_z, kron(spline_y_2, spline_x))) ...
                                %                         - kron(spline_t, kron(spline_z_2, kron(spline_y, spline_x)));
                                spline = kron(spline_t_2, spl_space) - kron(spline_t, spl_space_mixed);
                                
                                spline = reshape(spline, ...
                                    [length(space_x_indices), length(space_y_indices), ...
                                    length(space_z_indices), length(time_indices)]);
                                sol(space_x_indices, space_y_indices, space_z_indices, time_indices) = ...
                                    sol(space_x_indices, space_y_indices, space_z_indices, time_indices) + ...
                                    u_full(i,j,k,l) * spline;
                            end
                        end
                    end
                end
                toc
                
                if plotting
                    climits = [-1 1];
                    [X, Y, Z] = meshgrid(x, y, z);
                    
                    % figure
                    % for i=1:length(t)
                    %     subplot(1,2,1)
                    %     temp = reshape((sol(:,:,:,i) - climits(1)) / (climits(2) - climits(1)), [numel(sol(:,:,:,i)) 1]);
                    %     temp(1,3) = 0; % add two more colors, with are zero
                    %     temp(:,2) = 0.5;
                    %     temp(:,3) = 0.5;
                    %     scatter3(X(:), Y(:), Z(:), 100 * temp(:,1).^2 +0.01, temp, 'filled')
                    %     xlabel('x')
                    %     ylabel('y')
                    %     zlabel('z')
                    %     title(['Numerical solution at t = ' num2str(t(i))])
                    %     drawnow
                    %
                    %     subplot(1,2,2)
                    %     temp(:,1) = reshape((u_analytical(X, Y, Z, t(i)) - climits(1)) / (climits(2) - climits(1)), [numel(sol(:,:,:,i)) 1]);
                    %     scatter3(X(:), Y(:), Z(:), 100 * temp(:,1).^2 +0.01, temp, 'filled')
                    %     xlabel('x')
                    %     ylabel('y')
                    %     title(['Analytical solution at t = ' num2str(t(i))])
                    %     drawnow
                    %     pause(1)
                    % end
                    
                    [X, Y, Z, T] = ndgrid(x, y, z, t);
                    sol_ana = u_analytical(X,Y,Z,T);
                    figure
                    for i=1:length(t)
                        % x-y plane
                        subplot(2,3,1)
                        s = sol(:,:,length(z)/2,i);
                        surf(squeeze(s))
                        zlim([min(sol,[], 'all'), max(sol,[],'all')])
                        title(['Numerical solution in the x-y plane at z = 0.5, t = ' num2str(t(i))])
                        
                        subplot(2,3,4)
                        s = sol_ana(:,:,length(z)/2,i);
                        surf(squeeze(s))
                        zlim([min(sol,[], 'all'), max(sol,[],'all')])
                        title(['Analytical solution in the x-y plane at z = 0.5, t = ' num2str(t(i))])
                        
                        % x-z plane
                        subplot(2,3,2)
                        s = sol(:,length(y)/2,:,i);
                        surf(squeeze(s))
                        zlim([min(sol,[], 'all'), max(sol,[],'all')])
                        title(['Numerical solution in the x-z plane at y = 0.5, t = ' num2str(t(i))])
                        
                        subplot(2,3,5)
                        s = sol_ana(:,length(y)/2,:,i);
                        surf(squeeze(s))
                        zlim([min(sol,[], 'all'), max(sol,[],'all')])
                        title(['Analytical solution in the x-z plane at y = 0.5, t = ' num2str(t(i))])
                        
                        % y-z plane
                        subplot(2,3,3)
                        s = sol(length(x)/2,:,:,i);
                        surf(squeeze(s))
                        zlim([min(sol,[], 'all'), max(sol,[],'all')])
                        title(['Numerical solution in the x-z plane at x = 0.5, t = ' num2str(t(i))])
                        
                        subplot(2,3,6)
                        s = sol_ana(length(x)/2,:,:,i);
                        surf(squeeze(s))
                        zlim([min(sol,[], 'all'), max(sol,[],'all')])
                        title(['Analytical solution in the x-z plane at x = 0.5, t = ' num2str(t(i))])
                        
                        drawnow
                    end
                    
                    
                    figure
                    plot(t, u_analytical(0.5,0.5,0.5, t)), hold on
                    plot(t, squeeze(sol(length(x)/2,length(y)/2,length(z)/2,:)))
                    grid on
                    title('Core value')
                    legend('Analytical value', 'Numerical value')
                end
                
                if l2error
                    if has_analytical_solution
                        [X, Y, Z, T] = ndgrid(x, y, z, t);
                        sol_ana = u_analytical(X,Y,Z,T);
                    else
                        sol_ana = solRef;
                    end
                    switch ii
                        case 1
                            l2errCG_opt(refinementLevel_space, refinementLevel_time) = sqrt(mean( (sol-sol_ana).^2, 'all'));
                        case 2
                            l2errCG_lyap(refinementLevel_space, refinementLevel_time) = sqrt(mean( (sol-sol_ana).^2, 'all'));
                        case 3
                            l2err_galekin(refinementLevel_space, refinementLevel_time) = sqrt(mean( (sol-sol_ana).^2, 'all'));
                    end
                    
                    
                end
                
            end
        end
    end
end


errorCGopt = diag(l2errCG_opt);
% errorCGlyap = diag(l2errCG_lyap);
% errorGalerkin = diag(l2err_galekin);

solvingErrorCGopt = diag(solvingErrorCGopt);
% solvingErrorCGlyap = diag(solvingErrorCGlyap);
% solvingErrorGalerkin = diag(solvingErrorGalerkin);

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


%
% figure
% subplot(1,2,1)
% surf(l2err_galekin)
% xlabel('space refinement')
% ylabel('time refinement')
% title('L2 error')
%
% subplot(1,2,2)
% surf(solvingTimes_Galerkin)
% xlabel('space refinement')
% ylabel('time refinement')
% title('Galerkin CPU Time')
%
% figure
% subplot(1,2,1)
% semilogy(l2err_galekin', '*--', 'LineWidth', 2)
% grid on
% xlabel('Space refinements')
% legend('Time refinement = 1' , 'Time refinement = 2', 'Time refinement = 3',...
%     'Time refinement = 4','Time refinement = 5','Time refinement = 6', ...
%     'Time refinement = 7','Time refinement = 8','Time refinement = 9', ...
%     'Time refinement = 10')
% title('L2 error')
%
% subplot(1,2,2)
% semilogy(solvingTimes_Galerkin', '*--', 'LineWidth', 2)
% grid on
% xlabel('Space refinements')
% legend('Time refinement = 1' , 'Time refinement = 2', 'Time refinement = 3',...
%     'Time refinement = 4','Time refinement = 5','Time refinement = 6', ...
%     'Time refinement = 7','Time refinement = 8','Time refinement = 9', ...
%     'Time refinement = 10')
% title('Time to solve')


%
% loglog(reshape(l2errCG_opt, [], 1)', reshape(solvingTimes_cg_optimal,[], 1)', '*'), hold on
% %loglog(reshape(l2errCG_lyap, [], 1), reshape(solvingTimes_cg_lyap,[], 1), 's:')
% loglog(reshape(l2err_galekin, [], 1), reshape(solvingTimes_Galerkin,[], 1), 'd')
% legend('CG opt', 'Galerkin')
% xlabel('L2 Error')
% ylabel('CPU time')
% grid on

%% Write the whole data into files

timeCGopt = diag(solvingTimes_cg_opt);
% timeCGlyap = diag(solvingTimes_cg_lyap);
% timeGalerkin = diag(solvingTimes_Galerkin);

iterCGopt = diag(iter_opt);
% iterCGlyap = diag(iter_lyap);
% iterGalerkin = diag(iter_galerkin);


refinement = (1:4)';
% writetable(table(refinement,iterCGopt, timeCGopt, errorCGopt), ...
%     'data/3Dexample2-opt-exact-0','Delimiter',' ')
% writetable(table(refinement,iterCGlyap, timeCGlyap, errorCGlyap), ...
%     'data/3Dexample1-lyap-exact-0','Delimiter',' ')
% writetable(table(refinement,iterGalerkin, timeGalerkin, errorGalerkin), ...
%     'data/3Dexample1-galerkin','Delimiter',' ')