clear
close all
clc

addpath('../splines');
addpath('../splines/Utilities');
addpath('../solvers');



%% Specify all needed input data
problemConfiguration.bSplineOrder_time = 3;
problemConfiguration.bSplineOrder_space = 3;

problemConfiguration.offset_time_ansatz = [0, 2]; %[0, 2]
problemConfiguration.offset_time_test = [0, 2]; %[0, 2]
problemConfiguration.offset_space_ansatz = [1, 1];
problemConfiguration.offset_space_test = [1, 1];

problemConfiguration.d = 2; % Dimension

% Example 1 (stationarly waves, with analytical solution)
n = 2; m = 2;
problemConfiguration.f_time = {@(t) 2+ 0*t; @(t) t.^2*(n^2 + m^2)*pi^2;};
problemConfiguration.f_space_x = {@(x) sin(n*pi*x), @(x) sin(n*pi*x)};
problemConfiguration.f_space_y = {@(y) sin(m*pi*y), @(y) sin(m*pi*y)};
problemConfiguration.u_0_x = @(x) 0 * x;
problemConfiguration.u_0_y = @(y) 0 * y;
problemConfiguration.u_1_x = @(x) 0 * x;
problemConfiguration.u_1_y = @(y) 0 * y;

problemConfiguration.u_analytical = @(x,y,t) t.^2 .* sin(n*pi*x).*sin(m*pi*y);
problemConfiguration.has_analytical_solution = true;


% Example 2
% problemConfiguration.f_time = {@(t) 0*t};
% problemConfiguration.f_space_x = {@(x) 0*x};
% problemConfiguration.f_space_y = {@(y) 0*y};
% problemConfiguration.u_0_x = @(x) x* (x < 0.5) + (1-x).*(x > 0.5);
% problemConfiguration.u_0_y = @(y)  y* (y < 0.5) + (1-y).*(y > 0.5);
% problemConfiguration.u_1_x = @(x) 0*x;
% problemConfiguration.u_1_y = @(y) 0 *y;
% problemConfiguration.inital_conditions = true;
% problemConfiguration.has_analytical_solution = false;
% problemConfiguration.referenceSolutionPath = '2D-example2';
% load('2D-example2') % Reference solution



plotting = false;

% Specify the solution for plotting and the error calculation

resolution.x = 8;
resolution.y = 8;
resolution.t = 8;



l2error = true; % Calculate the l2 error

tolerance = 1e-8;
maxIt = 40;
exactFlag = true;

space_refinement = 2:5;
for refinementLevel_space = space_refinement
    for refinementLevel_time = refinementLevel_space
        
        problemConfiguration.refinementLevel_space = refinementLevel_space;
        problemConfiguration.refinementLevel_time = refinementLevel_time;
        
        fprintf('Creating problem...')
        problem = create2DWaveProblem(problemConfiguration);
        fprintf(' Done!\n')
        
        % Specify a function handle for the pcg methods
        funA=@(X)( problem.M_space * X * problem.Q_time' ...
            + problem.A_space' * X * problem.D_time' ...
            + problem.A_space * X * problem.D_time ...
             + problem.Q_space * X * problem.M_time');
         
         B = kron(problem.Q_time, problem.M_space) ...
        + kron(problem.D_time, problem.A_space') ...
        + kron(problem.D_time', problem.A_space) ...
        + kron(problem.M_time, problem.Q_space);
        g = problem.rhs(:); 
        
        temp = problem.D_time + problem.D_time';
        full(problem.D_time)
        full(temp)
        plot(eig(temp), refinementLevel_space*ones(size(problem.D_time,1)), '*'), grid on, hold on
        drawnow
  
        
%         

%         condest(problem.Q_time)
%         condest(problem.D_time)
%         condest(problem.M_time)
%         return
        
        for ii=4 % Test all three solutions
            
            %% Sove the linear equation system
            
            rhsfull=reshape(problem.rhs, [numel(problem.rhs) 1]);
            rhsfull=rhsfull/norm(rhsfull);
            
            switch ii
                case 1
                    % Product operator preconditioner
                    
                    problem.precond='optimal';
                    tt=tic;
                    [U, iterCGopt(refinementLevel_space, refinementLevel_time)]= ...
                        pcg_fun4(funA,rhsfull,0*rhsfull,problem,maxIt,tolerance,...
                        size(problem.M_space,1),size(problem.M_time,2),exactFlag);
                    U_cg_optimal=U(:)*norm(problem.rhs(:));
                    timeCGopt(refinementLevel_space, refinementLevel_time) = toc(tt);
                    
             %       [~, solvingErrorCGopt(refinementLevel_space, refinementLevel_time)] = ...
               %         calculate2DSolvingError(problem, U_cg_optimal);
                    
                    fprintf('Time refinement: %d, Space refinement: %d, Time to solve: %f\n', ...
                        refinementLevel_time, refinementLevel_space,  timeCGopt(refinementLevel_space, refinementLevel_time))
                    
                case 2
                    % Lyap operator preconditioner
                    problem.precond='lyap';
                    tt=tic;
                    [U, iterCGlyap(refinementLevel_space, refinementLevel_time)]=...
                        pcg_fun4(funA,rhsfull,0*rhsfull,problem,maxIt,tolerance,size(problem.M_space,1),size(problem.M_time,2),exactFlag);
                    timeCGlyap(refinementLevel_space, refinementLevel_time) = toc(tt);
                    U_cg_lyap=U(:)*norm(problem.rhs(:));
                    
                %    [~, solvingErrorCGlyap(refinementLevel_space, refinementLevel_time)] = ...
                %        calculate2DSolvingError(problem, U_cg_lyap);
                    
                    fprintf('Time refinement: %d, Space refinement: %d, Time to solve: %f\n', ...
                        refinementLevel_time, refinementLevel_space,  timeCGlyap(refinementLevel_space, refinementLevel_time))
                    
                case 3
                    % Galerkin
                    [uu,ss,vv]=svds(problem.rhs,1);
                    rhs1=uu(:,1)*sqrt(ss(1,1));
                    rhs2=vv(:,1)*sqrt(ss(1,1));
                    
                    
                    %      path(path,'./valeria/')
                    
                    info=1;
                    tt=tic;
                    [X1,X2,restot,iterGalerkin(refinementLevel_space, refinementLevel_time)]= ...
                        Galerkin3(problem.M_space,2*problem.A_space, ...
                        problem.Q_space,problem.Q_time, ...
                        (problem.D_time+problem.D_time')/2, ...
                        problem.M_time,rhs1,rhs2,maxIt,tolerance,info);
                    timeGalerkin(refinementLevel_space, refinementLevel_time) = toc(tt);
                    U=X1*X2';
                    U_galerkin = U(:);
               %     [~, solvingErrorGalerkin(refinementLevel_space, refinementLevel_time)] = ...
                  %      calculate2DSolvingError(problem, U_galerkin);
                    
                    fprintf('Galerkin: Time refinement: %d, Space refinement: %d, Time to solve: %f\n', ...
                        refinementLevel_time, refinementLevel_space,  timeGalerkin(refinementLevel_space, refinementLevel_time))
                case 4
                    tic
                    U_naiv = B \ g;
                    norm(B*U_naiv - g)
                    toc
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
                case 4
                    U = U_naiv;
                    fprintf('Testing the naiv solution\n');
            end
            
            %% Plot the solution
            if ~plotting && ~l2error
                continue
            end
            
            sol = get2Dsolution(problem, U, resolution);
            
            
            
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
                switch ii
                    case 1
                        errorCGopt(refinementLevel_space, refinementLevel_time) = calculate2DL2Error(problem, sol);
                    case 2
                        errorCGlyap(refinementLevel_space, refinementLevel_time) = calculate2DL2Error(problem, sol);
                    case 3
                        errorGalerkin(refinementLevel_space, refinementLevel_time) = calculate2DL2Error(problem, sol);
                    case 4
                        errorNaiv(refinementLevel_space, refinementLevel_time) = calculate2DL2Error(problem, sol);
                end
            end
            
        end
    end
end

diag(errorNaiv)
figure
semilogy(diag(errorNaiv), 'o--'), grid on
return

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

