clear
close all
clc

addpath('../../splines');
addpath('../../splines/Utilities');
addpath('../../solvers');

%% Specify all needed input data
problemConfiguration.bSplineOrder_time_ansatz = 1;
problemConfiguration.bSplineOrder_time_test = 3;
problemConfiguration.bSplineOrder_space_ansatz = 1;
problemConfiguration.bSplineOrder_space_test = 3;

problemConfiguration.offset_time_ansatz = [0, 0];
problemConfiguration.offset_time_test = [0, 2];

problemConfiguration.offset_space_ansatz = [0, 0];
problemConfiguration.offset_space_test = [1, 1];

% problemConfiguration.offset_time_ansatz = [0, 0];
% problemConfiguration.offset_time_test = [0, 2];
% 
% problemConfiguration.offset_space_ansatz = [0, 0];
% problemConfiguration.offset_space_test = [1, 1];



problemConfiguration.d = 1; % Dimension

% Example 1
% problemConfiguration.f_time = {@(t) 2+0*t, @(t) t.^2};
% problemConfiguration.f_space = {@(x) sin(2*pi*x), @(x) 4*pi^2*sin(2*pi*x)};
% problemConfiguration.u_0_x = @(x) 0.*x;
% problemConfiguration.u_1_x = @(x) 0.*x;
% 
% problemConfiguration.u_analytical = @(x,t) t.^2 .* sin(2*pi*x);
% problemConfiguration.dAlemebert = false;
% problemConfiguration.has_analytical_solution = true;

% Example 2
% f_time = {@(t) 0*t};
% f_space = {@(x) 0*x};
% u_0_x = @(x) x* (x < 0.5) + (1-x).*(x > 0.5);
% u_1_x = @(x) 0*x;
% inital_conditions = true;
% has_analytical_solution = false;
% load('1D-example2')

% Example 3
problemConfiguration.f_time = {@(t) 0*t};
problemConfiguration.f_space = {@(x) 0*x};
% problemConfiguration.u_0_x = @(x) sin(4*pi*x);
% problemConfiguration.u_0_x = @(x) x * (x < 0.5) + (1-x) * (x >= 0.5);
problemConfiguration.u_0_x = @(x)  1 * (x > 0.25 &&  x < 0.75);
problemConfiguration.u_1_x = @(x) 0*x;
problemConfiguration.inital_conditions = true;



plotting = false;

%% Specify the solution for plotting and the error calculation

ref_plotting = 11; 
x = linspace(0,1,2^ref_plotting+1);
t = linspace(0,1,2^ref_plotting+1);

resolution.x = ref_plotting;
resolution.t = ref_plotting;


l2error = true; % Calculate the l2 error

tolerance = 1e-10;
maxIt = 20000000;
exactFlag = true;

sol_ana = dAlembert1D(problemConfiguration.u_0_x,...
    problemConfiguration.u_1_x, 1, length(x), 1, length(t));

%  [X,T] = ndgrid(x,t);
%     sol_ana = problemConfiguration.u_analytical(X,T);




space_refinement = 2:6;
for refinementLevel_space = space_refinement
    for refinementLevel_time = refinementLevel_space
        
        problemConfiguration.refinementLevel_space = refinementLevel_space;
        problemConfiguration.refinementLevel_time = refinementLevel_time;
        
        fprintf('Creating problem...')
        problem = create1DWaveProblem(problemConfiguration);
        fprintf(' Done!\n')
        
        % Specify a function handle for the pcg methods
%         funA=@(X)( problem.M_space * X * problem.Q_time' ...
%             + problem.A_space' * X * problem.D_time' ...
%             + problem.A_space * X * problem.D_time ...
%             + problem.Q_space * X * problem.M_time');
%         
        for ii=[4] % Specify what to use for solving
            
            %% Sove the linear equation system
            
            rhsfull=reshape(problem.rhs, [numel(problem.rhs) 1]);
            rhsfull=rhsfull/norm(rhsfull);
            
            switch ii
                case 1
                    % Product operator preconditioner
                    
                    problem.precond='optimal';
                    tt=tic;
                    [U, iterCGopt(refinementLevel_space, refinementLevel_time)]= ...
                        pcg_fun4(funA,rhsfull,0*rhsfull,problem,maxIt,tolerance, 1e-2,...
                        size(problem.M_space,1),size(problem.M_time,2),exactFlag);
                    U_cg_optimal=U(:)*norm(problem.rhs(:));
                    timeCGopt(refinementLevel_space, refinementLevel_time) = toc(tt);
                    [~, solvingErrorCGopt(refinementLevel_space, refinementLevel_time)] = ...
                        calculate1DSolvingError(problem, U_cg_optimal)
                    
                    fprintf('Time refinement: %d, Space refinement: %d, Time to solve: %f\n', ...
                        refinementLevel_time, refinementLevel_space,  timeCGopt(refinementLevel_space, refinementLevel_time))
                    
                case 2
                    % Lyap operator preconditioner
                    problem.precond='lyap';
                    tt=tic;
                    [U, iterCGlyap(refinementLevel_space, refinementLevel_time)]=...
                        pcg_fun4(funA,rhsfull,0*rhsfull,problem,maxIt,tolerance,1e-2,size(problem.M_space,1),size(problem.M_time,2),exactFlag);
                    timeCGlyap(refinementLevel_space, refinementLevel_time) = toc(tt);
                    U_cg_lyap=U(:)*norm(problem.rhs(:));
%                     [~, solvingErrorCGlyap(refinementLevel_space, refinementLevel_time)] = ...
%                         calculate1DSolvingError(problem, U_cg_lyap);
                    
                    fprintf('Time refinement: %d, Space refinement: %d, Time to solve: %f\n', ...
                        refinementLevel_time, refinementLevel_space,  timeCGlyap(refinementLevel_space, refinementLevel_time))
                    
                case 3
                    % Galerkin
                    [uu,ss,vv]=svds(problem.rhs,1);
                    rhs1=uu(:,1)*sqrt(ss(1,1));
                    rhs2=vv(:,1)*sqrt(ss(1,1));
                    %                    path(path,'./valeria/')
                    tolG=tolerance;
                    info=1;
                    tt=tic;
                    [X1,X2,restot,iterGalerkin(refinementLevel_space, refinementLevel_time)]= ...
                        Galerkin3(problem.M_space,2*problem.A_space,problem.Q_space,problem.Q_time,...
                        (problem.D_time+problem.D_time')/2,problem.M_time,rhs1,rhs2,maxIt,tolG,1e-2,info);
                    timeGalerkin(refinementLevel_space, refinementLevel_time) = toc(tt);
                    U=X1*X2';
                    U_galerkin=U(:);
                    [~, solvingErrorGalerkin(refinementLevel_space, refinementLevel_time)] = ...
                        calculate1DSolvingError(problem, U_galerkin);
                    
                    fprintf('Galerkin: Time refinement: %d, Space refinement: %d, Time to solve: %f\n', ...
                        refinementLevel_time, refinementLevel_space,  timeGalerkin(refinementLevel_space, refinementLevel_time))
                case 4
                    % Backslash
                    
                    B = kron(problem.D_time, problem.M_space) + kron(problem.M_time, problem.A_space);
                    size(B)
                    
                    tt = tic;
                    U_backslash = B \ problem.rhs(:);
                    timeBackslash(refinementLevel_space, refinementLevel_time) = toc(tt);
                    norm(B*U_backslash - problem.rhs(:))
                    
            end
            
            % Set the resolution to the number of refinement, otherwise the
            % smoothing does not work
            res.x = refinementLevel_space;
            res.t = refinementLevel_time;
            
            
            switch ii
                case 1
                    U = U_cg_optimal;
                    fprintf('Testing the cg_opt solution\n');
                    solCG = get1Dsolution(problem, U, resolution);
                case 2
                    U = U_cg_lyap;
                    fprintf('Testing the cg_lyap solution\n');
                     solLyap = get1Dsolution(problem, U, resolution);
                  
                case 3
                    U = U_galerkin;
                    fprintf('Testing the galerkin solution\n');
                    solGalerkin = get1Dsolution(problem, U, resolution);
                case 4
                    U = U_backslash;
                    fprintf('Testing the backslash solution\n');
                     solMldive = get1Dsolution(problem, U, res);
            end
            
            %% Smoothing (just in time)
            [X, T] = meshgrid(linspace(0,1,2^res.x+1), linspace(0,1,2^res.t+1));
            [Xfine, Tfine] = meshgrid(x,t);
            solMldive = movmean(solMldive,2);
            solMldive = interp2(X,T,solMldive,Xfine,Tfine);
            
            
            %% Calculate the error  
    
            if l2error
                switch ii
                    case 1
                        %                         errorCGopt(refinementLevel_space, refinementLevel_time) = calculate1DL2Error(problem, sol);
                        errorVal(refinementLevel_space, refinementLevel_time) =sqrt(mean( (solCG-sol_ana).^2, 'all'))
                    case 2
                        %                         errorCGlyap(refinementLevel_space, refinementLevel_time) = calculate1DL2Error(problem, sol);
                        errorVal(refinementLevel_space, refinementLevel_time) =sqrt(mean( (solLyap-sol_ana).^2, 'all'))
                    case 3
                        %                         errorGalerkin(refinementLevel_space, refinementLevel_time) = calculate1DL2Error(problem, sol);
                        errorVal(refinementLevel_space, refinementLevel_time)  = sqrt(mean( (solGalerkin-sol_ana).^2, 'all'))
                    case 4
                        errorBackslash(refinementLevel_space, refinementLevel_time)  = sqrt(mean( (solMldive-sol_ana).^2, 'all'))
                end
            end
            
        end
    end
end



% Load the timestepping solution
load('../1D-example3-ts')
load('../1D-example3-ts-error')
load('../1D-example3-ts-time')


figure
subplot(2,2,1)
% semilogy(diag(errorVal), 'ro--', 'LineWidth', 2)
semilogy(diag(errorBackslash), 'g*--', 'LineWidth', 2), grid on, hold on
semilogy(errorTS(1:7), 'm*--', 'LineWidth', 2)
% Calculate the triangle coordinates
triangle_x = [length(errorBackslash)-1, length(errorBackslash), length(errorBackslash)-1,length(errorBackslash)-1];
triangle_y = [errorBackslash(end), 0.5 * errorBackslash(end), 0.5 * errorBackslash(end),errorBackslash(end)];
text(length(errorBackslash)-1.25, 0.75 * errorBackslash(end), '1')
text(length(errorBackslash)-0.5, 0.55 * errorBackslash(end), '1')
semilogy(triangle_x, triangle_y, 'k', 'LineWidth', 2)
xlabel('Refinement (space and time)')
ylabel('L_2 error')
title('1D error decay')
legend('Space Time (mldivide)', 'Timestepping')


subplot(2,2,3)
% semilogy(diag(timeCGlyap), 'ro--', 'LineWidth', 2)
semilogy(diag(timeBackslash), 'g*--', 'LineWidth', 2), grid on, hold on
semilogy(timeTS(1:7), 'm*--', 'LineWidth', 2)
xlabel('Refinement (space and time)')
ylabel('CPU time [s]')
title('1D solving time')
legend('Space Time (mldivide)', 'Timestepping')


subplot(2,2,2)
plot(x,sol_ana(:,1), 'LineWidth', 2), hold on
% plot(x,solLyap(:,1), 'r--', 'LineWidth', 2)
plot(x,solMldive(:,1), 'g:', 'LineWidth', 2)
plot(x,solInt(:,1), 'm--', 'LineWidth', 2) % Timestepping solution
xlabel('x')
ylabel('u(x,0)')
legend('Reference solution', 'Space Time (Mldivide)', 'Timestepping')
title('Solution for T = 0')


subplot(2,2,4)
plot(x,sol_ana(:,end), 'LineWidth', 2), hold on
% plot(x,solLyap(:,end), 'r--', 'LineWidth', 2)
plot(x,solMldive(:,end), 'g:', 'LineWidth', 2)
plot(x,solInt(:,end), 'm--', 'LineWidth', 2) % Timestepping solution
xlabel('x')
ylabel('u(x,0)')
legend('Reference solution', 'Space Time (Mldivide)', 'Timestepping')
title('Solution for T = 1')



pause
figure

% v = VideoWriter('smoothing-example')
% open(v)
for i=1:size(sol_ana,2)
    plot(x,sol_ana(:,i), 'LineWidth', 2), hold on
%     plot(x,solLyap(:,i), 'r--', 'LineWidth', 2)
    plot(x,solMldive(:,i), 'g:', 'LineWidth', 2)

     plot(x,solInt(:,i), 'm--', 'LineWidth', 2)
    xlabel('x')
    ylabel('u(x,0)')
    legend('Reference solution', 'Space Time (Mldivide)', 'Timestepping')
    title(['T = ' num2str(i / size(sol_ana,2))])
    drawnow
    hold off
%     frame = getframe(gcf);
%    writeVideo(v,frame);
end
% close(v)

%
%
% ylabel('x', 'Interpreter', 'LaTeX', 'FontSize', 40)
% xlabel('t', 'Interpreter', 'LaTeX', 'FontSize', 40)
% zlabel('u(x,t)', 'Interpreter', 'LaTeX', 'FontSize', 40)





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
%
% figure
% semilogy(l2err_galerkin, 'LineWidth', 3), grid on


%% Write the whole data into files

refinement = space_refinement';

errorVal = diag(errorVal);
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
semilogy(errorVal, '*--', 'LineWidth', 3), hold on, grid on
semilogy(errorCGlyap, 'd--', 'LineWidth', 3)
semilogy(errorGalerkin, 's--', 'LineWidth', 3)
title('L2 Error')
legend('CG opt', 'CG lyap', 'Galerkin')
xlabel('Space and time refinements')

subplot(3,2,4)
loglog(errorVal, timeCGopt, '*--', 'LineWidth', 3), hold on, grid on
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


% writetable(table(refinement,iterCGopt, timeCGopt, errorCGopt), ...
%     'data/1Dexample1-opt-exact-1-1e-10','Delimiter',' ')
% writetable(table(refinement,iterCGlyap, timeCGlyap, errorCGlyap), ...
%     'data/1Dexample1-lyap-exact-1-1e-10','Delimiter',' ')
% writetable(table(refinement,iterGalerkin, timeGalerkin, errorGalerkin), ...
%     'data/1Dexample1-galerkin-1e-10','Delimiter',' ')
