clear
close all
clc

addpath('..')
addpath('../splines');
addpath('../splines/Utilities');


%% Specify all needed input data
problemConfiguration.bSplineOrder_time_ansatz = 1;
problemConfiguration.bSplineOrder_time_test = 3;
problemConfiguration.bSplineOrder_space_ansatz = 1;
problemConfiguration.bSplineOrder_space_test = 3;

problemConfiguration.offset_time_ansatz = [0, 0];
problemConfiguration.offset_time_test = [0, 2];

problemConfiguration.offset_space_ansatz = [0, 0];
problemConfiguration.offset_space_test = [1, 1];


% problemConfiguration.offset_time_ansatz = [0, 2];
% problemConfiguration.offset_time_test = [0, 2];
% problemConfiguration.offset_space_ansatz = [1, 1];
% problemConfiguration.offset_space_test = [1, 1];

problemConfiguration.d = 1; % Dimension


%% Specify the example to test
% Example 1 (with analytical solution)
% problemConfiguration.f_time = {@(t) 2+0*t, @(t) t.^2};
% problemConfiguration.f_space = {@(x) sin(2*pi*x), @(x) 4*pi^2*sin(2*pi*x)};
% problemConfiguration.u_0_x = @(x) 0.*x;
% problemConfiguration.u_1_x = @(x) 0.*x;
% 
% problemConfiguration.u_analytical = @(x,t) t.^2 .* sin(2*pi*x);
% problemConfiguration.dAlemebert = false;
% problemConfiguration.analyticalSolution = true;


% Example 2 (with dAlemebert solution)
problemConfiguration.f_time = {@(t) 0*t};
problemConfiguration.f_space = {@(x) 0*x};
% problemConfiguration.u_0_x = @(x) (x > 0.25 && x < 0.75);
problemConfiguration.u_0_x = @(x) x .* (x < 0.5) + (1-x) .* (x>= 0.5);
problemConfiguration.u_1_x = @(x) 0*x;

problemConfiguration.dAlemebert = true;
problemConfiguration.analyticalSolution = false;

% Example 3 (without a reference solution)
% problemConfiguration.f_time = {@(t) t};
% problemConfiguration.f_space = {@(x) sin(2*pi*x)};
% problemConfiguration.u_0_x = @(x)  (x > 0.25 && x < 0.75);
% problemConfiguration.u_1_x = @(x) 0*x;
% 
% problemConfiguration.dAlemebert = false;
% problemConfiguration.analyticalSolution = false;


%% Specify if the space time solutions should be smoothend
problemConfiguration.smoothing = true;



%% Specify the solution for plotting and the error calculation

ref_plotting = 9;
x = linspace(0,1,2^ref_plotting+1);
t = linspace(0,1,2^ref_plotting+1);

resolution.x = ref_plotting;
resolution.t = ref_plotting;


refinements = 2:6;


methods = [4]; % Specify what to test:
% 1 - CG optimal (NOT Working for the non-optimal case!)
% 2 - CG Lyap (NOT Working for the non-optimal case!)
% 3 - Galerkin (NOT Working for the non-optimal case!)
% 4 - Backslash
% 5 - Timestepping

%% Specify the tolerance for the iterative methods
tolerance = 1e-10;
maxIt = 1e6;
exactFlag = true; % Exact or inexact preconditioning

tolTimestepping = 1e-7; % Tolerance for each time step

%% Specify the parameters for the time stepping reference solution
KTS = 2^(ref_plotting-2); % Number of time steps
space_refinementTS = ref_plotting -2;
splineOrderTS = 3;
toleranceTS = 1e-7;


%% Calculate the reference solution

if problemConfiguration.dAlemebert
    fprintf('Using dAlemebert1D to compute the reference solution')
    sol_ana = dAlembert1D(problemConfiguration.u_0_x,...
        problemConfiguration.u_1_x, 1, length(x), 1, length(t));
                        
                        
elseif problemConfiguration.analyticalSolution
    fprintf('Using the provided analytical solution')
    [X,T] = ndgrid(x,t);
    sol_ana = problemConfiguration.u_analytical(X,T);
else
    fprintf('Using time stepping to compute the reference solution')
    sol_ana = timeStepping1D(problemConfiguration, KTS, ...
        space_refinementTS, resolution, splineOrderTS, toleranceTS);
end


for refinementLevel_space = refinements
    for refinementLevel_time = refinementLevel_space
        
        problemConfiguration.refinementLevel_space = refinementLevel_space;
        problemConfiguration.refinementLevel_time  = refinementLevel_time;
        
        fprintf('Creating problem...')
        problem = create1DWaveProblem(problemConfiguration);
        fprintf(' Done!\n')
        
        
        % Specify a function handle for the pcg methods
        funA=@(X) problem.M_space * X * problem.D_time' ...
            + problem.A_space * X * problem.M_time';
        
        for ii=methods
            
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
                        calculate1DSolvingError(problem, U_cg_optimal);
                    
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
                    if size(B,1) ~= size(B,2)
                        warning('B is not a square matrix!')
                    end
                    
                    tt = tic;
                    U_backslash = B \ problem.rhs(:);
                    timeBackslash(refinementLevel_space, refinementLevel_time) = toc(tt);
                case 5
                    % Time stepping
                    [solTS, timeTS(refinementLevel_space, refinementLevel_time)] ...
                        = timeStepping1D(problemConfiguration, ...
                        2^refinementLevel_time + 1, ...
                        refinementLevel_space, ...
                        resolution, 3, tolTimestepping);
            end
            
            %% Get the solution and compute the error
            switch ii
                case 1
                    fprintf('Getting the cg_opt solution\n');
                    solCGOptimal = get1Dsolution(problem, U_cg_optimal, resolution);
                    errorCGOptimal(refinementLevel_space, refinementLevel_time) =sqrt(mean( (solCGOptimal-sol_ana).^2, 'all'));
                case 2
                    fprintf('Getting the cg_lyap solution\n');
                    solCGLyap = get1Dsolution(problem, U_cg_lyap, resolution);
                    errorCGLyap(refinementLevel_space, refinementLevel_time) =sqrt(mean( (solCGLyap-sol_ana).^2, 'all'));
                case 3
                    fprintf('Getting the galerkin solution\n');
                    solGalerkin = get1Dsolution(problem, U_galerkin, resolution);
                    errorGalerkin(refinementLevel_space, refinementLevel_time)  = sqrt(mean( (solGalerkin-sol_ana).^2, 'all'));
                case 4
                    fprintf('Getting the backslash solution\n');
                    if problemConfiguration.smoothing
                        res.x = refinementLevel_space;
                        res.t = refinementLevel_time;
                        
                        xCoarse = linspace(0,1,2^res.x+1);
%                         xCoarse = 0.5 * (xCoarse(1:end-1) + xCoarse(2:end));
                        tCoarse = linspace(0,1,2^res.t+1);
%                         tCoarse = 0.5 * (tCoarse(1:end-1) + tCoarse(2:end));
                        [solBackslash] = get1Dsolution(problem, U_backslash, res);
%                         
%                         % Compute the cell average
                        solSmooth = zeros(size(solBackslash,1)+1, size(solBackslash,2)+1);
                        for i = 2:size(solSmooth,1)-1 
                            for j = 2:size(solSmooth,2)-1 
                               solSmooth(i,j) = 0.25 * ...
                                   (solBackslash(i-1,j-1) ...
                                   + solBackslash(i-1,j) ...
                                   + solBackslash(i,j-1) ...
                                   + solBackslash(i,j));
                            end
                        end
                        
              % Interp2 does not extrapolate, so we have to do it
                        % by hand before we feed the data into the method
                        solSmooth(:,1) = solSmooth(:,2) + (solSmooth(:,2) - solSmooth(:,3));
%                         
                         solSmooth(:,end) = solSmooth(:,end-1) + (solSmooth(:,end-1) - solSmooth(:,end-2));

                        
                        
%                         [X, T] = ndgrid(0.5 * (xCoarse(1:end-1) + xCoarse(2:end)), ...
%                             0.5 * (tCoarse(1:end-1) + tCoarse(2:end)));
[X, T] = meshgrid(xCoarse, tCoarse);
% [X, T] = meshgrid(xCoarse(2:end-1), tCoarse(2:end-1));
%                         surf(X,T, solBackslash);
                        [Xfine, Tfine] = meshgrid(x,t);
%                         solBackslash(2:end,:) = solBackslash(1:end-1,:);
                         
                         
%                         solBackslash = movmean(solBackslash,2);
%                          solBackslash(end,:) = [];
%                           solBackslash(1,:) = 0 * solBackslash(1,:);
%                            solBackslash(end,:) = 0 * solBackslash(end,:);
                         
%                         F = griddedInterpolant(X, T, solSmooth, 'next');
%                         solBackslash = F(Xfine, Tfine);


                    
                        solBackslash = interp2(X,T,solSmooth,Xfine,Tfine, 'linear');
                        
                       
                        
                    else
                        solBackslash = get1Dsolution(problem, U_backslash, resolution);
                    end
                    
                    errorBackslash(refinementLevel_space, refinementLevel_time)  = sqrt(mean( (solBackslash-sol_ana).^2, 'all', 'omitnan'));
                case 5
                    errorTS(refinementLevel_space, refinementLevel_time) = sqrt(mean((solTS - sol_ana).^2, 'all'));
            end
            
        end
    end
end


err = diag(errorBackslash)

conv_rate = log(err(end)/ err(end-1)) / log(1/2)


%% Plotting

% Build the legend

minError = Inf; % For the slope triangle

counter = 1;
labelsLong{1} = 'Reference solution';
if ismember(1,methods)
    labels{counter} = 'CG optimal';
    labelsLong{counter+1} = 'CG optimal';
    counter = counter + 1;
end
if ismember(2,methods)
    labels{counter} = 'CG lyap';
    labelsLong{counter+1} = 'CG lyap';
    counter = counter + 1;
end
if ismember(3,methods)
    labels{counter} = 'Galerkin';
    labelsLong{counter+1} = 'Galerkin';
    counter = counter + 1;
end
if ismember(4,methods)
    labels{counter} = 'Backslash';
    labelsLong{counter+1} = 'Backslash';
    counter = counter + 1;
end
if ismember(5,methods)
    labels{counter} = 'Time stepping';
    labelsLong{counter+1} = 'Time stepping';
    counter = counter + 1;
end

figure
subplot(2,3,2)
plot(x,sol_ana(:,1), 'LineWidth', 2), hold on, grid on

subplot(2,3,3)
index = round(size(sol_ana,1)/3);
plot(x,sol_ana(:,index), 'LineWidth', 2), hold on, grid on

subplot(2,3,5)
index = round(size(sol_ana,1)/1.5);
plot(x,sol_ana(:,index), 'LineWidth', 2), hold on, grid on

subplot(2,3,6)
plot(x,sol_ana(:,end), 'LineWidth', 2), hold on, grid on

for m=methods
    subplot(2,3,1)
    switch m
        case 1
            semilogy(diag(errorCGopt), '*--', 'LineWidth', 2)
            minError = min(minError, min(errorCGopt, [], 'all'));
        case 2
            semilogy(diag(errorCGlyap), '*--', 'LineWidth', 2)
            minError = min(minError, min(errorCGLyap, [], 'all'));
        case 3
            semilogy(diag(errorGalerkin), '*--', 'LineWidth', 2)
            minError = min(minError, min(errorGalerkin, [], 'all'));
        case 4
            semilogy(diag(errorBackslash), '*--', 'LineWidth', 2)
            minError = min(minError, min(diag(errorBackslash), [], 'all'));
        case 5
            semilogy(diag(errorTS), '*--', 'LineWidth', 2)
            minError = min(minError, min(diag(errorTS), [], 'all'));
    end
    grid on, hold on
    
    
    subplot(2,3,4)
    switch m
        case 1
            semilogy(diag(timeCGOptimal), '*--', 'LineWidth', 2)
        case 2
            semilogy(diag(timeCGLyap), '*--', 'LineWidth', 2)
        case 3
            semilogy(diag(timeGalerkin), '*--', 'LineWidth', 2)
        case 4
            semilogy(diag(timeBackslash), '*--', 'LineWidth', 2)
        case 5
            semilogy(diag(timeTS), '*--', 'LineWidth', 2)
    end
    grid on, hold on
    
    subplot(2,3,2)
    switch m
        case 1
            plot(x,solCGOptimal(:,1), 'LineWidth', 2)
        case 2
            plot(x,solCGLyap(:,1), 'LineWidth', 2)
        case 3
            plot(x,solGalerkin(:,1), 'LineWidth', 2)
        case 4
            plot(x,solBackslash(:,1), 'LineWidth', 2)
        case 5
            plot(x,solTS(:,1), 'LineWidth', 2)
    end
    
    subplot(2,3,3)
    switch m
        case 1
            index = round(size(solCGOptimal,1)/3);
            plot(x,solCGOptimal(:,index), 'LineWidth', 2)
        case 2
            index = round(size(solCGLyap,1)/3);
            plot(x,solCGLyap(:,index), 'LineWidth', 2)
        case 3
            index = round(size(solGalerkin,1)/3);
            plot(x,solGalerkin(:,index), 'LineWidth', 2)
        case 4
            index = round(size(solBackslash,1)/3);
            plot(x,solBackslash(:,index), 'LineWidth', 2)
        case 5
            index = round(size(solTS,1)/3);
            plot(x,solTS(:,index), 'LineWidth', 2)
    end
    
    subplot(2,3,5)
    switch m
        case 1
            index = round(size(solCGOptimal,1)/1.5);
            plot(x,solCGOptimal(:,index), 'LineWidth', 2)
        case 2
            index = round(size(solCGLyap,1)/1.5);
            plot(x,solCGLyap(:,index), 'LineWidth', 2)
        case 3
            index = round(size(solGalerkin,1)/1.5);
            plot(x,solGalerkin(:,index), 'LineWidth', 2)
        case 4
            index = round(size(solBackslash,1)/1.5);
            plot(x,solBackslash(:,index), 'LineWidth', 2)
        case 5
            index = round(size(solTS,1)/1.5);
            plot(x,solTS(:,index), 'LineWidth', 2)
    end
    
    subplot(2,3,6)
    switch m
        case 1
            plot(x,solCGOptimal(:,end), 'LineWidth', 2)
        case 2
            plot(x,solCGLyap(:,end), 'LineWidth', 2)
        case 3
            plot(x,solGalerkin(:,end), 'LineWidth', 2)
        case 4
            plot(x,solBackslash(:,end), 'LineWidth', 2)
        case 5
            plot(x,solTS(:,end), 'LineWidth', 2)
    end
    
    
    
end



subplot(2,3,4)
legend(labels)
title('CPU time')
xlabel('Refinements')
ylabel('CPU time [s]')

subplot(2,3,2)
legend(labelsLong)
title('Solution for t=0')
xlabel('x')
ylabel('u(x,0)')

subplot(2,3,3)
legend(labelsLong)
title('Solution for t=0.33')
xlabel('x')
ylabel('u(x,0)')

subplot(2,3,5)
legend(labelsLong)
title('Solution for t=0.66')
xlabel('x')
ylabel('u(x,0)')

subplot(2,3,6)
legend(labelsLong)
title('Solution for t=1')
xlabel('x')
ylabel('u(x,0)')


labels{end+1} = 'Linear convergence';
subplot(2,3,1)
% Plot the slope triangle
triangle_x = [refinements(end)-1, refinements(end), refinements(end)-1,refinements(end)-1];
triangle_y = [minError, 0.5 * minError, 0.5 * minError,minError];
semilogy(triangle_x, triangle_y, 'k', 'LineWidth', 2)
text(refinements(end)-1.25, 0.75 * minError, '1')
text(refinements(end)-0.5, 0.55 * minError, '1')
legend(labels)
title('Error decay')
xlabel('Refinements')
ylabel('L_2 error')




