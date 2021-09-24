% Benchmarks the space-time methods with all three solvers againts a
% time-stepping scheme.

clear
close all
clc

addpath(genpath('../source'))

resolution.x = 10; 
resolution.t = 10;

refinements = 1:8; 

splineOrder = 3;

[X, T] = ndgrid(linspace(0,1, 2^resolution.x + 1), ...
    linspace(0,1, 2^resolution.t + 1));

%% Define the problem

dimension = 1;
f_time = [];
f_space = [];
u_0 = @(x) x * (x < 0.5) + (1-x).*(x >= 0.5);
u_1 = [];
mu = 1;

solutionAnalytical = dAlembert1D(u_0, @(x) 0 * x, sqrt(mu), ...
    2^resolution.x + 1, 1, 2^resolution.t +1);
name = '1D-continuous.dat';


%% Compute and test the numerical solutions
for refinement = refinements
    fprintf("\n######################\n")
    fprintf("Computing refinement %d\n", refinement)
    fprintf("######################\n")
    % Space time
    problemConfiguration = defineProblem(dimension, f_time, f_space, u_0, ...
        u_1, mu, refinement, refinement, splineOrder, splineOrder);
    
    fprintf('Creating problem\n')
    problem = createProblem(problemConfiguration);
    
    % Galerkin
    fprintf('Computing Galerkin Solution\n')
    tic;
    [UGalerkin, iterGalerkin(refinement)] = solveProblem(problem, 'galerkin');
    timeGalerkin(refinement) = toc;
    fprintf('Getting the solution\n')
    solutionGalerkin = getSolution(problem, UGalerkin, resolution);
    unknowns(refinement) = length(UGalerkin);
    
    % CG Lyapunov
    fprintf('Computing CG Lyapunov\n')
    tic;
    [UCGLyap, iterCGLyap(refinement)] = solveProblem(problem, 'cg-lyap');
    timeCGLyap(refinement) = toc;
    fprintf('Getting the solution\n')
    solutionCGLyap = getSolution(problem, UCGLyap, resolution);
    
    % CG Optimal
    fprintf('Computing CG Optimal\n')
    tic;
    [UCGOpt, iterCGOpt(refinement)] = solveProblem(problem, 'cg-opt');
    timeCGOpt(refinement) = toc;
    fprintf('Getting the solution\n')
    solutionCGOpt = getSolution(problem, UCGOpt, resolution);
    
    % Backslash
    tic;
    UBackslash = solveProblem(problem, 'backslash');
    timeBackslash(refinement) = toc;
    solutionBackslash = getSolution(problem, UBackslash, resolution);
    
    % Time-Stepping
    fprintf('Creating Time-Stepping problem\n')
    problemTSConfiguration = defineTSProblem(dimension, f_time, f_space, ...
        u_0, u_1, mu, refinement, refinement, splineOrder);
    
    problemTS = createTSProblem(problemTSConfiguration);
    
    fprintf('Computing Time-Stepping\n')
    tic;
    UTS = solveTSProblem(problemTS);
    timeTS(refinement) = toc;
    fprintf('Getting the solution\n')
    solutionTS = getTSSolution(problemTS,UTS, resolution);
    
    % Compute the errors
    fprintf('Computing the L2 errors\n')
    errorGalerkin(refinement) = sqrt(mean( (solutionGalerkin-solutionAnalytical).^2, 'all'));
    errorCGLyap(refinement) = sqrt(mean( (solutionCGLyap-solutionAnalytical).^2, 'all'));
    errorCGOpt(refinement) = sqrt(mean( (solutionCGOpt-solutionAnalytical).^2, 'all'));
    errorBackslash(refinement) = sqrt(mean( (solutionBackslash-solutionAnalytical).^2, 'all'));
    errorTS(refinement) = sqrt(mean( (solutionTS-solutionAnalytical).^2, 'all'));
    
    
end


%% Write the data into a file

% Transpose the data
unknowns = unknowns';
refinements = refinements';
timeGalerkin = timeGalerkin';
errorGalerkin = errorGalerkin';
iterGalerkin = iterGalerkin';

timeCGLyap = timeCGLyap';
errorCGLyap = errorCGLyap';
iterCGLyap = iterCGLyap';

timeCGOpt = timeCGOpt';
errorCGOpt = errorCGOpt';
iterCGOpt = iterCGOpt';

timeBackslash = timeBackslash';
errorBackslash = errorBackslash';

timeTS = timeTS';
errorTS = errorTS';

t = table(refinements, unknowns, timeGalerkin, errorGalerkin, iterGalerkin, ...
    timeCGLyap, errorCGLyap, iterCGLyap, ...
    timeCGOpt, errorCGOpt, iterCGOpt, ...
    timeBackslash, errorBackslash, ...
    timeTS, errorTS)

writetable(t, name, 'Delimiter', '\t')

