clear
close all
clc


addpath('../source');
addpath('../splines');

resolution.x = 11;
resolution.t = 11;

refinements = 1:9;

[X, T] = ndgrid(linspace(0,1, 2^resolution.x + 1), ...
        linspace(0,1, 2^resolution.t + 1));

%% Define the problem

dimension = 1;
f_time = {@(t) 2+0*t, @(t) t.^2};
f_space = {@(x) sin(2*pi*x), @(x) 4*pi^2*sin(2*pi*x)};
u_0 = [];
u_1 = [];
mu = 1;
name = '1D-smooth-5th-order-splines.dat';
u_analytical = @(x,t) t.^2 .* sin(2*pi*x);
solutionAnalytical = u_analytical(X,T);
splineOrder = 5;

problemConfiguration.u_analytical = @(x,t) t.^2 .* sin(2*pi*x);

% dimension = 1;
% f_time = [];
% f_space = [];
% u_0 = @(x) sin(3*pi*x);
% u_1 = [];
% mu = 1;
% name = '1D-smooth';






%% Compute and test the numerical solutions
for refinement = refinements
fprintf("Computing refinement %d\n", refinement)    
% Space time   
problemConfiguration = defineProblem(dimension, f_time, f_space, u_0, ...
    u_1, mu, refinement, refinement, splineOrder, splineOrder); 

problem = createProblem(problemConfiguration);

% Galerkin
tic;
[UGalerkin, iterGalerkin(refinement)] = solveProblem(problem, 'galerkin');
timeGalerkin(refinement) = toc;
solutionGalerkin = getSolution(problem, UGalerkin, resolution);

% CG Lyapunov
tic;
[UCGLyap, iterCGLyap(refinement)] = solveProblem(problem, 'cg-lyap');
timeCGLyap(refinement) = toc;
solutionCGLyap = getSolution(problem, UCGLyap, resolution);

% CG Optimal
tic;
[UCGOpt, iterCGOpt(refinement)] = solveProblem(problem, 'cg-opt');
timeCGOpt(refinement) = toc;
solutionCGOpt = getSolution(problem, UCGOpt, resolution);

% Backslash
tic;
UBackslash = solveProblem(problem, 'backslash');
timeBackslash(refinement) = toc;
solutionBackslash = getSolution(problem, UBackslash, resolution);

% Time-Stepping 
problemTSConfiguration = defineTSProblem(dimension, f_time, f_space, ...
    u_0, u_1, mu, refinement, refinement, splineOrder); 

problemTS = createTSProblem(problemTSConfiguration);

tic;
UTS = solveTSProblem(problemTS);
timeTS(refinement) = toc;

solutionTS = getTSSolution(problemTS,UTS, resolution);

% Compute the errors
errorGalerkin(refinement) = sqrt(mean( (solutionGalerkin-solutionAnalytical).^2, 'all'));
errorCGLyap(refinement) = sqrt(mean( (solutionCGLyap-solutionAnalytical).^2, 'all'));
errorCGOpt(refinement) = sqrt(mean( (solutionCGOpt-solutionAnalytical).^2, 'all'));
errorBackslash(refinement) = sqrt(mean( (solutionBackslash-solutionAnalytical).^2, 'all'));
errorTS(refinement) = sqrt(mean( (solutionTS-solutionAnalytical).^2, 'all'));


end

%% Write the data into a file

% Transpose the data
unknowns = (4.^refinements)';
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
    timeCGOpt, errorCGOpt, iterCGOpt, ...,
    timeBackslash, errorBackslash, timeTS, errorTS)

writetable(t, name, 'Delimiter', '\t')
