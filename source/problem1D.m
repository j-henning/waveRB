clear
close all
clc


addpath('../splines');
addpath('../splines/Utilities');
% addpath('../solvers');

%% 1D space time problem
problemConfiguration = defineProblem(1, ... % Dimension
    [], ... % f_time
    [], ... % f_space
    @(x) sin(3*pi*x), ... % u_0
    @(x) x .* (x-1), ... % u_1
    1, ... % Wave speed, mu = c^2
    5, ... % refinement_time
    5); % refinement_space

problem = createProblem(problemConfiguration);

problem = changeWaveSpeed(problem,1);

U = solveProblem(problem);

resolution.x = 6;
resolution.t = 6;

[solution, splines] = getSolution(problem, U, resolution);

figure
plotSolution(solution, resolution);


%% 1D time stepping problem
problemTSConfiguration = defineTSProblem(1, ... % Dimension
    [], ... % f_time
    [], ... % f_space
    @(x) sin(3*pi*x), ... % u_0
    @(x) x .* (x-1), ... % u_1
    1, ... % Wave speed, mu = c^2
    5, ... % refinement_time
    5); % refinement_space

problemTS = createTSProblem(problemTSConfiguration);

UTS = solveTSProblem(problemTS);

resolution.x = 6;
resolution.t = 6;
solutionTS = getTSSolution(problemTS,UTS, resolution);
figure
plotSolution(solutionTS, resolution)

norm(solution-solutionTS) / norm(solution)
figure
surf(solution-solutionTS)
