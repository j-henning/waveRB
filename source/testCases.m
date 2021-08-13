clear
close all
clc


addpath('../splines');
addpath('../splines/Utilities');
% addpath('../solvers');

%% 1D space time problem
problemConfiguration = defineProblem(1, ... % Dimension
    {@(t) t.^2}, ... % f_time
    {@(x) sin(pi*x)}, ... % f_space
    @(x) sin(3*pi*x), ... % u_0
    @(x) x .* (x-1), ... % u_1
    2, ... % Wave speed, mu = c^2
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


%% 3D space time problem
funR = @(x,y,z) sqrt((x-0.5).^2 + (y-0.5).^2 + (z-0.5).^2);  % Returns the radius
problemConfiguration = defineProblem(3, ... % Dimension
    [], ... % f_time
    [], ... % f_space
    @(x,y,z) ( funR(x,y,z) <= 0.2), ... % u_0
    @(x,y,z) ( funR(x,y,z) <= 0.2), ... % u_1
    .1, ... % Wave speed, mu = c^2
    5, ... % refinement_time
    5); % refinement_space

problem = createProblem(problemConfiguration);
U = solveProblem(problem);
resolution.x = 6;
resolution.y = 6;
resolution.z = 6;
resolution.t = 6;

solution = getSolution(problem, U, resolution, splines);
plotSolution(solution, resolution);



%% 1D time stepping problem
problemTSConfiguration = defineTSProblem(1, ... % Dimension
    {@(t) t.^2}, ... % f_time
    {@(x) sin(pi*x)}, ... % f_space
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

