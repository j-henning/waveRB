clear
close all
clc


addpath('../source');
addpath('../splines');

problemConfiguration = defineProblem(1, ... % Dimension
    {@(t) t.^2}, ... % f_time
    {@(x) sin(pi*x)}, ... % f_space
    @(x) sin(3*pi*x), ... % u_0
    @(x) x .* (x-1), ... % u_1
    2, ... % Wave speed, mu = c^2
    9, ... % refinement_time
    9, ... % refinement_space
    4, ...
    4 ...
    ); 

problem = createProblem(problemConfiguration);


U = solveProblem(problem);

resolution.x = 10;
resolution.t = 10;

[solution, splines] = getSolution(problem, U, resolution);

figure
plotSolution(solution, resolution);

fprintf('Passed the 1D space time test\n')