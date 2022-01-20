% Benchmarks the space-time methods with all three solvers againts a
% time-stepping scheme.

clear
close all
clc

addpath(genpath('../source'))

resolution.x = 6; %6
resolution.y = 6;
resolution.z = 6;
resolution.t = 6;

refinementTime = 10; %10
refinementSpace = 6; %6


splineOrder = 3;



%% Define the problem

dimension = 3;
f_time = [];
f_space = [];

cutoff = sqrt(2) / 5;

funR = @(x,y,z) sqrt((x-0.5).^2 + (y-0.5).^2 + (z-0.5).^2);  % Returns the radius
u0r = @(r) (abs(r) <= cutoff);
u_0 = @(x,y,z) u0r(funR(x,y,z));
u_1 = [];
mu = 2;



problemTSConfiguration = defineTSProblem(dimension, f_time, f_space, ...
    u_0, u_1, mu, refinementTime, refinementSpace, splineOrder);

problemTS = createTSProblem(problemTSConfiguration);

fprintf('Computing Time-Stepping\n')
UTS = solveTSProblem(problemTS, 100);
fprintf('Getting the solution\n')
solutionReference = getTSSolution(problemTS,UTS, resolution);
save('3d-discontinuous-reference-mu-2', 'solutionReference', 'resolution', '-v7.3')
