% Benchmarks the space-time methods with all three solvers againts a
% time-stepping scheme.

clear
close all
clc

addpath(genpath('../source'))



refinements = 1:8;
nref = length(refinements);

splineOrder = 3;

errorGalerkin = zeros(nref,1);
errorCGLyap = zeros(nref,1);
errorCGOpt = zeros(nref,1);
errorTS = zeros(nref,1);

timeGalerkin = zeros(nref,1);
timeCGLyap = zeros(nref,1);
timeCGOpt = zeros(nref,1);
timeTS = zeros(nref,1);

iterGalerkin = zeros(nref,1);
iterCGLyap = zeros(nref,1);
iterCGOpt = zeros(nref,1);



%% Define the problem

dimension = 3;
f_time = [];
f_space = [];

cutoff = sqrt(2) / 5;

funR = @(x,y,z) sqrt((x-0.5).^2 + (y-0.5).^2 + (z-0.5).^2);  % Returns the radius
u0r = @(r) (1 - abs(r) / cutoff) .* (abs(r) <= cutoff);
u_0 = @(x,y,z) u0r(funR(x,y,z));
u_1 = [];
mu = 4; % equals c^2


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
    [UGalerkin, iterGalerkin(refinement)] = ...
        solveProblem(problem, 'galerkin', 1e6);
    timeGalerkin(refinement) = toc;
    unknowns(refinement) = length(UGalerkin)
    timeGalerkin
    iterGalerkin

    %     % CG Lyapunov
    fprintf('Computing CG Lyapunov\n')
    tic;
    [UCGLyap, iterCGLyap(refinement)] = solveProblem(problem, 'cg-lyap', 1e6);
    timeCGLyap(refinement) = toc;
    timeCGLyap
    iterCGLyap

    % CG Optimal
    fprintf('Computing CG Optimal\n')
    tic;
    [UCGOpt, iterCGOpt(refinement)] = solveProblem(problem, 'cg-opt', 1e6);
    timeCGOpt(refinement) = toc;
    timeCGOpt
    iterCGOpt


    % Time-Stepping
    fprintf('Creating Time-Stepping problem\n')
    problemTSConfiguration = defineTSProblem(dimension, f_time, f_space, ...
        u_0, u_1, mu, refinement, refinement, splineOrder);

    problemTS = createTSProblem(problemTSConfiguration);

    fprintf('Computing Time-Stepping\n')
    tic;
    UTS = solveTSProblem(problemTS, 100);
    timeTS(refinement) = toc;
    timeTS

    save('onlyResiduum', 'timeGalerkin', 'iterGalerkin', ...
        'timeCGLyap', 'iterCGLyap', ...
        'timeCGOpt', 'iterCGOpt', ...
        'timeTS')
end



% timeGalerkin
% timeCGLyap
% timeCGOpt
% timeTS

iterGalerkin
% iterCGLyap
% iterCGOpt

% save('continuous3D', 'errorGalerkin', 'timeGalerkin', 'iterGalerkin', ...
%     'errorCGLyap', 'timeCGLyap', 'iterCGLyap', ...
%     'errorCGOpt', 'timeCGOpt', 'iterCGOpt', ...
%     'errorTS', 'timeTS')


% refinements = refinements';
% t = table(refinements, errorGalerkin, timeGalerkin, iterGalerkin, ...
%       errorCGLyap, timeCGLyap, iterCGLyap, ...
%       errorCGOpt, timeCGOpt, iterCGOpt, ...
%       errorTS, timeTS)
% writetable(t, 'continuous.dat')
