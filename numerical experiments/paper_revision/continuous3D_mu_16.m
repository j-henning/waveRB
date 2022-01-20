% Benchmarks the space-time methods with all three solvers against a
% time-stepping scheme.

clear
close all
clc

addpath(genpath('../source'))

resolution.x = 6;
resolution.y = 6;
resolution.z = 6;
resolution.t = 6;


refinements = 1:5;
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
mu = 16; % equals c^2


load(strcat(['3d-continuous-reference-mu-' num2str(mu)]) )

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
        solveProblem(problem, 'galerkin', 100, 1e-5,1e-2);
    timeGalerkin(refinement) = toc;
    fprintf('Getting the solution\n')
    solutionGalerkin = getSolution(problem, UGalerkin, resolution);
    unknowns(refinement) = length(UGalerkin);

    % CG Lyapunov
    fprintf('Computing CG Lyapunov\n')
    tic;
    [UCGLyap, iterCGLyap(refinement)] = solveProblem(problem, 'cg-lyap', 100);
    timeCGLyap(refinement) = toc;
    fprintf('Getting the solution\n')
    solutionCGLyap = getSolution(problem, UCGLyap, resolution);

    % CG Optimal
    fprintf('Computing CG Optimal\n')
    tic;
    [UCGOpt, iterCGOpt(refinement)] = solveProblem(problem, 'cg-opt', 100);
    timeCGOpt(refinement) = toc;
    fprintf('Getting the solution\n')
    solutionCGOpt = getSolution(problem, UCGOpt, resolution);

    % Time-Stepping
    fprintf('Creating Time-Stepping problem\n')
    problemTSConfiguration = defineTSProblem(dimension, f_time, f_space, ...
        u_0, u_1, mu, refinement, refinement, splineOrder);

    problemTS = createTSProblem(problemTSConfiguration);

    fprintf('Computing Time-Stepping\n')
    tic;
    UTS = solveTSProblem(problemTS, 100);
    timeTS(refinement) = toc;
    fprintf('Getting the solution\n')
    solutionTS = getTSSolution(problemTS,UTS, resolution);

    % Compute the errors
    fprintf('Computing the L2 errors\n')
    errorGalerkin(refinement) = sqrt(mean( (solutionGalerkin-solutionReference).^2, 'all', 'omitnan'));
    errorCGLyap(refinement) = sqrt(mean( (solutionCGLyap-solutionReference).^2, 'all', 'omitnan'));
    errorCGOpt(refinement) = sqrt(mean( (solutionCGOpt-solutionReference).^2, 'all', 'omitnan'));
    errorTS(refinement) = sqrt(mean( (solutionTS-solutionReference).^2, 'all', 'omitnan'));

    errorGalerkin
    errorCGLyap
    errorCGOpt
    errorTS

    timeGalerkin
    timeCGLyap
    timeCGOpt
    timeTS

    iterGalerkin
    iterCGLyap
    iterCGOpt

end




save(strcat(['continuous3D-mu-' num2str(mu)]), 'errorGalerkin', 'timeGalerkin', 'iterGalerkin', ...
    'errorCGLyap', 'timeCGLyap', 'iterCGLyap', ...
    'errorCGOpt', 'timeCGOpt', 'iterCGOpt', ...
    'errorTS', 'timeTS')

refinements = refinements';
t = table(refinements, errorGalerkin, timeGalerkin, iterGalerkin, ...
    errorCGLyap, timeCGLyap, iterCGLyap, ...
    errorCGOpt, timeCGOpt, iterCGOpt, ...
    errorTS, timeTS)
writetable(t, strcat(['continuous-mu-' num2str(mu)]))


% subplot(1,3,1)
% semilogy(errorGalerkin, 'o--'), hold on, grid on
% semilogy(errorTS, '*:')
% legend('Galerkin', 'Timestepping')
% xlabel('Refinement')
% ylabel('L2 Error')
% title('Error')
%
%
% subplot(1,3,2)
% semilogy(timeGalerkin, 'o--'), hold on, grid on
% semilogy(timeTS, '*:')
% legend('Galerkin', 'Timestepping')
% xlabel('Refinement')
% ylabel('Walltimes')
% title('Walltimes')
%
% subplot(1,3,3)
% semilogy(timeGalerkin, errorGalerkin, 'o--'), hold on, grid on
% semilogy(timeTS, errorTS, '*:')
% legend('Galerkin', 'Timestepping')
% xlabel('Walltime')
% ylabel('L2 Error')
% title('Error over Walltime')



