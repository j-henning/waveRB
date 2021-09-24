% Benchmarks the space-time methods with all three solvers againts a
% time-stepping scheme.

clear
close all
clc

addpath(genpath('../source'))

resolution.x = 7; 
resolution.y = 7;
resolution.z = 7;
resolution.t = 7;


refinements = 1:5;

splineOrder = 5;



%% Define the problem

dimension = 3;
f_time = [];
f_space = [];
funR = @(x,y,z) sqrt((x-0.5).^2 + (y-0.5).^2 + (z-0.5).^2);  % Returns the radius
u0r = @(r) (1 - 5 * abs(r)).*(abs(r) <= 0.2);
u_0 = @(x,y,z) u0r(funR(x,y,z));
u_1 = [];
mu = .01;
name = '3D-continuousOrder5.dat';

% Compute the analytical solution
x = linspace(0,1, 2^resolution.x + 1);
y = linspace(0,1, 2^resolution.y + 1);
z = linspace(0,1, 2^resolution.z + 1);
t = linspace(0,1, 2^resolution.t + 1);

[X, Y, Z, T] = ndgrid(x,y,z,t);



solutionAnalytical = zeros(2^resolution.x+1, 2^resolution.y+1, ...
    2^resolution.z+1, 2^resolution.t+1);

for i=1:2^resolution.x+1
    for j=1:2^resolution.y+1
        for k=1:2^resolution.z+1
            r = funR(x(i), y(j), z(k));
            solutionAnalytical(i,j,k,:) = dAlembertSphere(r,t,sqrt(mu), u0r);
        end
    end
end



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
    errorGalerkin(refinement) = sqrt(mean( (solutionGalerkin-solutionAnalytical).^2, 'all', 'omitnan'));
    errorCGLyap(refinement) = sqrt(mean( (solutionCGLyap-solutionAnalytical).^2, 'all', 'omitnan'));
    errorCGOpt(refinement) = sqrt(mean( (solutionCGOpt-solutionAnalytical).^2, 'all', 'omitnan'));
    errorTS(refinement) = sqrt(mean( (solutionTS-solutionAnalytical).^2, 'all', 'omitnan'));
    
    
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


timeTS = timeTS';
errorTS = errorTS';

t = table(refinements, unknowns, timeGalerkin, errorGalerkin, iterGalerkin, ...
    timeCGLyap, errorCGLyap, iterCGLyap, ...
    timeCGOpt, errorCGOpt, iterCGOpt, ...
    timeTS, errorTS)

writetable(t, name, 'Delimiter', '\t')


plotSolution(solutionGalerkin , resolution)
