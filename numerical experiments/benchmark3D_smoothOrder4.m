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

splineOrder = 4;

[X, Y, Z, T] = ndgrid(linspace(0,1, 2^resolution.x + 1), ...
    linspace(0,1, 2^resolution.y + 1), ...
    linspace(0,1, 2^resolution.z + 1), ...
    linspace(0,1, 2^resolution.t + 1));

%% Define the problem

dimension = 3;
mm = 3;
nn = 2;
oo = 1;
f_time = {@(t) 2+ 0.*t; @(t) (t.^2)*(mm^2 + nn^2 + oo^2)*(pi^2)};
% f_space = {@(x,y,z) sin(mm*pi*x).*sin(nn*pi*y).*sin(oo*pi*z); ...
%     @(x,y,z) sin(mm*pi*x).*sin(nn*pi*y).*sin(oo*pi*z)};

f_space = {@(x) sin(mm*pi*x), @(y) sin(nn*pi*y), @(z) sin(oo*pi*z); ...
    @(x) sin(mm*pi*x), @(y) sin(nn*pi*y), @(z) sin(oo*pi*z)};

u_0 = [];
u_1 = [];
u_analytical = @(x,y,z,t) t.^2 .* sin(mm*pi*x).*sin(nn*pi*y).*sin(oo*pi*z);
mu = 1;
solutionAnalytical = u_analytical(X,Y,Z,T);
name = '3D-smoothOrder4.dat';


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
    errorTS(refinement) = sqrt(mean( (solutionTS-solutionAnalytical).^2, 'all'));
    
    
end


%% Write the data into a file

% Transpose the data
unknowns = unknowns';
refinements = refinements';
timeGalerkin = timeGalerkin';
errorGalerkin = errorGalerkin';
iterGalerkin = iterGalerkin';


timeTS = timeTS';
errorTS = errorTS';

t = table(refinements, unknowns, timeGalerkin, errorGalerkin, iterGalerkin, ...
    timeTS, errorTS)

writetable(t, name, 'Delimiter', '\t')

