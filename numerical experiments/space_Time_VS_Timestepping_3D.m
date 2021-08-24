clear
close all
clc


addpath('../source');
addpath('../splines');

resolution.x = 6; %6
resolution.y = 6;
resolution.z = 6;
resolution.t = 6;


refinements = 1:5;

splineOrder = 5;

[X, Y, Z, T] = ndgrid(linspace(0,1, 2^resolution.x + 1), ...
    linspace(0,1, 2^resolution.y + 1), ...
    linspace(0,1, 2^resolution.z + 1), ...
    linspace(0,1, 2^resolution.t + 1));

%% Define the problem

dimension = 3;
mm = 1;
nn = 1;
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
name = '3D-smooth-5th-order.dat';

totalTime = tic;
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
    %     tic;
    %     UBackslash = solveProblem(problem, 'backslash');
    %     timeBackslash(refinement) = toc;
    %     solutionBackslash = getSolution(problem, UBackslash, resolution);
    
    % Time-Stepping
    fprintf('Creating Time-Stepping problem\n')
    problemTSConfiguration = defineTSProblem(dimension, f_time, f_space, ...
        u_0, u_1, mu, refinement, refinement, splineOrder);
    
    problemTS = createTSProblem(problemTSConfiguration);
    
    fprintf('Computing Time-Stepping\n')
    tic;
    UTS = solveTSProblem(problemTS);
    timeTS(refinement) = toc;
    fprintf('Getting the solution Solution\n')
    solutionTS = getTSSolution(problemTS,UTS, resolution);
    
    % Compute the errors
    fprintf('Computing the L2 errors\n')
    errorGalerkin(refinement) = sqrt(mean( (solutionGalerkin-solutionAnalytical).^2, 'all'));
    errorCGLyap(refinement) = sqrt(mean( (solutionCGLyap-solutionAnalytical).^2, 'all'));
    errorCGOpt(refinement) = sqrt(mean( (solutionCGOpt-solutionAnalytical).^2, 'all'));
%     errorBackslash(refinement) = sqrt(mean( (solutionBackslash-solutionAnalytical).^2, 'all'));
    errorTS(refinement) = sqrt(mean( (solutionTS-solutionAnalytical).^2, 'all'));
    
    
end

totalTime = toc(totalTime)

% errorTS
% 
% figure
% subplot(2,2,1)
% surf(squeeze(solutionTS(32,:,:,end) - solutionAnalytical(32,:,:,end)))
% subplot(2,2,2)
% surf(squeeze(solutionTS(:,32,:,end) - solutionAnalytical(:,32,:,end)))
% subplot(2,2,3)
% surf(squeeze(solutionTS(:,:,32,end) - solutionAnalytical(:,:,32,end)))
% subplot(2,2,4)
% surf(squeeze(solutionTS(:,:,32,1) - solutionAnalytical(:,:,32,1)))


% figure
% subplot(1,2,1)
% surf(squeeze(solutionAnalytical(:,:,32,end)))
% subplot(1,2,2)
% surf(squeeze(solutionTS(:,:,32,end)))







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



