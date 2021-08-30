% Creates a 3D problem which is parameter dependant and builds an hp-RB
% space for that problem. Afterwords, the space as well an error estimator
% is tested on a set of test points.

clear
close all
clc

addpath(genpath('../source'))

% Define the parameter range
muMin = .01;
muMax = 4;

% Define the resolution
resolution.x = 4;
resolution.y = 4;
resolution.z = 4;
resolution.t = 4;

Nh = 5; % Number of snapshots during the h refinement
Np = 10; % Number of snapshots during the p refinement (Np >= Nh)
hTolerance = 1e-1; % Tolerance for the h refinement
height = 2; % Tree height

NhEst = 10; % Number of snapshots during the h refinement for the Estimator
NpEst = 15; % Number of snapshots during the p refinement (Np >= Nh) for the Estimator
hToleranceEst = 1e-2; % Tolerance for the h refinement for the Estimator
heightEst = 5; % Tree height for the Estimator

Xi = rand(1,Nh) * (muMax - muMin) + muMin;
Xi = sort(Xi, 'ascend');
XiTest = rand(1,200) * (muMax - muMin) + muMin;
XiTest = sort(XiTest, 'ascend');

solver = 'galerkin';
maxIt = 100;
tolerance1 = 1e-10;
tolerance2 = 1e-2;

mm = 3;
nn = 2;
oo = 1;
% TODO: Make this wave speed dependent
f_time = {@(t) 2+ 0.*t; @(t) (t.^2)*(mm^2 + nn^2 + oo^2)*(pi^2)};
f_space = {@(x) sin(mm*pi*x), @(y) sin(nn*pi*y), @(z) sin(oo*pi*z); ...
    @(x) sin(mm*pi*x), @(y) sin(nn*pi*y), @(z) sin(oo*pi*z)};

u_0 = [];
u_1 = [];
mu = 1;
dimension = 3;
refinement = 3;
splineOrder = 3;


problemConfiguration = defineProblem(dimension, f_time, f_space, u_0, ...
    u_1, mu, refinement, refinement, splineOrder, splineOrder);




%% Offline phase

% Calculate the hp-RB and the error estimator
t = tic;
[tree, treeEstimator] = rbOffline(problemConfiguration, ...
    resolution, muMin, muMax, Nh, Np, hTolerance, height,...
    NhEst, NpEst,hToleranceEst, heightEst, Xi, ...
    solver, maxIt, tolerance1, tolerance2);
toc(t)
% TODO: Return splines and pOne for the theta Calculation

% Calculate the theta value for the error estimator
t = tic;
theta = thetaCalculation(problemConfiguration, muMin, muMax, 1000, 100, 0.05, tree, treeEstimator, resolution, solver, maxIt, tolerance1, tolerance2);
toc(t)


%% Online phase

solTest = cell(length(XiTest),1);



pOne = createProblem(problemConfiguration);
splines = computeSplines(pOne, resolution);
fprintf('Testing all solutions')

timesRB = zeros(length(XiTest),1);
timesNormal = zeros(length(XiTest),1);
parfor j=1:length(XiTest)
    
    mu = XiTest(j);
    
    p = changeWaveSpeed(pOne, mu);
    
    tic;
    u_N_rec = getRBhpSolutionVector(tree, mu, p);
    timesRB(j) = toc;
    u_M_rec = getRBhpSolutionVector(treeEstimator, mu, p);
    
    
    sol_rec_N = getSolution(p, u_N_rec,  resolution, ...
        splines);
    
    sol_rec_M = getSolution(p, u_M_rec,  resolution, ...
        splines);
    
    tic;
    Utest = solveProblem(p, solver, maxIt, tolerance1, tolerance2);
    timesNormal(j) = toc;
    
    solTest{j} = getSolution(p, Utest,  resolution, ...
        splines);
    errorTest(j) = sqrt(mean( (solTest{j}-sol_rec_N).^2, 'all'));
    errorEstTest(j) = sqrt(mean( (sol_rec_M-sol_rec_N).^2, 'all'));
end

fprintf('Time RB:          %f\n', sum(timesRB))
fprintf('Time normal:      %f\n', sum(timesNormal))
fprintf('Speed up facteor: %f\n',sum(timesNormal) / sum(timesRB))
fprintf(' Done!\n')
fprintf('Maximum test error: %e\n', max(errorTest))

errorEstTest = errorEstTest ./ (1 - theta);

% Plotting
figure
subplot(1,2,1)
semilogy(XiTest, errorTest, '*-'), hold on, grid on
semilogy(XiTest, errorEstTest, 'd:')
% semilogy(S, 0.5*min(errorTest)*ones(size(S)), 'k*')
xlabel('\mu')
ylabel('Error')
legend('Error', 'Error estimator', 'Snapshot parameters')
title(['Error over test set (theta = ' num2str(theta) ', factor = ' num2str(1/(1-theta)) ')'])

subplot(1,2,2)
semilogy(XiTest,errorEstTest./errorTest), hold on
semilogy(XiTest(errorEstTest./errorTest < 1), ...
    errorEstTest(errorEstTest./errorTest < 1) ...
    ./errorTest(errorEstTest./errorTest < 1), 'ro')
title('Error Estimator / Error')
legend('Estimator / Error', 'Values smaller than 1')
grid on