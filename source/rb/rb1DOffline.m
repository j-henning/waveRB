% Nh = 5; % Number of snapshots during the h refinement
% Np = 10; % Number of snapshots during the p refinement (Np >= Nh)
% hTolerance = 1e-1; % Tolerance for the h refinement
% 
% NhEst = 10; % Number of snapshots during the h refinement for the Estimator
% NpEst = 15; % Number of snapshots during the p refinement (Np >= Nh) for the Estimator
% hToleranceEst = 1e-2; % Tolerance for the h refinement for the Estimator
% 
% Xi = rand(1,Nh) * (muMax - muMin) + muMin;
% Xi = sort(Xi, 'ascend');
% XiTest = rand(1,200) * (muMax - muMin) + muMin;
% XiTest = sort(XiTest, 'ascend');

% TODO: Remove the 1D name
function [tree, treeEstimator] = rb1DOffline(problemConfiguration, ...
    resolution, muMin, muMax, Nh, Np, hTolerance, height, NhEst, NpEst, ...
    hToleranceEst, heightEst, Xi, solver, maxIt, tolerance1, tolerance2)

% Create one problem which can be scaled 
pOne = createProblem(problemConfiguration);
if pOne.mu ~= 1 % mu does not need to be one, but we scale it to be one for completness
   pOne = changeWaveSpeed(pOne, 1);
end

splines = computeSplines(pOne, resolution); 


% Compute the inital solutions
U = zeros(2^problemConfiguration.refinementLevel_time * ...
    (2^problemConfiguration.dimension)^problemConfiguration.refinementLevel_space, length(Xi));
sol = cell(length(Xi),1);

fprintf("Calculating the solutions for mu = %e, ..., %e", Xi(1), Xi(end));
parfor i=1:length(Xi) % TODO: Par for
    mu = Xi(i);
    
    problem = changeWaveSpeed(pOne, mu)
    
    
    U(:,i) = solveProblem(problem, solver, maxIt, tolerance1, tolerance2); % TODO: Make this variable
    sol{i} = getSolution(problem, U(:,i), resolution, splines);
    
    
end
fprintf(' Done!\n')



S = Xi(1);
Y = U(:,1) / norm(U(:,1));
% fprintf("N = %d (%d)\n", 1, maxIt)
fprintf("Choosing %f as inital mu\n", S(1))

% Q_tM_h = kron(pOne.Q_time, pOne.M_space);

%% Greedy algorithm



tree = cell(2^height-1,1);
for i = 1:2^height - 1
    tree{i} = struct('muMin', NaN, 'muMax', NaN, 'U', NaN, 'Y', NaN, ...
        'Xi', NaN,'S', NaN, 'solutions', NaN, 'maxError', Inf, 'center', NaN);
end

% Set the root node
tree{1}.muMin = muMin;
tree{1}.muMax = muMax;
tree{1}.U = U;
tree{1}.Xi = Xi';
tree{1}.S = Xi(1);
tree{1}.Y = U(:,1) / norm(U(:,1));
tree{1}.solutions = sol;
tree{1}.maxError = Inf;

fprintf('Computing h refinement ...')
tree = applyHGreedy(tree, hTolerance, Nh, pOne, resolution, splines, solver, maxIt, tolerance1, tolerance2);
fprintf(' Done!\n')


% Refine the estimator tree
fprintf('Computing h refinement for the error estimator...')
[treeEstimator] = resizeTree(tree, heightEst);
[treeEstimator] = applyHGreedy(treeEstimator, hToleranceEst, NhEst, pOne, resolution, splines, solver, maxIt, tolerance1, tolerance2);


 maxErrors = [];

% Enrich the spaces of the leave nodes
fprintf('Enriching leaf nodes ...')
[tree] = enrich(tree, Np, pOne, splines, resolution, solver, maxIt, tolerance1, tolerance2);
fprintf(' Done!\n')



tree = pGreedy(tree, Np, 1e-3,  pOne, splines, resolution);

figure
subplot(1,2,1)
plotIntervals(tree)
title('RB model')


fprintf('Error estimator ...')

fprintf('Computing h refinement for the error estimator...')


[treeEstimator] = enrich(treeEstimator, Np*3, pOne, splines, resolution, solver, maxIt, tolerance1, tolerance2);
% 
% [treeEstimator] = enrich(tree, Np*3, pOne, splines, resolution);


[treeEstimator] = pGreedy(treeEstimator, Np*3, 1e-7,pOne, splines, resolution);
fprintf(' Done!\n')


subplot(1,2,2)
plotIntervals(treeEstimator)
title('Estimator Estimator')





end

