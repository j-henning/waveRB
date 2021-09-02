% Creates and 'fits' a hp Reduced Basis Model for a problem defined by
% problemConfiguration and a parameter mu in [muMin, muMax] as well as an
% hierarchical error estimator (which needs a theta value to be calculated
% in a later step to properly work)
% Input:
% problemConfiguration - Definition of the problem to solve
% resolution           - Evaluation resolution
% muMin                - Left boundary of the parameter interval
% muMax                - Right boundary of the parameter interval
% Nh                   - Maximum number of snapshots for the h refinement
% Np                   - Maximum number of snapshots for the p refinement
% hTolerance           - Tolerance for the h refinement
% pTolerance           - Tolerance for the p refinement
% height               - Maximum height of the tree
% NhEst                - Nh value for the error estimator
% NpEst                - Np value for the error estimator
% hToleranceEst        - hTolerance value for the error estimator
% pToleranceEst        - pTolerance value for the error estimator
% heightEst            - height parameter for the estimator
% Xi                   - Initial paramter set
% solver               - Name of the solver
% maxIt                - Maximum number of iterations
% tolerance1           - First tolerance of the solver
% tolerance2           - Second tolerance of the solver
% Output:
% tree                 - Reduced Basis Model
% treeEstimator        - Reduced Basis Model for the error estimation
% pOne                 - Problem that can be copied
% splines              - Computed spline data
function [tree, treeEstimator, pOne, splines] = rbOffline(problemConfiguration, ...
    resolution, muMin, muMax, Nh, Np, hTolerance, pTolerance, height, NhEst, NpEst, ...
    hToleranceEst, pToleranceEst, heightEst, Xi, solver, maxIt, tolerance1, tolerance2)

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
    U(:,i) = solveProblem(problem, solver, maxIt, tolerance1, tolerance2);
    sol{i} = getSolution(problem, U(:,i), resolution, splines);
end
fprintf(' Done!\n')
fprintf("Choosing %f as inital mu\n", Xi(1))


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

% Enrich the spaces of the leave nodes
fprintf('Enriching leaf nodes ...')
[tree] = enrich(tree, Np, pOne, splines, resolution, solver, maxIt, tolerance1, tolerance2);
fprintf(' Done!\n')


tree = pGreedy(tree, Np, pTolerance,  pOne, splines, resolution);

figure
subplot(1,2,1)
plotIntervals(tree)
title('RB model')


fprintf('Error estimator ...\n')
fprintf('Computing h refinement for the error estimator...\n')
[treeEstimator] = enrich(treeEstimator, NpEst, pOne, splines, resolution, solver, maxIt, tolerance1, tolerance2);
[treeEstimator] = pGreedy(treeEstimator, Np*3, pToleranceEst,pOne, splines, resolution);
fprintf(' Done!\n')

subplot(1,2,2)
plotIntervals(treeEstimator)
title('Estimator Estimator')
end