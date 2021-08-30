% Approximates the theta value needed for the hierarchical error estimator
% by every grid point of the grid muMin:...:muMax x 0:...1 of dimension
% Nmu x Ntheta
% Input:
% problemConfiguration - Problem definition
% muMin                - Left boundary of the parameter interval
% muMax                - Right boundary of the parameter interval
% Nmu                  - Number of query points in the parameter interval
% Ntheta               - Number of query points in the theta interval
% percentage           - Percentage of points that get exactly calculated
% tree                 - Reduced Basis Model
% treeEstimator        - Reduced Basis Model for the error estimation
% resolution           - Evaluation resolution
% solver               - Name of the solver
% maxIt                - Maximum number of iterations
% tolerance1           - First tolerance for the solver
% tolerance2           - Second tolerance for the solver
% Ouput:
% theta                - Estimation of the theta value
function theta = thetaCalculation(problemConfiguration, muMin, muMax, Nmu, Ntheta, percentage, tree, treeEstimator, resolution,  solver, maxIt, tolerance1, tolerance2)
pOne = createProblem(problemConfiguration);
splines = computeSplines(pOne, resolution);
fprintf('Calculating theta ...')

% Split the interval from muMin to muMax into several subintervals and find
% for each one the best theta
muInt = linspace(muMin, muMax, Nmu);
thetaInt = linspace(0, 1, Ntheta);
theta = 0;

differenceVectors = zeros(size(muInt));
for i=1:length(muInt)
    
    problem = changeWaveSpeed(pOne, muInt(i));
    
    u_N_rec = getRBhpSolutionVector(tree, muInt(i), problem);
    u_M_rec = getRBhpSolutionVector(treeEstimator, muInt(i), problem);
    differenceVectors(i) = norm(u_N_rec - u_M_rec);
end

[~, indices] = sort(differenceVectors, 'ascend');
indices = indices(1:ceil(percentage * length(muInt)));
muInt = muInt(indices);
thetaLoc = zeros(size(muInt));

for i=1:length(muInt)
    problem = changeWaveSpeed(pOne, muInt(i));
    
    u_N_rec = getRBhpSolutionVector(tree, muInt(i), problem);
    u_M_rec = getRBhpSolutionVector(treeEstimator, muInt(i), problem);
    
    sol_rec_N = getSolution(problem, u_N_rec,  resolution, ...
        splines);
    sol_rec_M = getSolution(problem, u_M_rec,  resolution, ...
        splines);
    
    U = solveProblem(problem, solver, maxIt, tolerance1, tolerance2);
    sol = getSolution(problem, U,  resolution);
    
    fmu = sqrt(mean( (sol-sol_rec_M).^2, 'all'));
    gmu = sqrt(mean( (sol-sol_rec_N).^2, 'all'));
    
    val = abs(fmu - thetaInt * gmu);
    [~, minIndex] = min(val);
    
    thetaLoc(i) = thetaInt(minIndex);
    if thetaLoc(i) > theta
        theta = thetaLoc(i);
    end
    
end
fprintf(' Done!\n')
fprintf('Theta = %f\n', theta)
end