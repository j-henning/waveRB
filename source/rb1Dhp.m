% hp RB implementation based on the paper of Eftang, Patera and Ronquist


clear
close all
clc

format short

addpath('../splines');
addpath('../splines/Utilities');
addpath('../solvers');

%% Specify all needed input data
problemConfiguration.bSplineOrder_time = 3;
problemConfiguration.bSplineOrder_space = 3;

problemConfiguration.offset_time_ansatz = [0, 2];
problemConfiguration.offset_time_test = [0, 2];
problemConfiguration.offset_space_ansatz = [1, 1];
problemConfiguration.offset_space_test = [1, 1];


% Example
problemConfiguration.d = 1; % Dimension
problemConfiguration.f_time = {@(t) 0*t};
problemConfiguration.f_space = {@(x) 0*x};
problemConfiguration.u_0_x = @(x) 100 * (x > 0.25 && x < 0.75) ...
    .* ((x-0.25).^2 .* (x -0.75).^2);
problemConfiguration.u_1_x = @(x) 0*x;
problemConfiguration.inital_conditions = true;
problemConfiguration.has_analytical_solution = false;


problemConfiguration.refinementLevel_time = 5;
problemConfiguration.refinementLevel_space = 5;

resolution.x = 7;
resolution.t = 7;


problemConfiguration.mu = 1;

maxIt = 100;
% tolerance = 1e-6;
% toleranceEst = 1e-6; % Error estimator

muMin = .1;
muMax = 2;

Nh = 5; % Number of snapshots during the h refinement
Np = 10; % TODO: Find a good name and value
hTolerance = 1e-1; % Tolerance for the h refinement

NhEst = 10;
NpEst = 15;
hToleranceEst = 1e-2;

Xi = rand(1,Nh) * (muMax - muMin) + muMin;
Xi = sort(Xi, 'ascend');
XiTest = rand(1,200) * (muMax - muMin) + muMin;
XiTest = sort(XiTest, 'ascend');





%% Offline phase

% Create one problem which can be copied and reused
problemConfiguration.mu = 1;
pOne = create1DWaveProblem(problemConfiguration);



U = zeros(2^problemConfiguration.refinementLevel_time * ...
    2^problemConfiguration.refinementLevel_space, length(Xi));
sol = cell(length(Xi),1);


% Precompute the needed splines
[splines.splines_time, splines.splines_2_time, splines.splines_space, ...
    splines.splines_2_space] = precompute1DSplines(pOne,  resolution);

fprintf("Calculating the solutions for mu = %e, ..., %e", Xi(1), Xi(end));
parfor i=1:length(Xi)
    mu = Xi(i);
    
    pC = problemConfiguration;
    pC.mu = mu;
    
    
    problem = pOne;
    problem.mu = mu;
    problem.A_space = mu * problem.A_space;
    problem.Q_space = mu^2 * problem.Q_space;
    
    
    U(:,i) = solveProblem(problem);
    sol{i} = get1DsolutionSmart(problem, U(:,i), resolution, ...
        splines.splines_time, splines.splines_2_time, ...
        splines.splines_space, splines.splines_2_space);
    
    
end
fprintf(' Done!\n')



S = Xi(1);
Y = U(:,1) / norm(U(:,1));
% fprintf("N = %d (%d)\n", 1, maxIt)
fprintf("Choosing %f as inital mu\n", S(1))

% Q_tM_h = kron(pOne.Q_time, pOne.M_space);

%% Greedy algorithm


height = 3;
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
tree = applyHGreedy(tree, hTolerance, Nh, pOne, resolution, splines);
fprintf(' Done!\n')


% Refine the estimator tree
fprintf('Computing h refinement for the error estimator...')
[treeEstimator] = resizeTree(tree, 5);
[treeEstimator] = applyHGreedy(treeEstimator, hToleranceEst, NhEst, pOne, resolution, splines);


 maxErrors = [];

% Enrich the spaces of the leave nodes
fprintf('Enriching leaf nodes ...')
[tree] = enrich(tree, Np, pOne, splines, resolution);
fprintf(' Done!\n')



tree = pGreedy(tree, Np, 1e-3,  pOne, splines, resolution);

figure
subplot(1,2,1)
plotIntervals(tree)
title('RB model')


fprintf('Error estimator ...')

fprintf('Computing h refinement for the error estimator...')


[treeEstimator] = enrich(treeEstimator, Np*3, pOne, splines, resolution);
% 
% [treeEstimator] = enrich(tree, Np*3, pOne, splines, resolution);


[treeEstimator] = pGreedy(treeEstimator, Np*3, 1e-7,pOne, splines, resolution);
fprintf(' Done!\n')


subplot(1,2,2)
plotIntervals(treeEstimator)
title('Estimator Estimator')




%% Calculation of the theta value for the online error estimator
fprintf('Calculating theta ...')
% minProb = @(x) minProblem(x, pOne, Y_N, Y_M,resolution,splines_time, splines_2_time, splines_space, splines_2_space);

% Split the interval from muMin to muMax into several subintervals and find
% for each one the best theta



muInt = linspace(muMin, muMax, 1000);
thetaInt = linspace(0, 1, 200);
theta = 0;


for i=1:length(muInt)
    minVal = Inf;
    
        
        problem = pOne;
        problem.mu = muInt(i);
        problem.A_space = muInt(i) * problem.A_space;
        problem.Q_space = muInt(i)^2 * problem.Q_space;
        
        u_N_rec = getRBhpSolutionVector(tree, muInt(i), problem);
        u_M_rec = getRBhpSolutionVector(treeEstimator, muInt(i), problem);
        differenceVectors(i) = norm(u_N_rec - u_M_rec);
end


[val, indices] = sort(differenceVectors, 'ascend');

indices = indices(1:ceil(1 * length(muInt))); % Only use the biggest 5 percent

% figure
% plot(muInt, differenceVectors), hold on
% plot(muInt(indices), differenceVectors(indices), 'ro')


% muInt = muInt(indices);

for i=1:length(muInt)
    minVal = Inf;
    
        
        problem = pOne;
        problem.mu = muInt(i);
        problem.A_space = muInt(i) * problem.A_space;
        problem.Q_space = muInt(i)^2 * problem.Q_space;
        
        u_N_rec = getRBhpSolutionVector(tree, muInt(i), problem);
        u_M_rec = getRBhpSolutionVector(treeEstimator, muInt(i), problem);
        
        
        sol_rec_N = get1DsolutionSmart(problem, u_N_rec,  resolution, ...
            splines.splines_time, splines.splines_2_time, ...
            splines.splines_space, splines.splines_2_space);
        
        sol_rec_M = get1DsolutionSmart(problem, u_M_rec,  resolution, ...
            splines.splines_time, splines.splines_2_time, ...
            splines.splines_space, splines.splines_2_space);



U = solveProblem(problem);
sol = get1DsolutionSmart(problem, U,  resolution,... 
    splines.splines_time, splines.splines_2_time, splines.splines_space, splines.splines_2_space);


fmu = sqrt(mean( (sol-sol_rec_M).^2, 'all'));
gmu = sqrt(mean( (sol-sol_rec_N).^2, 'all'));

% differenceVectors(i) = norm(u_N_rec - u_M_rec);
% differenceSolutions(i) = sqrt(mean( (sol_rec_M-sol_rec_N).^2, 'all'));




val = abs(fmu - thetaInt * gmu);
[minVal, minIndex] = min(val);


thetaLoc(i) = thetaInt(minIndex);
    if thetaLoc(i) > 0.9
        fprintf('Theta local = %f, mu = %f\n', thetaLoc(i), muInt(i));
    end

        if thetaLoc(i) > theta
            theta = thetaLoc(i);
        end

end

% figure
% plot(muInt, thetaLoc)
% xlabel('mu')
% ylabel('theata')
% drawnow
% 


fprintf(' Done!\n')

% for i=1:length(muInt) - 1
%     [~,fval,exitflag,output] = fmincon(minProb,[(muInt(i) + muInt(i+1)) / 2, 0],[],[],[],[],[muInt(i)-eps, 0-eps], [muInt(i+1)+eps, 2+eps]);
%     x = output.bestfeasible.x;
%     theta = max(theta, x(2));
%
%     if theta > 1
%         warning('Theta is larger than 1! Stopping theta calculation!')
%         break;
%     end
% end

% differenceSolutions = 1./differenceSolutions;
% differenceSolutions = differenceSolutions ./ max(differenceSolutions);
% 
% differenceVectors = 1./differenceVectors;
% differenceVectors = differenceVectors ./ max(differenceVectors);
% 
% thetaLoc = thetaLoc ./ max(thetaLoc);

% figure
% 
% [muInt, perm] = sort(muInt);
% thetaLoc = thetaLoc(perm);
% 

figure
subplot(2,1,1)
plot(muInt, thetaLoc, 'r'), grid on 
xlabel('\mu')
title('Theta')
% semilogy(muInt, differenceSolutions)
subplot(2,1,2)
semilogy(muInt, 1./differenceVectors), grid on, hold on
indicesTheta = find(thetaLoc > 0);
plot(muInt(indicesTheta), 1./differenceVectors(indicesTheta), 'ro')
xlabel('\mu')
title('1 / norm(U_{est} - U)_{Euklidean}')
legend('1 / norm(U_{est} - U)_{Euklidean}', '\theta > 0')
drawnow
% title('Vector difference')
% 
% figure
% plotError(tree, 'b')
% hold on
% title('Max error')
% 
% plotError(tree, 'b')
% plotError(treeEstimator, 'g')
% title('Max error estimator')
% grid on
% drawnow
% 


fprintf('Theta = %f\n', theta)

% R = corrcoef([thetaLoc; differenceSolutions; differenceVectors]')


%% Online phase: Check error and hierarchical error estimator for different
% parameters
fprintf('Starting online phase\n')
solTest = cell(length(XiTest),1);



fprintf("Calculating the test solutions for mu = %e, ..., %e", XiTest(1), XiTest(end));
parfor i=1:length(XiTest)
    mu = XiTest(i);
    
    pC = problemConfiguration;
    pC.mu = mu;
    
    
    problem = pOne;
    problem.mu = mu;
    problem.A_space = mu * problem.A_space;
    problem.Q_space = mu^2 * problem.Q_space;
    
    
    Utest = solveProblem(problem);
    solTest{i} = get1DsolutionSmart(problem, Utest,  resolution, ...
        splines.splines_time, splines.splines_2_time, ...
        splines.splines_space, splines.splines_2_space);
    
    
end
fprintf(' Done!\n')



fprintf('Testing all solutions')
parfor j=1:length(XiTest)
    
    mu = XiTest(j);
    
    p = pOne;
    p.A_space = mu * p.A_space;
    p.Q_space = mu^2 * p.Q_space;
    p.mu = mu;
    
    u_N_rec = getRBhpSolutionVector(tree, mu, p);
    u_M_rec = getRBhpSolutionVector(treeEstimator, mu, p);

    
    sol_rec_N = get1DsolutionSmart(p, u_N_rec,  resolution, ...
        splines.splines_time, splines.splines_2_time, ...
        splines.splines_space, splines.splines_2_space);
    
    sol_rec_M = get1DsolutionSmart(p, u_M_rec,  resolution, ...
        splines.splines_time, splines.splines_2_time, ...
        splines.splines_space, splines.splines_2_space);

    errorTest(j) = sqrt(mean( (solTest{j}-sol_rec_N).^2, 'all'));
    errorEstTest(j) = sqrt(mean( (sol_rec_M-sol_rec_N).^2, 'all'));
end
fprintf(' Done!\n')
% fprintf('Desired tolerance : %e\n', tolerance)
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
semilogy(XiTest(find(errorEstTest./errorTest < 1)), ...
    errorEstTest(find(errorEstTest./errorTest < 1)) ...
    ./errorTest(find(errorEstTest./errorTest < 1)), 'ro')
title('Error Estimator / Error')
legend('Estimator / Error', 'Values smaller than 1')
grid on








%% Help functions


% Create a function for the minimization problem
function val = minProblem(x, pOne, Y_N, Y_M,resolution,splines_time, splines_2_time, splines_space, splines_2_space)
mu = x(1);
q = x(2);

% fprintf('Testing mu = %f, q = %f', mu, q);

p = pOne;
p.mu = mu;
p.A_space = mu * p.A_space;
p.Q_space = mu^2 * p.Q_space;


% Compute the reduced system

B_N = Y_N' * (kron(p.Q_time, p.M_space) ...
    + kron(p.D_time, p.A_space') ...
    + kron(p.D_time', p.A_space) ...
    + kron(p.M_time, p.Q_space)) * Y_N;
f_N = Y_N' * p.rhs(:);
u_N = B_N \ f_N;
u_N_rec = Y_N * u_N;

B_M = Y_M' * (kron(p.Q_time, p.M_space) ...
    + kron(p.D_time, p.A_space') ...
    + kron(p.D_time', p.A_space) ...
    + kron(p.M_time, p.Q_space)) * Y_M;
f_M = Y_M' * p.rhs(:);
u_M = B_M \ f_M;
u_M_rec = Y_M * u_M;


sol_rec_N = get1DsolutionSmart(p, u_N_rec,  resolution,splines_time, splines_2_time, splines_space, splines_2_space);
sol_rec_M = get1DsolutionSmart(p, u_M_rec,  resolution,splines_time, splines_2_time, splines_space, splines_2_space);


% Compute the reference solutions

problem = pOne;
problem.mu = mu;
problem.A_space = mu * problem.A_space;
problem.Q_space = mu^2 * problem.Q_space;


U = solveProblem(problem);
sol = get1DsolutionSmart(problem, U,  resolution,splines_time, splines_2_time, splines_space, splines_2_space);


fmu = sqrt(mean( (sol-sol_rec_M).^2, 'all'));
gmu = sqrt(mean( (sol-sol_rec_N).^2, 'all'));


val = fmu - q * gmu;
% fprintf(' - val = %e\n', val)
end


