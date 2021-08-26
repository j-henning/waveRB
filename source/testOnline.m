function [] = testOnline(problemConfiguration, resolution, XiTest)

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

end