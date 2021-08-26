function theta = thetaCalculation(problemConfiguration, muMin, muMax, Nmu, Ntheta, percentage, tree, treeEstimator, resolution,  solver, maxIt, tolerance1, tolerance2)
pOne = createProblem(problemConfiguration);
splines = computeSplines(pOne, resolution);
fprintf('Calculating theta ...')
% minProb = @(x) minProblem(x, pOne, Y_N, Y_M,resolution,splines_time, splines_2_time, splines_space, splines_2_space);

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

indices = indices(1:ceil(percentage * length(muInt))); % Only use the biggest 5 percent

% figure
% plot(muInt, differenceVectors), hold on
% plot(muInt(indices), differenceVectors(indices), 'ro')


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
    
    % differenceVectors(i) = norm(u_N_rec - u_M_rec);
    % differenceSolutions(i) = sqrt(mean( (sol_rec_M-sol_rec_N).^2, 'all'));
    
    
    
    
    val = abs(fmu - thetaInt * gmu);
    [~, minIndex] = min(val);
    
    
    thetaLoc(i) = thetaInt(minIndex);
    
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

% figure
% subplot(2,1,1)
% plot(muInt, thetaLoc, 'r'), grid on
% xlabel('\mu')
% title('Theta')
% % semilogy(muInt, differenceSolutions)
% subplot(2,1,2)
% semilogy(muInt, 1./differenceVectors), grid on, hold on
% indicesTheta = find(thetaLoc > 0);
% plot(muInt(indicesTheta), 1./differenceVectors(indicesTheta), 'ro')
% xlabel('\mu')
% title('1 / norm(U_{est} - U)_{Euklidean}')
% legend('1 / norm(U_{est} - U)_{Euklidean}', '\theta > 0')
% drawnow
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

end