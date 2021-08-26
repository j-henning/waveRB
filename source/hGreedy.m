function [muNew, muNewIndex, errorMax, Xi, solutions, U] = hGreedy(Xi, solutions, U, Y, N, muMin, muMax, pOne, resolution, splines, solver, maxIt, tolerance1, tolerance2)
% Sample new parameters to get barN parameters in the interval muMin, muMax
newN = N - length(Xi);
XiNew = rand(newN, 1) * (muMax - muMin) + muMin;
snapshotsNew = cell(newN,1);
UNew = zeros(size(U,1), newN);


% Compute the new snapshots
for i = 1:newN
    mu = XiNew(i);
    problem = pOne;
    
    problem.mu = mu;
    problem.A_space = mu * problem.A_space;
    problem.Q_space = mu^2 * problem.Q_space;
    
    
    UNew(:,i) = solveProblem(problem, solver, maxIt, tolerance1, tolerance2);
    snapshotsNew{i} = getSolution(problem, UNew(:,i), ...
        resolution, splines);
end

% Concatinate the parameters and the snapshots
Xi = [Xi; XiNew];
solutions = [solutions; snapshotsNew];
U = horzcat(U, UNew);


% Search the next parameter snapshot
errorMax = 0;
muNewIndex = -1;
for i=1:length(Xi)
    problem = pOne;
    mu = Xi(i);
    problem.mu = mu;
    problem.A_space = mu * problem.A_space;
    problem.Q_space = mu^2 * problem.Q_space;
    B = Y' * (kron(problem.Q_time, problem.M_space) ...
        + kron(problem.D_time, problem.A_space') ...
        + kron(problem.D_time', problem.A_space) ...
        + kron(problem.M_time, problem.Q_space)) * Y;
    f = Y' * problem.rhs(:);
    
    % Compute the reduces solution u_N_rec
    u_N = B \ f;
    
    u_N_rec = Y * u_N;
    sol_rec = getSolution(problem, u_N_rec, ...
        resolution, splines);
    err = sqrt(mean( (solutions{i}-sol_rec).^2, 'all'));
    
    if err > errorMax
        errorMax = err;
        muNewIndex = i;
    end
    
end
muNew = Xi(muNewIndex);



end