function [V] = weakGreedy(N_max, eps_g, Xi, mu_1, pC)

Xi_g = [];
V = [];
N = 1;
delta(1) = eps_g + 1;
mu(1) = mu_1;

while N <= N_max && delta(N) > eps_g
    % Create and solve a new problem for the truth solution
    pC.c = mu(N);
    problem = create1DWaveProblem(pC);
    u_new = solveProblem(problem);
    
    
    % Add Gram-Schmidt if wanted
    if isempty(V)
        z = u_new(:);
    else
        z = u_new(:) - problem.V * problemV' * problem.X * u;
    end
    z = z / sqrt(z' * problem.X * z); 
    
    V = [V z];
    
    Xi_g = [Xi_g mu(N)];
    
    N = N + 1;
    [delta(N), mu(N)] = maxError(Xi, Xi_g, problem);
end
end

% TODO
function [delta, mu] = maxError(Xi, Xi_g, problem)
deltaMax = 0;
mu = inf;
for i = 1:length(Xi) % Test every mu
    if ismember(Xi(i), Xi_q)
        % The parameter is already in our training set
        continue;
    end
    
    % TODO: Fast error estimator
    
    % Get the reduced solution for the current mu
    
    % Iterative over every v in V_h to calculate the error in the dual norm
    dualError = 0;
    for j=1:numel(problem.rhs) % TODO Itera
        currentError = 0;
        g = problem.rhs(j);
        
        % Get the reduced solution u_N(mu)
        
        
        
        b = evaluateLSH(mu, p, index)
        
        dualError = dualError + currentError;
    end
    
    if dualError > deltaMax
        deltaMax = dualError;
        mu = Xi(i);
    end
end

end

function sol = solveProblem(problem)


% Use Galerkin for solving
[uu,ss,vv]=svds(problem.rhs,1);
rhs1=uu(:,1)*sqrt(ss(1,1));
rhs2=vv(:,1)*sqrt(ss(1,1));
tolG=1e-5;
maxIt = 100;
info = 1;
[X1,X2]= ...
    Galerkin3(problem.M_space,2*problem.A_space,problem.Q_space,problem.Q_time,...
    (problem.D_time+problem.D_time')/2,problem.M_time,rhs1,rhs2,maxIt,tolG,info);
U=X1*X2';
U=U(:);

resolution.x = 8;
resolution.t = 8;
sol = get1Dsolution(problem, U, resolution);

end



function b = evaluateLSH(mu, p, index)
% todo
b = 0;
end