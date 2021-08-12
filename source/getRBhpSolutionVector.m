% Return the coefficient vector for a given binary tree and a parameter mu
% and a problem p
function U = getRBhpSolutionVector(tree, mu, p)
    index = findInterval(tree, mu);
    
    B = tree{index}.Y' * (kron(p.Q_time, p.M_space) ...
        + kron(p.D_time, p.A_space') ...
        + kron(p.D_time', p.A_space) ...
        + kron(p.M_time, p.Q_space)) * tree{index}.Y;
    f = tree{index}.Y' * p.rhs(:);
    
    U = B \ f;
    U = tree{index}.Y * U;
end