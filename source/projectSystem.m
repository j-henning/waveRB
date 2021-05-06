function [A_Nq, f_Nq] = projectSystem(A_hq, f_hq, X, Y)

% Compute A_N
A_Nq = cell(numel(A_hq), 1);
for q = 1:numel(A_hq)
    A_Nq{q} = X' * A_hq{q} * Y; 
end

% Compute f_N
f_Nq = cell(numel(f_hq), 1);
for q = 1:numel(f_hq)
    f_Nq{q} = Y' * f_hq{q};
end

end