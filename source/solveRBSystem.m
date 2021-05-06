function u_N = solveRBSystem(A_Nq, f_Nq,theta_Aq, theta_fq, mu)

% Compute A(mu)
A_N = zeros(size(A_Nq{1}));
for q = 1:nummel(A_Nq)
    A_N = A_N + theta_Aq{q}(mu) * A_Nq{q};
end


% Compute f(mu)
f_N = zeros(size(f_Nq{1}));
for q = 1:numel(f_Nq)
   f_N = f_N + theta_fq{q}(mu) * f_Nq{q}; 
end

u_N = A_N \ f_N;
end