% Solves
% (A^t (x) M^s + M^t (x) A^s - N^t (x) (N^s)^T - (N^t)^T (x) N^s) x = b
function [U, errCurrent] = restarting_GMRES(A_time, M_space, M_time, A_space, N_time, ...
    N_space, F, max_iterations, threshold, GMRES_k)
x_sol = zeros(size(F,1)*size(F,2),1); % Initial guess
b = reshape(F, [size(F,1)*size(F,2) 1]);
K = size(F,2);

errCurrent = Inf;
iterations_left = max_iterations;

while errCurrent > threshold
    [x_sol, r] = gmresKron(A_time, M_space, M_time, A_space, ...
        -N_time, N_space', -N_time', N_space, b, x_sol,...
        min(GMRES_k,iterations_left), threshold);
    errCurrent = r(end); 
    iterations_left = iterations_left - GMRES_k;
    if iterations_left <= 0
        break;
    end
end


if errCurrent > threshold
    warning('GMRES did NOT reach the target threshold!')
end


U = reshape(x_sol, [size(M_space,1) K]);
end