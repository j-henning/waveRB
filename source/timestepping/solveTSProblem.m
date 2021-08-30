% Solves the time-stepping problem. For every timestep, a linear system is
% solved with a CG method using MATLABs PCG implementation
% Input:
% problemTS - Time-stepping problem
% maxIt     - Maximum number of iterations per time-step (default = 100)
% tolerance - Tolerance for the CG method (default = 1e-7)
% Output:
% U         - Solution vector
function [U] = solveTSProblem(problemTS, maxIt, tolerance)
if nargin < 2
    maxIt = 100;
end

if nargin < 3
    tolerance = 1e-7;
end

U = zeros(size(problemTS.M,2), problemTS.K);
% Save U0
[U(:,1), ~] = pcg(problemTS.M,problemTS.U0,tolerance, maxIt);

% Calculate the first step:
[temp1, ~] = pcg(problemTS.M,problemTS.U1, tolerance, maxIt);
[temp2, ~] = pcg(problemTS.M,(problemTS.Fvals(:,1) ...
    - problemTS.A*U(:,1)), tolerance, maxIt);

U(:,2) = U(:,1) + problemTS.tau* temp1 ...
    + (problemTS.tau*problemTS.tau/2) * temp2;

ts4 = problemTS.tau*problemTS.tau/4;

% Calculate the other steps:
for k=3:problemTS.K
    Fpp = problemTS.Fvals(:,k-2); % penultimate time step
    Fp = problemTS.Fvals(:,k-1); % last time step
    F = problemTS.Fvals(:,k); % current time step
    
    % Crank-Nicolson
    [U(:, k), ~] = pcg(problemTS.M + ts4 * problemTS.A, ( ts4* (F + 2*Fp+ Fpp) ...
        - ts4 * problemTS.A * (2*U(:, k-1) + U(:, k-2)) ...
        + problemTS.M * (2*U(:,k-1) - U(:, k-2))), tolerance, maxIt, [], [], U(:, k-1));
end
end