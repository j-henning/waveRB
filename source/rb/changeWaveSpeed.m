% Copies an existing problem with a wave speed mu_old and alters all
% matrices (and the wave speed parameter) so they represent the new mu
% Input:
% problem - old problem with an arbitrary wave speed
% mu      - new wave speed
% Output:
% problem - copy of the old problem with the new wave speed
function [problem] = changeWaveSpeed(problem, mu)
problem.Q_space = mu^2 / problem.mu^2 *  problem.Q_space;
problem.A_space = mu / problem.mu * problem.A_space;
problem.mu = mu;
end