% Evaluates the solution of a problem on with the desired resolution at all
% query points in space and time on [0,1]^1 x [0,1]
% Input:
% p          - Problem
% U          - Solution vector
% resolution - Resolution
% splines    - Splines calculated by computeSplines
% Output:
% sol        - Solution
function [sol] = get1Dsolution(p, U,  resolution, splines)

u = reshape(full(U), numel(U) / p.dim_time, p.dim_time);

% Preparing boundary for plot
u_full = zeros(p.nsol_space, p.nsol_time);

u_full(1+p.offset_space_ansatz(1):p.nsol_space-p.offset_space_ansatz(2),...
    1+p.offset_time_ansatz(1):p.nsol_time-p.offset_time_ansatz(2)) = u;

sol = zeros(2^resolution.x + 1, 2^resolution.t + 1);

sol = splines.space' * u_full * splines.time2 ...
    - p.mu * splines.space2' * u_full * splines.time;

end

