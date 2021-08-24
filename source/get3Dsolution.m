% Evaluates the solution of a problem on with the desired resolution at all
% query points in space and time on [0,1]^3 x [0,1]
function [sol] = get3Dsolution(p, U,  resolution, splines)
if resolution.x ~= resolution.y || resolution.x ~= resolution.z
    error('x, y and z resolution should be the same for performance reasons!')
end

% x = linspace(0, 1, 2^resolution.x + 1);
% y = linspace(0, 1, 2^resolution.y + 1);
% z = linspace(0, 1, 2^resolution.z + 1);
% t = linspace(0, 1, 2^resolution.t + 1);

% if smoothing
% x = 0.5 * (x(2:end) + x(1:end-1));
% y = 0.5 * (y(2:end) + y(1:end-1));
% z = 0.5 * (z(2:end) + z(1:end-1));
% t = 0.5 * (t(2:end) + t(1:end-1));
% end


% Todo: Check transpose
u = reshape(full(U), round((numel(U) / p.dim_time)^(1/3)), ...
    round((numel(U) / p.dim_time)^(1/3)), ...
    round((numel(U) / p.dim_time)^(1/3)), ...
    p.dim_time);

% Preparing boundary for plot
u_full = zeros(p.nsol_space, p.nsol_space, p.nsol_space, p.nsol_time);


u_full(1+p.offset_space_ansatz(1):p.nsol_space-p.offset_space_ansatz(2),...
    1+p.offset_space_ansatz(1):p.nsol_space-p.offset_space_ansatz(2), ...
    1+p.offset_space_ansatz(1):p.nsol_space-p.offset_space_ansatz(2), ...
    1+p.offset_time_ansatz(1):p.nsol_time-p.offset_time_ansatz(2)) = u;


timeDim = size(splines.time,1);
spaceDim = size(splines.space,1);

% Define some shortcuts
s = splines.space';
s2 = splines.space2';
t = splines.time';
t2 = splines.time2';

U = reshape(u_full, [spaceDim^2, spaceDim * timeDim]);


sol = kron(s,s) * U * kron(t2, s)' ...
    - p.mu * ( kron(s,s2) * U * kron(t,s)' ...
    + kron(s2,s) * U * kron(t,s)' ...
    + kron(s,s) * U * kron(t,s2)');

sol = reshape(sol, [size(splines.space,2), size(splines.space,2), size(splines.space,2) size(splines.time,2)]);

end

