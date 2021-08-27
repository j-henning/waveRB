% Evaluates the solution of a problem on with the desired resolution at all
% query points in space and time on [0,1]^2 x [0,1]
function [sol] = get2Dsolution(p, U,  resolution, splines)
if resolution.x ~= resolution.y 
    error('x, y and z resolution should be the same for performance reasons!')
end

% Todo: Check transpose
u = reshape(full(U), round((numel(U) / p.dim_time)^(1/2)), ...
    round((numel(U) / p.dim_time)^(1/2)), ...
    p.dim_time);

% Preparing boundary for plot
u_full = zeros(p.nsol_space, p.nsol_space, p.nsol_time);


u_full(1+p.offset_space_ansatz(1):p.nsol_space-p.offset_space_ansatz(2),...
    1+p.offset_space_ansatz(1):p.nsol_space-p.offset_space_ansatz(2), ...
    1+p.offset_time_ansatz(1):p.nsol_time-p.offset_time_ansatz(2)) = u;


timeDim = size(splines.time,1);
spaceDim = size(splines.space,1);

% Define some shortcuts
s = splines.space';
s2 = splines.space2';
t = splines.time';
t2 = splines.time2';

U = reshape(u_full, [spaceDim * timeDim, spaceDim]);


sol = kron(s,s) * U * t2' ...
    - p.mu * ((kron(s,s2) + kron(s2,s)) * U * t');

sol = reshape(sol, [size(splines.space,2), size(splines.space,2), size(splines.time,2)]);

end

