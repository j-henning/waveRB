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


sol = zeros(2^resolution.x + 1, 2^resolution.y + 1, ...
    2^resolution.z + 1, 2^resolution.t + 1);




for i=1:size(u_full,1) % every space node in x
    spline_x = splines.space(i,:);
    spline_x_2 = splines.space2(i,:);
    space_x_indices = min(find(spline_x, 1, 'first'), find(spline_x_2, 1, 'first')):...
        max(find(spline_x, 1, 'last'), find(spline_x_2, 1, 'last'));
    spline_x = spline_x(space_x_indices);
    spline_x_2 = spline_x_2(space_x_indices);
    
    for j=1:size(u_full, 2) % every space node in y
        
        spline_y = splines.space(j,:);
        spline_y_2 = splines.space2(j,:);
        space_y_indices = min(find(spline_y, 1, 'first'), find(spline_y_2, 1, 'first')):...
            max(find(spline_y, 1, 'last'), find(spline_y_2, 1, 'last'));
        spline_y = spline_y(space_y_indices);
        spline_y_2 = spline_y_2(space_y_indices);
        
        for k=1:size(u_full, 3) % every space node in z
            spline_z = splines.space(k,:);
            spline_z_2 = splines.space2(k,:);
            space_z_indices = min(find(spline_z, 1, 'first'), find(spline_z_2, 1, 'first')):...
                max(find(spline_z, 1, 'last'), find(spline_z_2, 1, 'last'));
            spline_z = spline_z(space_z_indices);
            spline_z_2 = spline_z_2(space_z_indices);
            
            spl_space = kron(spline_z, kron(spline_y, spline_x));
            spl_space_mixed = kron(spline_z, kron(spline_y, spline_x_2))...
                + kron(spline_z, kron(spline_y_2, spline_x)) ...
                + kron(spline_z_2, kron(spline_y, spline_x));
            
            for l=1:size(u_full,4) % every time node
                
                if abs(u_full(i,j,k,l)) < eps
                    continue;
                end
                
                spline_t = splines.time(l,:);
                spline_t_2 = splines.time2(l,:);
                time_indices = min(find(spline_t, 1, 'first'), find(spline_t_2, 1, 'first')):...
                    max(find(spline_t, 1, 'last'), find(spline_t_2, 1, 'last'));
                spline_t = spline_t(time_indices);
                spline_t_2 = spline_t_2(time_indices);
                
                % Build the whole spline
                %                     spline = kron(spline_t_2, kron(spline_z, kron(spline_y, spline_x))) ...
                %                         - kron(spline_t, kron(spline_z, kron(spline_y, spline_x_2))) ...
                %                         - kron(spline_t, kron(spline_z, kron(spline_y_2, spline_x))) ...
                %                         - kron(spline_t, kron(spline_z_2, kron(spline_y, spline_x)));
                spline = kron(spline_t_2, spl_space) - p.mu * kron(spline_t, spl_space_mixed);
                
                spline = reshape(spline, ...
                    [length(space_x_indices), length(space_y_indices), ...
                    length(space_z_indices), length(time_indices)]);
                sol(space_x_indices, space_y_indices, space_z_indices, time_indices) = ...
                    sol(space_x_indices, space_y_indices, space_z_indices, time_indices) + ...
                    u_full(i,j,k,l) * spline;
            end
        end
    end
end

end

