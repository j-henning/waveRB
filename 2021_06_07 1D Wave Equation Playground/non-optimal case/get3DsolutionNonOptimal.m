% Evaluates the solution of a problem on with the desired resolution at all
% query points in space and time on [0,1]^3 x [0,1]
function [sol] = get3DsolutionNonOptimal(p, U,  resolution, smoothing)
if resolution.x ~= resolution.y || resolution.x ~= resolution.z
    error('x, y and z resolution should be the same for performance reasons!')
end

if nargin < 4
    smoothing = false;
end

c = 1;
if isfield(p, 'c')
    c = p.c;
end

x = linspace(0, 1, 2^resolution.x + 1);
y = linspace(0, 1, 2^resolution.y + 1);
z = linspace(0, 1, 2^resolution.z + 1);
t = linspace(0, 1, 2^resolution.t + 1);

if smoothing
x = 0.5 * (x(2:end) + x(1:end-1));
y = 0.5 * (y(2:end) + y(1:end-1));
z = 0.5 * (z(2:end) + z(1:end-1));
t = 0.5 * (t(2:end) + t(1:end-1));
end


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


sol = zeros(length(x), length(y), length(z), length(t));

% Precalculate the splines

splines_space = zeros(size(u_full,1), length(x));

for i=1:size(splines_space,1)
     if i < 5 || i > size(splines_space,1) - 5
        for j=1:size(splines_space,2)
            
%             if j > 3 && splines_space(i,j-1) == 0 ...
%                     && splines_space(i,j-2) == 0 ...
%                     && splines_space(i,j-3) ~= 0 ...
%                     && splines_2_space(i,j-1) == 0 ...
%                     && splines_2_space(i,j-2) == 0 ...
%                     && splines_2_space(i,j-3) ~= 0
%                 % Don't calculate all the zeros values of spline
%                 break;
%             end
            
            g = @(s) Ndiff(p.T_space_ansatz, p.bSplineOrder_space_ansatz,0,i,(s==0)*(s+eps) + (s==1)*(s-eps) + (s > 0 && s < 1)*s);
            
            splines_space(i,j) = g(x(j));
            
            
         end
    else
        offset = 2^(resolution.x - p.refinementLevel_space);
        splines_space(i, offset+1:end) = splines_space(i-1, 1:end-offset);
    end
end

if false % (p.T_space == p.T_time)
    splines_time = splines_space;
else
    
    splines_time = zeros(size(u_full,4), length(t));
    
    for i=1:size(splines_time,1)
        if i < 5 || i > size(splines_time,1) - 5
            for j=1:size(splines_time,2)
%                 if j > 3 && splines_time(i,j-1) == 0 ...
%                         && splines_time(i,j-2) == 0 ...
%                         && splines_time(i,j-3) ~= 0 ...
%                         && splines_2_time(i,j-1) == 0 ...
%                         && splines_2_time(i,j-2) == 0 ...
%                         && splines_2_time(i,j-3) ~= 0
%                     % Don't calculate all the zeros values of spline
%                     break;
%                 end
                
                g = @(s) Ndiff(p.T_time_ansatz, p.bSplineOrder_time_ansatz,0,i,(s==0)*(s+eps) + (s==1)*(s-eps) + (s > 0 && s < 1)*s);
                
                splines_time(i,j) = g(t(j));
                
            end
        else
            offset = 2^(resolution.t - p.refinementLevel_time);
            splines_time(i, offset+1:end) = splines_time(i-1, 1:end-offset);
        end
    end
end



for i=1:size(u_full,1) % every space node in x
    if mod(i, size(u_full,1)/8) == 0
        disp(['Getting solution: ' num2str(i / size(u_full,1) * 100) '%'])
    end
    spline_x = splines_space(i,:);
    space_x_indices = find(spline_x, 1, 'first'):find(spline_x, 1, 'last');
    spline_x = spline_x(space_x_indices);
    
    for j=1:size(u_full, 2) % every space node in y
        
        spline_y = splines_space(j,:);
        space_y_indices = find(spline_y, 1, 'first'):find(spline_y, 1, 'last');
        spline_y = spline_y(space_y_indices);
        
        for k=1:size(u_full, 3) % every space node in z
            spline_z = splines_space(k,:);
            space_z_indices = find(spline_z, 1, 'first'):find(spline_z, 1, 'last');
            spline_z = spline_z(space_z_indices);
            
            spl_space = kron(spline_z, kron(spline_y, spline_x));

            for l=1:size(u_full,4) % every time node
                
                if abs(u_full(i,j,k,l)) < eps
                    continue;
                end
                
                spline_t = splines_time(l,:);
                time_indices = find(spline_t, 1, 'first'):find(spline_t, 1, 'last');
                spline_t = spline_t(time_indices);
                
                % Build the whole spline
                spline = kron(spline_t, spl_space);
                
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

