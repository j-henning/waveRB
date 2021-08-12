% Evaluates the solution of a problem on with the desired resolution at all
% query points in space and time on [0,1] x [0,1]
function [sol] = get1DsolutionSmart(p, U,  resolution,splines_time, splines_2_time, splines_space, splines_2_space)
mu = 1;
if isfield(p, 'mu')
    mu = p.mu;
end


x = linspace(0, 1, 2^resolution.x + 1);
t = linspace(0, 1, 2^resolution.t + 1);



u = reshape(full(U), numel(U) / p.dim_time, p.dim_time);

% Preparing boundary for plot
u_full = zeros(p.nsol_space, p.nsol_time);


u_full(1+p.offset_space_ansatz(1):p.nsol_space-p.offset_space_ansatz(2),...
    1+p.offset_time_ansatz(1):p.nsol_time-p.offset_time_ansatz(2)) = u;


sol = zeros(length(x), length(t));

for i=1:size(u_full,1) % every space node in x
    %     if mod(i, size(u_full,1)/8) == 0
    %         disp(['Plotting preparation: ' num2str(i / size(u_full,1) * 100) '%'])
    %     end
    spline_x = splines_space(i,:);
    spline_x_2 = splines_2_space(i,:);
    space_x_indices = min(find(spline_x, 1, 'first'), find(spline_x_2, 1, 'first')):...
        max(find(spline_x, 1, 'last'), find(spline_x_2, 1, 'last'));
    spline_x = spline_x(space_x_indices);
    spline_x_2 = spline_x_2(space_x_indices);
    
    for k=1:size(u_full,2) % every time node
        if abs(u_full(i,k)) < eps
            continue;
        end
        
        spline_t = splines_time(k,:);
        spline_t_2 = splines_2_time(k,:);
        time_indices = min(find(spline_t, 1, 'first'), find(spline_t_2, 1, 'first')):...
            max(find(spline_t, 1, 'last'), find(spline_t_2, 1, 'last'));
        spline_t = spline_t(time_indices);
        spline_t_2 = spline_t_2(time_indices);
        
        % Build the whole spline
        spline =  kron(spline_t_2, spline_x) ...
            -  kron(spline_t, mu * spline_x_2);
        
        spline = reshape(spline, ...
            [length(space_x_indices),  length(time_indices)]);
        sol(space_x_indices, time_indices) = ...
            sol(space_x_indices, time_indices) + ...
            u_full(i,k) * spline;
    end
end

end

