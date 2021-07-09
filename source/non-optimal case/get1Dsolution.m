% Evaluates the solution of a problem on with the desired resolution at all
% query points in space and time on [0,1] x [0,1]
function [sol] = get1Dsolution(p, U,  resolution)
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

% Precalculate the splines

splines_space = zeros(size(u_full,1), length(x));


for i=1:size(splines_space,1)
    %     if i < 5 || i > size(splines_space,1) - 5
    for j=1:size(splines_space,2)
        
        %             if j > 3 && splines_space(i,j-1) == 0 ...
        %                         && splines_space(i,j-2) == 0 ...
        %                         && splines_space(i,j-3) ~= 0 ...
        %                         && splines_2_space(i,j-1) == 0 ...
        %                         && splines_2_space(i,j-2) == 0 ...
        %                         && splines_2_space(i,j-3) ~= 0
        %                     % Don't calculate all the zeros values of spline
        %                     break;
        %             end
        
        g = @(s) Ndiff(p.T_space_ansatz, p.bSplineOrder_space_ansatz,0,i,(s==0)*(s+eps) + (s==1)*(s-eps) + (s > 0 && s < 1)*s);
   
        
        splines_space(i,j) = g(x(j));

        
    end
    %     else
    %         offset = 2^(resolution.x - p.refinementLevel_space);
    %         splines_space(i, offset+1:end) = splines_space(i-1, 1:end-offset);
    %         splines_2_space(i, offset+1:end) = splines_2_space(i-1, 1:end-offset);
    %     end
end

% if (p.T_space == p.T_time)
%     splines_time = splines_space;
%     splines_2_time = splines_2_space;
% else
    
    splines_time = zeros(size(u_full,2), length(t));

    for i=1:size(splines_time,1)
        %         if i < 5 || i > size(splines_time,1) - 5
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
        %         else
        %             offset = 2^(resolution.t - p.refinementLevel_time);
        %             splines_time(i, offset+1:end) = splines_time(i-1, 1:end-offset);
        %             splines_2_time(i, offset+1:end) = splines_2_time(i-1, 1:end-offset);
        %         end
    end
% end



for i=1:size(u_full,1) % every space node in x
    %     if mod(i, size(u_full,1)/8) == 0
    %         disp(['Plotting preparation: ' num2str(i / size(u_full,1) * 100) '%'])
    %     end
    spline_x = splines_space(i,:);

    space_x_indices = find(spline_x, 1, 'first'):find(spline_x, 1, 'last');
    spline_x = spline_x(space_x_indices);
 
    
    for k=1:size(u_full,2) % every time node
        if abs(u_full(i,k)) < eps
            continue;
        end
        
        spline_t = splines_time(k,:);
        
        time_indices = find(spline_t, 1, 'first'):find(spline_t, 1, 'last');
        spline_t = spline_t(time_indices);
     
        
        % Build the whole spline
        spline = kron(spline_t, spline_x);
        
        spline = reshape(spline, ...
            [length(space_x_indices),  length(time_indices)]);
        sol(space_x_indices, time_indices) = ...
            sol(space_x_indices, time_indices) + ...
            u_full(i,k) * spline;
    end
end

end

