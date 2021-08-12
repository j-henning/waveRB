% Evaluates the solution of a problem on with the desired resolution at all
% query points in space and time on [0,1] x [0,1]
function [splines_time, splines_2_time, splines_space, splines_2_space] = ...
    precompute1DSplines(p,  resolution)

x = linspace(0, 1, 2^resolution.x + 1);
t = linspace(0, 1, 2^resolution.t + 1);


% Precalculate the splines



splines_space = zeros(p.nsol_space, length(x));
splines_2_space = zeros(p.nsol_space, length(x));

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
        
        g = @(s) Ndiff(p.T_space, p.bSplineOrder_space,0,i,(s==0)*(s+eps) + (s==1)*(s-eps) + (s > 0 && s < 1)*s);
        h = @(s) Ndiff(p.T_space, p.bSplineOrder_space,2,i,(s==0)*(s+eps) + (s==1)*(s-eps) + (s > 0 && s < 1)*s);
        
        splines_space(i,j) = g(x(j));
        
        splines_2_space(i,j) = h(x(j));
        
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
    
    splines_time = zeros(p.nsol_time, length(t));
    splines_2_time = zeros(p.nsol_time, length(t));
    
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
            
            g = @(s) Ndiff(p.T_time, p.bSplineOrder_time,0,i,(s==0)*(s+eps) + (s==1)*(s-eps) + (s > 0 && s < 1)*s);
            h = @(s) Ndiff(p.T_time, p.bSplineOrder_time,2,i,(s==0)*(s+eps) + (s==1)*(s-eps) + (s > 0 && s < 1)*s);
            
            splines_time(i,j) = g(t(j));
            splines_2_time(i,j) = h(t(j));
            
        end
        %         else
        %             offset = 2^(resolution.t - p.refinementLevel_time);
        %             splines_time(i, offset+1:end) = splines_time(i-1, 1:end-offset);
        %             splines_2_time(i, offset+1:end) = splines_2_time(i-1, 1:end-offset);
        %         end
    end
% end





end

