% Precompute the time and space splines as well as their second derivatives
% Input:
% problem    - A problem which specifies all dimension
% resolution - Resolution of the space and time splines. The splines get
%              evaluated on x = linspace(0, 1, 2^resolution.x + 1) and
%              t = linspace(0, 1, 2^resolution.t + 1);
% Output:
% splines    - Struct consisting of splines.space, splines.space2 as well
%              as splines.time, splines.time2
function [splines] = computeSplines(problem, resolution)
% Compute the splines if not provided
x = linspace(0, 1, 2^resolution.x + 1);
t = linspace(0, 1, 2^resolution.t + 1);

splines.space = zeros(problem.nsol_space, length(x));
splines.space2 = zeros(problem.nsol_space, length(x));

for i=1:size(splines.space,1)
    % Find the interval in which the splines 'lives'
    intervalStart = problem.T_space(i);
    intervalEnd = problem.T_space(i+problem.bSplineOrder_space);
    
    indices = find(x >= intervalStart, 1, 'first'):find(x <= intervalEnd,1, 'last');
    
    
    for j=indices
        
        g = @(s) Ndiff(problem.T_space, problem.bSplineOrder_space,0,i,(s==0)*(s+eps) + (s==1)*(s-eps) + (s > 0 && s < 1)*s);
        h = @(s) Ndiff(problem.T_space, problem.bSplineOrder_space,2,i,(s==0)*(s+eps) + (s==1)*(s-eps) + (s > 0 && s < 1)*s);
        
        splines.space(i,j) = g(x(j));
        
        splines.space2(i,j) = h(x(j));
        
    end
end

if length(problem.T_space) == length(problem.T_time)
    if problem.T_space == problem.T_time
    splines.time = splines.space;
    splines.time2 = splines.space2;
    end
else
    
    splines.time = zeros(problem.nsol_time, length(t));
    splines.time2 = zeros(problem.nsol_time, length(t));
    
    for i=1:size(splines.time,1)
        for j=1:size(splines.time,2)
            g = @(s) Ndiff(problem.T_time, problem.bSplineOrder_time,0,i,(s==0)*(s+eps) + (s==1)*(s-eps) + (s > 0 && s < 1)*s);
            h = @(s) Ndiff(problem.T_time, problem.bSplineOrder_time,2,i,(s==0)*(s+eps) + (s==1)*(s-eps) + (s > 0 && s < 1)*s);
            
            splines.time(i,j) = g(t(j));
            splines.time2(i,j) = h(t(j));
            
        end
        
    end
end

splines.time = sparse(splines.time);
splines.space = sparse(splines.space);
splines.time2 = sparse(splines.time2);
splines.space2 = sparse(splines.space2);
end