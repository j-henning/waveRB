% Evaluates the solution of a time-stepping problem on a grid
% Input:
% problemTS   - Time-stepping problem
% UTS         - Solution vector
% resolution  - Resolution struct
% splines     - Spline data
% Output:
% solutionsTS - Solution
function [solutionTS] = getTSSolution(problemTS, UTS, resolution, splines)
if nargin < 4
    % Compute the splines if not provided
    p = problemTS;
    x = linspace(0, 1, 2^resolution.x + 1);
    
    switch problemTS.dimension
        case 1
            spaceDimension = size(UTS,1)+2;
        case 2
            spaceDimension = round(size(UTS,1)^0.5+2);
        case 3
            spaceDimension = round(size(UTS,1)^(1/3)+2);
    end
    
    splines.space = zeros(spaceDimension, length(x));
    
    for i=1:size(splines.space,1)
        for j=1:size(splines.space,2)
            % TODO: Optimize this
            g = @(s) Ndiff(p.T_space, p.bSplineOrder,0,i,(s==0)*(s+eps) + (s==1)*(s-eps) + (s > 0 && s < 1)*s);
            splines.space(i,j) = g(x(j));
        end
    end
end


switch problemTS.dimension
    case 1
        solutionTS = get1DTSSolution(problemTS, UTS, resolution, splines);
    case 2
        solutionTS = get2DTSSolution(problemTS, UTS, resolution, splines);
    case 3
        solutionTS = get3DTSSolution(problemTS, UTS, resolution, splines); 
end
end