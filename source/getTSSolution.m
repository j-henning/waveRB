function [solutionTS] = getTSSolution(problemTS, UTS, resolution, splines)
if nargin < 4
    % Compute the splines if not provided
    p = problemTS;
    x = linspace(0, 1, 2^resolution.x + 1);

    
    splines.space = zeros(size(UTS,1)+2, length(x));
    splines.space2 = zeros(size(UTS,1)+2, length(x));
    
    for i=1:size(splines.space,1)
        for j=1:size(splines.space,2)
            
            g = @(s) Ndiff(p.T_space, p.bSplineOrder,0,i,(s==0)*(s+eps) + (s==1)*(s-eps) + (s > 0 && s < 1)*s);
            h = @(s) Ndiff(p.T_space, p.bSplineOrder,2,i,(s==0)*(s+eps) + (s==1)*(s-eps) + (s > 0 && s < 1)*s);
            
            splines.space(i,j) = g(x(j));
            
            splines.space2(i,j) = h(x(j));
            
        end
    end    
end


switch problemTS.dimension
    case 1
        solutionTS = get1DTSSolution(problemTS, UTS, resolution, splines);
    case 2
        error('Not yet implemented')
    case 3
        error('Not yet implemented')
    
end
end