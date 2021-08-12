function [solution, splines] = getSolution(problem, U, resolution, splines)
if nargin < 4
    % Compute the splines if not provided
    p = problem;
    x = linspace(0, 1, 2^resolution.x + 1);
    t = linspace(0, 1, 2^resolution.t + 1);
    
    splines.space = zeros(p.nsol_space, length(x));
    splines.space2 = zeros(p.nsol_space, length(x));
    
    for i=1:size(splines.space,1)
        for j=1:size(splines.space,2)
            
            g = @(s) Ndiff(p.T_space, p.bSplineOrder_space,0,i,(s==0)*(s+eps) + (s==1)*(s-eps) + (s > 0 && s < 1)*s);
            h = @(s) Ndiff(p.T_space, p.bSplineOrder_space,2,i,(s==0)*(s+eps) + (s==1)*(s-eps) + (s > 0 && s < 1)*s);
            
            splines.space(i,j) = g(x(j));
            
            splines.space2(i,j) = h(x(j));
            
        end
    end
    
    if (problem.T_space == problem.T_time)
        splines.time = splines.space;
        splines.time2 = splines.space2;
    else
        
        splines.time = zeros(p.nsol_time, length(t));
        splines.time2 = zeros(p.nsol_time, length(t));
        
        for i=1:size(splines.time,1)
            for j=1:size(splines.time,2)
                g = @(s) Ndiff(p.T_time, p.bSplineOrder_time,0,i,(s==0)*(s+eps) + (s==1)*(s-eps) + (s > 0 && s < 1)*s);
                h = @(s) Ndiff(p.T_time, p.bSplineOrder_time,2,i,(s==0)*(s+eps) + (s==1)*(s-eps) + (s > 0 && s < 1)*s);
                
                splines.time(i,j) = g(t(j));
                splines.time2(i,j) = h(t(j));
                
            end
            
        end
    end
    
end

% Call the method for the right dimension
switch problem.dimension
    case 1
        solution = get1Dsolution(p, U,  resolution, splines);
    case 2
        error('Not yet implemented')
    case 3
        error('Not yet implemented')
end
end