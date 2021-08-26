function [solution, splines] = getSolution(problem, U, resolution, splines)

if nargin < 4
    splines = computeSplines(problem, resolution);
end


% 
% splines.space = sparse(splines.space);
% splines.space2 = sparse(splines.space2);
% splines.time = sparse(splines.time);
% splines.time2 = sparse(splines.time2);

% Call the method for the right dimension
switch problem.dimension
    case 1
        solution = get1Dsolution(problem, U,  resolution, splines);
    case 2
        error('Not yet implemented')
    case 3
        solution = get3Dsolution(problem, U,  resolution, splines);
end
end