% Evaluates the solution of a problem on with the desired resolution at all
% query points in space and time
% Input:
% problem    - Problem
% U          - Solution vector
% resolution - Resolution
% splines    - Splines calculated by computeSplines (optional)
% Output:
% solution   - Solution
% splines    - Spline data
function [solution, splines] = getSolution(problem, U, resolution, splines)

if nargin < 4
    splines = computeSplines(problem, resolution);
end

% Call the method for the right dimension
switch problem.dimension
    case 1
        solution = get1Dsolution(problem, U,  resolution, splines);
    case 2
        solution = get2Dsolution(problem, U,  resolution, splines);
    case 3
        solution = get3Dsolution(problem, U,  resolution, splines);
end
end