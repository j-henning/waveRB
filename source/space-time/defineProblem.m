% Defines a problem which can than be created through the createProblem
% function
% Input:
% dimension            - Dimension of the problem (1, 2 or 3)
% f_time               - Cell array containing the rhs time component
% f_space              - Cell array containing the rhs space component
% u_0                  - Function handle of u0
% u_1                  - Function handle of u1
% mu                   - Wave speed parameter
% refinement_time      - Number of refinements in time
% refinement_space     - Number of refinements in space
% bSplineOrder_time    - Spline order in time (default = 3)
% bSplineOrder_space   - Spline order in space (default = 3)
% offset_time_ansatz   - Number of splines to skip (default = [0 2])
% offset_time_test     - Number of splines to skip (default = [0 2])
% offset_space_ansatz  - Number of splines to skip (default = [1 1])
% offset_space_test    - Number of splines to skip (default = [1 1])
% Output:
% problemConfiguration - A struct containing a the data need for the
%                        createProblem function
function [problemConfiguration] = defineProblem(dimension, f_time, ...
    f_space, u_0, u_1, mu, refinement_time, refinement_space, ...
    bSplineOrder_time, bSplineOrder_space, offset_time_ansatz, ...
    offset_time_test, offset_space_ansatz, offset_space_test)
% Copy the first seven inputs
problemConfiguration.dimension = dimension;
problemConfiguration.f_time = f_time;
problemConfiguration.f_space = f_space;
problemConfiguration.u_0 = u_0;
problemConfiguration.u_1 = u_1;
problemConfiguration.mu = mu;
problemConfiguration.refinementLevel_time = refinement_time;
problemConfiguration.refinementLevel_space = refinement_space;

% Set the default parameters for the rest
problemConfiguration.bSplineOrder_time = 3;
problemConfiguration.bSplineOrder_space = 3;
problemConfiguration.offset_time_ansatz = [0 2];
problemConfiguration.offset_time_test = [0 2];
problemConfiguration.offset_space_ansatz = [1 1];
problemConfiguration.offset_space_test = [1 1];

% Overwrite the default parameters if needed
if nargin > 8
    problemConfiguration.bSplineOrder_time = bSplineOrder_time;
end

if nargin > 9
    problemConfiguration.bSplineOrder_space = bSplineOrder_space;
end

if nargin > 10
    problemConfiguration.offset_time_ansatz = offset_time_ansatz;
end

if nargin > 11
    problemConfiguration.offset_time_test = offset_time_test;
end

if nargin > 12
    problemConfiguration.offset_space_ansatz = offset_space_ansatz;
end

if nargin > 13
    problemConfiguration.offset_space_test = offset_space_test;
end

end
