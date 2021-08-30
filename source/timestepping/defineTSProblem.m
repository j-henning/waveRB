% Define a time-stepping problem
% Input:
% dimension              - Space dimension (1, 2 or 3)
% f_time                 - Struct for the rhs time functions
% f_space                - Struct for the rhs space functions
% u_0                    - Function handle for u0
% u_1                    - Function handle for u1
% mu                     - Wave speed paremeter
% refinement_time        - Number of refinements in time
% refinement_space       - Number of refinements in space
% bSplineOrder           - Spline order in space (default = 3)
% offset_space           - Number of splines to skipt (default = [1 1])
% Output:
% problemTSConfiguration - Definition of a time-stepping problem
function [problemTSConfiguration] = defineTSProblem(dimension, f_time, ...
    f_space, u_0, u_1, mu, refinement_time, refinement_space, ...
    bSplineOrder, offset_space)
% Copy the first seven inputs
problemTSConfiguration.dimension = dimension;
problemTSConfiguration.f_time = f_time;
problemTSConfiguration.f_space = f_space;
problemTSConfiguration.u_0 = u_0;
problemTSConfiguration.u_1 = u_1;
problemTSConfiguration.mu = mu;
problemTSConfiguration.refinementLevel_time = refinement_time;
problemTSConfiguration.refinementLevel_space = refinement_space;

% Set the default parameters for the rest
problemTSConfiguration.bSplineOrder = 3;
problemTSConfiguration.offset_space = [1 1];

% Overwrite the default parameters if needed
if nargin > 8
    problemTSConfiguration.bSplineOrder = bSplineOrder;
end

if nargin > 9
    problemTSConfiguration.offset_space = offset_space;
end
end
