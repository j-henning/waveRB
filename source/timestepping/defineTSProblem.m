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
    