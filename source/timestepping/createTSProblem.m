% Creates a time-stepping problem
% Input:
% problemTSConfiguration - Configuration containing all the needed data
% Output:
% problemTS              - Time-stepping problem containing all matrices,
%                          initial values as well as a rhs
function [problemTS] = createTSProblem(problemTSConfiguration)
switch problemTSConfiguration.dimension
    case 1
        problemTS = create1DTSWaveProblem(problemTSConfiguration);
    case 2
        problemTS = create2DTSWaveProblem(problemTSConfiguration);
    case 3
        problemTS = create3DTSWaveProblem(problemTSConfiguration);
end
end