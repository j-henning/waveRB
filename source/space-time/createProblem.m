% This functions returns a 'problem' struct consisting of of all matrices 
% as well a rhs for a given problemConfiguration. Internally,
% create1DWaveProblem, create2DWaveProblem or create3DWaveProblem is
% called.
% Input:
% - problemConfiguration: Can be created with defineProblem(...)
% Output:
% - problem: Contains all data from the configuration as well the matrices
%            and the rhs
function [problem] = createProblem(problemConfiguration)
    switch problemConfiguration.dimension
        case 1
            problem = create1DWaveProblem(problemConfiguration);
        case 2
            problem = create2DWaveProblem(problemConfiguration);
        case 3
            problem = create3DWaveProblem(problemConfiguration);
        otherwise 
            error([num2str(problemConfiguration.dimension) ...
                ' dimensions are not supported'])
    end
end