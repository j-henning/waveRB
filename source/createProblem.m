function [problem] = createProblem(problemConfiguration)
    switch problemConfiguration.dimension
        case 1
            problem = create1DWaveProblem(problemConfiguration);
        case 2
            problem = create2DWaveProblemImproved(problemConfiguration);
        case 3
            problem = create3DWaveProblem(problemConfiguration);
    end
end