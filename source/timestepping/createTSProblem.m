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