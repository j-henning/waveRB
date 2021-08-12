function [problemTS] = createTSProblem(problemTSConfiguration)
    switch problemTSConfiguration.dimension
        case 1
            problemTS = create1DTSWaveProblem(problemTSConfiguration);
        case 2
            error('Not yet implemented')
        case 3
            error('Not yet implemented')
    end
end