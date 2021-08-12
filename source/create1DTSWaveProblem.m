function [problemTS] = create1DTSWaveProblem(problemConfiguration)
problemTS = problemConfiguration;

problemTS.T_space = InitUniNodes(problemTS.refinementLevel_space, ...
    problemTS.bSplineOrder);

% Compute the mass and stiffness matrix
problemTS.M = StiffMat(problemTS.T_space, ... % Ansatz discretization
    problemTS.T_space, ... % Test discretization
    problemTS.bSplineOrder, ... % Ansatz spline order
    problemTS.bSplineOrder, ... % Test spline order
    problemTS.offset_space,... % Ansatz offset
    problemTS.offset_space, ... % Test offset
    [0, 0] ... % diff
    );

problemTS.A = problemTS.mu * StiffMat(problemTS.T_space, ... % Ansatz discretization
    problemTS.T_space, ... % Test discretization
    problemTS.bSplineOrder, ... % Ansatz spline order
    problemTS.bSplineOrder, ... % Test spline order
    problemTS.offset_space,... % Ansatz offset
    problemTS.offset_space, ... % Test offset
    [1, 1] ... % diff
    );


problemTS.K = 2^ problemTS.refinementLevel_time; % Number of time steps
problemTS.tau = 1 / (problemTS.K-1);

% Compute the right hand side
problemTS.ntest = length(problemTS.T_space) - problemTS.bSplineOrder;
problemTS.Fvals = zeros(problemTS.ntest-2,problemTS.K);

if ~isempty(problemTS.f_time)
    for k=1:problemTS.K
        for i=1:length(problemTS.f_time)
            temp = zeros(problemTS.ntest-2,1);
            for j=1+problemTS.offset_space(1):problemTS.ntest-problemTS.offset_space(2)
                g = @(x) problemTS.f_space{i}(x) * Ndiff(problemTS.T_space, problemTS.bSplineOrder,0,j,x);
                for step = 0:problemTS.bSplineOrder-1
                    temp(j-1) = temp(j-1) ...
                        + Gaussq(problemTS.T_space(j+step),problemTS.T_space(j+step+1),g,3);
                end
            end
            problemTS.Fvals(:, k) = problemTS.Fvals(:,k) + temp .* problemTS.f_time{i}((k-1)/(problemTS.K-1));
        end
    end
end

% Compute the inital values
problemTS.U0 = zeros(problemTS.ntest-2,1);
problemTS.U1 = zeros(problemTS.ntest-2,1);

if ~isempty(problemTS.u_0)
    for j=1+problemTS.offset_space(1):problemTS.ntest-problemTS.offset_space(2) 
        g = @(x) problemTS.u_0(x) .* Ndiff(problemTS.T_space, problemTS.bSplineOrder,0,j,x);
        for step = 0:problemTS.bSplineOrder-1
            problemTS.U0(j-1) = problemTS.U0(j-1) ...
                + Gaussq(problemTS.T_space(j+step),problemTS.T_space(j+step+1),g,3);
        end
    end
end


if ~isempty(problemTS.u_1)
    for j=1+problemTS.offset_space(1):problemTS.ntest-problemTS.offset_space(2)
        h = @(x) problemTS.u_1(x) .* Ndiff(problemTS.T_space, problemTS.bSplineOrder,0,j,x);
        for step = 0:problemTS.bSplineOrder-1
            problemTS.U1(j-1) = problemTS.U1(j-1) ...
                + Gaussq(problemTS.T_space(j+step),problemTS.T_space(j+step+1),h,3);
        end
    end
end

end