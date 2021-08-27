function [problemTS] = create2DTSWaveProblem(problemConfiguration)
problemTS = problemConfiguration;
problemTS.T_space = InitUniNodes(problemTS.refinementLevel_space, ...
    problemTS.bSplineOrder);

% Compute the mass and stiffness matrix
M1D = StiffMat(problemTS.T_space, ... % Ansatz discretization
    problemTS.T_space, ... % Test discretization
    problemTS.bSplineOrder, ... % Ansatz spline order
    problemTS.bSplineOrder, ... % Test spline order
    problemTS.offset_space,... % Ansatz offset
    problemTS.offset_space, ... % Test offset
    [0, 0] ... % diff
    );

A = problemTS.mu * StiffMat(problemTS.T_space, ... % Ansatz discretization
    problemTS.T_space, ... % Test discretization
    problemTS.bSplineOrder, ... % Ansatz spline order
    problemTS.bSplineOrder, ... % Test spline order
    problemTS.offset_space,... % Ansatz offset
    problemTS.offset_space, ... % Test offset
    [1, 1] ... % diff
    );

problemTS.M = kron(M1D, M1D);
problemTS.A = kron(A, M1D) + kron(M1D,A);

problemTS.K = 2^ problemTS.refinementLevel_time; % Number of time steps
problemTS.tau = 1 / (problemTS.K-1);

% Compute the right hand side
problemTS.ntest = length(problemTS.T_space) - problemTS.bSplineOrder;
problemTS.Fvals = zeros((problemTS.ntest-sum(problemTS.offset_space))^2,problemTS.K);

if ~isempty(problemTS.f_time)
    if size(problemTS.f_space,2) == 2
        % The right hand side is provided in a nested form, so we can use
        % the much faster integration
        
        for i=1:length(problemTS.f_time)
            
            
            % first space component
            rhs_space_x = zeros(problemTS.ntest-(problemTS.offset_space(1)+problemTS.offset_space(2)),1);
            
            for j=1+problemTS.offset_space(1):problemTS.ntest-problemTS.offset_space(2)
                g = @(x) problemTS.f_space{i,1}(x) * Ndiff(problemTS.T_space, problemTS.bSplineOrder,0,j,(x==0)*(x+eps) + (x==1)*(x-eps) + (x > 0 && x < 1)*x);
                for step = 0:problemTS.bSplineOrder-1
                    rhs_space_x(j-problemTS.offset_space(1)) = rhs_space_x(j-problemTS.offset_space(1)) ...
                        + Gaussq(problemTS.T_space(j+step),problemTS.T_space(j+step+1),g,5);
                end
            end
            
            % second space component
            problemTS.ntest_space = length(problemTS.T_space) - problemTS.bSplineOrder;
            rhs_space_y = zeros(problemTS.ntest-(problemTS.offset_space(1)+problemTS.offset_space(2)),1);
            
            for j=1+problemTS.offset_space(1):problemTS.ntest-problemTS.offset_space(2)
                g = @(y) problemTS.f_space{i,2}(y) * Ndiff(problemTS.T_space, problemTS.bSplineOrder,0,j,(y==0)*(y+eps) + (y==1)*(y-eps) + (y > 0 && y < 1)*y);
                for step = 0:problemTS.bSplineOrder-1
                    rhs_space_y(j-problemTS.offset_space(1)) = rhs_space_y(j-problemTS.offset_space(1)) ...
                        + Gaussq(problemTS.T_space(j+step),problemTS.T_space(j+step+1),g,5);
                end
            end
            
            
            for k=1:problemTS.K
                problemTS.Fvals(:, k) = problemTS.Fvals(:,k) + kron(rhs_space_y, rhs_space_x).* problemTS.f_time{i}((k-1)/(problemTS.K-1));
            end
        end
        
    else % The right hand side is in nested form, so we have to do more work
        
        for i=1:length(problemTS.f_time)
            
            % Integrate over the space
            F_space = zeros(problemTS.ntest-sum(problemTS.offset_space), ...
                problemTS.ntest-sum(problemTS.offset_space));
            
            
            for j=1+problemTS.offset_space(1):problemTS.ntest-problemTS.offset_space(2)
                for k=1+problemTS.offset_space(1):problemTS.ntest-problemTS.offset_space(2)
                    
                    
                    
                    g = @(z)@(y)@(x) problemTS.f_space{i}(x,y,z) .* ...
                        Ndiff(problemTS.T_space, problemTS.bSplineOrder,0,j,x) ...
                        .* Ndiff(problemTS.T_space, problemTS.bSplineOrder,0,k,y) ...
                        .* Ndiff(problemTS.T_space, problemTS.bSplineOrder,0,l,z);
                    
                    % Integrate
                    for stepX = 0:problemTS.bSplineOrder-1
                        for stepY = 0:problemTS.bSplineOrder-1
                            
                            
                            a = [problemTS.T_space(j+stepX), problemTS.T_space(k+stepY)];
                            b = [problemTS.T_space(j+stepX+1), problemTS.T_space(k+stepY+1)];
                            
                            % Check if f is zero in the center of the
                            % inteval, as well as on the boundary
                            % if yes assume that the integral is zero.
                            % This is done for performance reasons
                            center = 0.5 * (a+b);
                            if abs(problemTS.f_space{i}(center(1),center(2))) < eps && ...
                                    abs(problemTS.f_space{i}(a(1),a(2))) < eps && ...
                                    abs(problemTS.f_space{i}(a(1),b(2))) < eps && ...
                                    abs(problemTS.f_space{i}(b(1),b(2))) < eps && ...
                                    abs(problemTS.f_space{i}(b(1),a(2))) < eps
                                continue;
                            end
                            
                            % Integrate from a to b
                            F_space(j-problemTS.offset_space(1), ...
                                k-problemTS.offset_space(1)) = ...
                                F_space(j-problemTS.offset_space(1), ...
                                k-problemTS.offset_space(1))  + ...
                                Gaussq(a,b,g,3);
                            
                        end
                    end
                    
                    
                end
            end
            
            for k=1:problemTS.K
                problemTS.Fvals(:, k) = problemTS.Fvals(:,k) + F_space(:) .* problemTS.f_time{i}((k-1)/(problemTS.K-1));
            end
        end
    end
end


% Compute the inital values
U0 = zeros(problemTS.ntest-sum(problemTS.offset_space), ...
    problemTS.ntest-sum(problemTS.offset_space), 1);
U1 = U0;

if ~isempty(problemTS.u_0) || ~isempty(problemTS.u_1)
    for j=1+problemTS.offset_space(1):problemTS.ntest-problemTS.offset_space(2)
        for k=1+problemTS.offset_space(1):problemTS.ntest-problemTS.offset_space(2)
            g = @(y)@(x) problemTS.u_0(x,y,z) .* ...
                Ndiff(problemTS.T_space, problemTS.bSplineOrder,0,j,x) .* ...
                Ndiff(problemTS.T_space, problemTS.bSplineOrder,0,k,y);
            h = @(y)@(x) problemTS.u_1(x,y,z) .* ...
                Ndiff(problemTS.T_space, problemTS.bSplineOrder,0,j,x) .* ...
                Ndiff(problemTS.T_space, problemTS.bSplineOrder,0,k,y);
            
            for stepX = 0:problemTS.bSplineOrder-1
                for stepY = 0:problemTS.bSplineOrder-1
                    % Integrate from a to b
                    a = [problemTS.T_space(j+stepX); problemTS.T_space(k+stepY)];
                    b = [problemTS.T_space(j+stepX+1); problemTS.T_space(k+stepY+1)];
                    center = 0.5 * (a+b);
                    if ~isempty(problemTS.u_0)
                        if abs(problemTS.u_0(center(1),center(2))) > eps || ...
                                abs(problemTS.u_0(a(1),a(2))) > eps || ...
                                abs(problemTS.u_0(a(1),b(2))) > eps || ...
                                abs(problemTS.u_0(b(1),b(2))) > eps || ...
                                abs(problemTS.u_0(b(1),a(2))) > eps
                            U0(j-problemTS.offset_space(1), ...
                                k-problemTS.offset_space(1)) = ...
                                U0(j-problemTS.offset_space(1), ...
                                k-problemTS.offset_space(1)) + Gaussq(a,b,g,5);
                        end
                    end
                    
                    if ~isempty(problemTS.u_1)
                        if abs(problemTS.u_1(center(1),center(2))) > eps || ...
                                abs(problemTS.u_1(a(1),a(2))) > eps || ...
                                abs(problemTS.u_1(a(1),b(2))) > eps || ...
                                abs(problemTS.u_1(b(1),b(2))) > eps || ...
                                abs(problemTS.u_1(b(1),a(2))) > eps
                            U1(j-problemTS.offset_space(1), ...
                                k-problemTS.offset_space(1)) = ...
                                U1(j-problemTS.offset_space(1), ...
                                k-problemTS.offset_space(1)) + Gaussq(a,b,h,5);
                        end
                    end
                end
            end
        end
    end
end

problemTS.U0 = U0(:);
problemTS.U1 = U1(:);
end