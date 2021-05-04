clear
close all
clc

addpath('../splines');
addpath('../splines/Utilities');
addpath('../solvers');


% 
% load('2DreferenceSol.mat');


%% Define the example
% Example 1
n = 3; m = 2;
f_time = {@(t) 2+ 0*t; @(t) t.^2*(n^2 + m^2)*pi^2;};
f_space_x = {@(x) sin(n*pi*x), @(x) sin(n*pi*x)};
f_space_y = {@(y) sin(m*pi*y), @(y) sin(m*pi*y)};
u_analytical = @(x,y,t) t.^2 .* sin(n*pi*x).*sin(m*pi*y);

initial_conditions = false;

% Example 2
% f_time = {@(t) 0*t};
% f_space_x = {@(x) 0*x};
% f_space_y = {@(y) 0*y};
% u0_x = @(x) x* (x < 0.5) + (1-x).*(x > 0.5);
% u0_y = @(y)  y* (y < 0.5) + (1-y).*(y > 0.5);
% u1_x = @(x) 0*x; u1_y = @(y) 0 *y;
% path = '2D-example2';
% initial_conditions = true;

% Example 3
% f_time = {@(t) 0*t};
% f_space_x = {@(x) 0*x};
% f_space_y = {@(y) 0*y};
% u0_x = @(x) 1 * (x > 0.25) .* (x < 0.75);
% u0_y = @(y) 1 * (y > 0.25) .* (y < 0.75);
% u1_x = @(x) 0*x; u1_y = @(y) 0 *y;
% path = '2D-example3';
% initial_conditions = true;



%% Define the ansatz and trial splines
laplace = false;
splineOrder = 3; % at least four if laplace is true, otherwise at least 2

%% Define the resolutions to test
refinements = (1:8)';


% Define the resolution for the L2 error calculation (and plotting)
x = linspace(0,1, 2^9+1);
y = linspace(0,1, 2^9+1);
Kref = 2^9+1;

% Compute the analytical solution
% load(path)

[Xq, Yq, Tq] = ndgrid(x, y, linspace(0,1, Kref));
sol_ana = u_analytical(Xq, Yq, Tq);

times = zeros(size(refinements));
l2error = zeros(size(refinements));
        

for refinement = refinements'
        fprintf('Refinement: %d\n', refinement);
        
        K = 2^ refinement; % Number of time steps
        tau = 1 / (K-1);
        sol = zeros(length(x), length(y), K);
        
        T_space_test = InitUniNodes(refinement, splineOrder);
        T_space_ansatz = InitUniNodes(refinement, splineOrder);
        
        %% Calculation of the space matrices
        
        ansatz_offset = [1 1];
        test_offset = [1 1];
        
        % 1D Matrices
        M1D = StiffMat(T_space_ansatz, ... % Ansatz discretization
            T_space_test, ... % Test discretization
            splineOrder, ... % Ansatz spline order
            splineOrder, ... % Test spline order
            ansatz_offset,... % Ansatz offset
            test_offset, ... % Test offset
            [0, 0] ... % diff
            );
        
        Q = StiffMat(T_space_ansatz, ... % Ansatz discretization
            T_space_test, ... % Test discretization
            splineOrder, ... % Ansatz spline order
            splineOrder, ... % Test spline order
            ansatz_offset,... % Ansatz offset
            test_offset, ... % Test offset
            [1, 1] ... % diff
            );
        
        if laplace
            O = StiffMat(T_space_ansatz, ... % Ansatz discretization
                T_space_test, ... % Test discretization
                splineOrder, ... % Ansatz spline order
                splineOrder, ... % Test spline order
                ansatz_offset,... % Ansatz offset
                test_offset, ... % Test offset
                [2, 0] ... % diff
                );
            
            
            
            S = StiffMat(T_space_ansatz, ... % Ansatz discretization
                T_space_test, ... % Test discretization
                splineOrder, ... % Ansatz spline order
                splineOrder, ... % Test spline order
                ansatz_offset,... % Ansatz offset
                test_offset, ... % Test offset
                [3, 1] ... % diff
                );
        end
        % 2D Matrces
        dimM = size(M1D,1);
        if laplace
            
%             funM = @(x) reshape(M1D * reshape(x, [dimM dimM]) * O' ...
%                 + O * reshape(x, [dimM dimM]) * M1D', [], 1);
%             
%             funA = @(x) reshape(M1D * reshape(x, [dimM dimM]) * S' ...
%                 + Q * reshape(x, [dimM dimM]) * O' ...
%                 + S * reshape(x, [dimM dimM]) * M1D' ...
%                 + O * reshape(x, [dimM dimM]) * Q', [], 1);
            M = kron(O, M1D) + kron(M1D, O);
            A = kron(S, M1D) + kron(O, Q) ...
                + kron(M1D, S) + kron (Q, O);
        else
%             funM = @(x) reshape(M1D * reshape(x, [dimM dimM]) * M1D', [], 1);
%             
%             funA = @(x)  reshape(M1D * reshape(x, [dimM dimM]) * Q' ...
%                 + Q * reshape(x, [dimM dimM]) * M1D', [], 1);
            M = kron(M1D, M1D);
            A = kron(Q, M1D) + kron(M1D,Q);
        end

        
        %% Calculation of the right hand side
        ntest = length(T_space_test) - splineOrder;
        Fvals = zeros((ntest-sum(test_offset)).^2,K);
        
        
        for k=1:K
            for i=1:length(f_time)
                % First space dimension
                tempX = zeros(ntest-sum(test_offset),1);
                for j=1+test_offset(1):ntest-test_offset(2) % Use a offset of [1 1]
                    g = @(x) f_space_x{i}(x) * Ndiff(T_space_test, splineOrder,0,j,x);
                    for step = 0:splineOrder-test_offset(1)%bSplineOrder_space_test-1
                        tempX(j-test_offset(1)) = tempX(j-test_offset(1)) ...
                            + Gaussq(T_space_test(j+step),T_space_test(j+step+1),g,3);
                    end
                end
                tempY = zeros(ntest-sum(test_offset),1);
                for j=1+test_offset(1):ntest-test_offset(2) % Use a offset of [1 1]
                    g = @(y) f_space_y{i}(y) * Ndiff(T_space_test, splineOrder,0,j,y);
                    for step = 0:splineOrder-test_offset(1)%bSplineOrder_space_test-1
                        tempY(j-test_offset(1)) = tempY(j-test_offset(1)) ...
                            + Gaussq(T_space_test(j+step),T_space_test(j+step+1),g,3);
                    end
                end
                
                Fvals(:, k) = Fvals(:,k) + kron(tempY, tempX) .* f_time{i}((k-1)/(K-1));
            end
        end
        
        
        %% Calculation of u0 and u1
        
        if initial_conditions
            tempXu0 = zeros(ntest-sum(test_offset),1);
            tempXu1 = zeros(ntest-sum(test_offset),1);
            for j=1+test_offset(1):ntest-test_offset(2)
                g = @(x) u0_x(x) * Ndiff(T_space_test, splineOrder,0,j,x);
                h = @(x) u1_x(x) * Ndiff(T_space_test, splineOrder,0,j,x);
                
                for step = 0:splineOrder-1
                    tempXu0(j-test_offset(1)) = tempXu0(j-test_offset(1)) ...
                        + Gaussq(T_space_test(j+step),T_space_test(j+step+1),g,3);
                    tempXu1(j-test_offset(1)) = tempXu1(j-test_offset(1)) ...
                        + Gaussq(T_space_test(j+step),T_space_test(j+step+1),h,3);
                end
            end
            
            tempYu0 = zeros(ntest-sum(test_offset),1);
            tempYu1 = zeros(ntest-sum(test_offset),1);
            for j=1+test_offset(1):ntest-test_offset(2)
                g = @(y) u0_y(y) * Ndiff(T_space_test, splineOrder,0,j,y);
                h = @(y) u1_y(y) * Ndiff(T_space_test, splineOrder,0,j,y);
                
                for step = 0:splineOrder-1
                    tempYu0(j-test_offset(1)) = tempYu0(j-test_offset(1)) ...
                        + Gaussq(T_space_test(j+step),T_space_test(j+step+1),g,3);
                    tempYu1(j-test_offset(1)) = tempYu1(j-test_offset(1)) ...
                        + Gaussq(T_space_test(j+step),T_space_test(j+step+1),h,3);
                end
            end
            U0 = kron(tempXu0, tempYu0);
            U1 = kron(tempXu1, tempYu1);
        else
            U0 = zeros((ntest-sum(test_offset)).^2,1);
            U1 = zeros((ntest-sum(test_offset)).^2,1);
        end
        
        
        
         %% Time stepping
        fprintf('Time stepping...')
        tic
        U = zeros(size(U0,1), K);
        % Save U0
        U(:,1) = pcg(M, U0, ...
            1e-7, ... % Tolerance
            min(5000, size(U0,1))); % maxit
        
        % Calculate the first step:
        [temp, flag] = pcg(M, U1, ...
            1e-7, ... % Tolerance
            min(5000, size(U0,1))); % maxit
        
        if flag ~= 0
            warning('CG did not converge')
        end
        
        [temp2, flag] = pcg(M, Fvals(:,1)-A*U(:,1), ...
            1e-7, ... % Tolerance
            min(5000, size(U0,1))); % maxit
        
        if flag ~= 0
            warning('CG did not converge')
        end
        
        U(:,2) = U(:,1) + tau*temp  + (tau*tau/2) * temp2;
        ts4 = tau*tau/4;
        funMA = @(x) M*x + ts4 * A*x;
        
        % Calculate the other steps:
        for k=3:K
            disp(num2str(k))
            Fpp = Fvals(:,k-2); % penultimate time step
            Fp = Fvals(:,k-1); % last time step
            F = Fvals(:,k); % current time step
            
            % Crank-Nicolson
            %             U(:, k) = (M + ts4 * A) \  ( ts4* (F + 2*Fp+ Fpp) ...
            %                 - ts4 * A * (2*U(:, k-1) + U(:, k-2)) ...
            %                 + M * (2*U(:,k-1) - U(:, k-2)));
            
            [U(:, k), flag] = pcg(funMA,  ts4* (F + 2*Fp+ Fpp) ...
                - ts4 * A*(2*U(:, k-1) + U(:, k-2)) ...
                + M*(2*U(:,k-1) - U(:, k-2)),...
                1e-7, ... % Tolerance
                min(5000, size(U0,1)), ... % maxIt
                [], [], ... %Preconditioner
                U(:,k-1)); % Inital value
            
            if flag ~= 0
                warning('CG did not converge')
            end
            if true%mod(k, K/4) == 0
                fprintf('Time stepping progress: %.2f%%\n', 100*k/K)
            end
        end
        fprintf(' Done.\n')
        times(refinement) =  toc;
        
        %% Get the solution
        tic;
        % Add zeros coefficients for the boundaries
        u_full = zeros(sqrt(size(U,1)) + sum(test_offset),sqrt(size(U,1)) + sum(test_offset), size(U,2));
        u_full(1+test_offset(1):end-test_offset(2), ...
            1+test_offset(1):end-test_offset(2),:) = ...
            reshape(U, [sqrt(size(U,1)),sqrt(size(U,1)), size(U,2)]);
        
        splines = zeros(size(u_full,1), length(x));
        splines_2 = zeros(size(u_full,1), length(x));
        
        for i=1:size(splines,1)
            for j=1:size(splines,2)
                g = @(s) Ndiff(T_space_ansatz, splineOrder,0,i,(s==1)*(s-eps) + (s < 1)*s);
                h = @(s) Ndiff(T_space_ansatz, splineOrder,2,i,(s==1)*(s-eps) + (s < 1)*s);
                
                splines(i,j) = g(x(j));
                splines_2(i,j) = h(x(j));
            end
        end
        
        for i=1:size(u_full,1) % every space node in x
            spline_x = splines(i,:);
            
            if laplace
                spline_x_2 = splines_2(i,:);
                
                space_x_indices = min(find(spline_x, 1, 'first'), find(spline_x_2, 1, 'first')):...
                    max(find(spline_x, 1, 'last'), find(spline_x_2, 1, 'last'));
            else
                space_x_indices = find(spline_x, 1, 'first'):find(spline_x, 1, 'last');
            end
            
            for j=1:size(u_full,2) % every space node in y
                spline_y = splines(j,:);
                
                if laplace
                    spline_y_2 = splines_2(j,:);
                    
                    space_y_indices = min(find(spline_y, 1, 'first'), find(spline_y_2, 1, 'first')):...
                        max(find(spline_y, 1, 'last'), find(spline_y_2, 1, 'last'));
                else
                    space_y_indices = find(spline_y, 1, 'first'):find(spline_y, 1, 'last');
                end
                
                
                for k=1:K % Every time step     
                
                    if u_full(i,j,k) < eps && u_full(i,j,k) > - eps
                        continue;
                    end
                    
                    
                    if laplace
                        spline = kron(spline_y(space_y_indices), spline_x_2(space_x_indices)) ...
                            + kron(spline_y_2(space_y_indices), spline_x(space_x_indices));
                    else
                        spline = kron(spline_y(space_y_indices), spline_x(space_x_indices));
                    end
                    
                    spline = reshape(spline, ...
                        [length(space_x_indices), length(space_y_indices)]);
                    sol(space_x_indices, space_y_indices, k) = ...
                        sol(space_x_indices, space_y_indices, k) + ...
                        u_full(i,j,k) * spline;
                    
                end
            end
        end
        
        toc
        
        
        % Linearaly interpolate the solution to match Kref time steps
        [X, Y, T] = ndgrid(x, y, linspace(0,1, K));
        [Xq, Yq, Tq] = ndgrid(x, y, linspace(0,1, Kref));
        solInt = interpn(X,Y,T,sol,Xq,Yq,Tq);
       
        l2error(refinement) = sqrt(mean((solInt-sol_ana).^2, 'all'));
        
        t = table(refinements, times, l2error)
            writetable(t, 'data/2Dexample1-Timestepping', 'Delimiter',' ');
end

