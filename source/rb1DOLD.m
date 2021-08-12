clear
close all
clc

addpath('../splines');
addpath('../splines/Utilities');
addpath('../solvers');

%% Specify all needed input data
problemConfiguration.bSplineOrder_time = 3;
problemConfiguration.bSplineOrder_space = 3;

problemConfiguration.offset_time_ansatz = [0, 2];
problemConfiguration.offset_time_test = [0, 2];
problemConfiguration.offset_space_ansatz = [1, 1];
problemConfiguration.offset_space_test = [1, 1];


%% 1D Examples
% Example 1
% problemConfiguration.d = 1; % Dimension
% problemConfiguration.f_time = {@(t) 2+0*t, @(t) t.^2};
% problemConfiguration.f_space = {@(x) sin(2*pi*x), @(x) 4*pi^2*sin(2*pi*x)};
% problemConfiguration.u_0_x = @(x) 0.*x;
% problemConfiguration.u_1_x = @(x) 0.*x;
%
% problemConfiguration.u_analytical = @(x,t) t.^2 .* sin(2*pi*x);
% problemConfiguration.dAlemebert = false;
% problemConfiguration.has_analytical_solution = true;

% Example 2
problemConfiguration.d = 1; % Dimension
problemConfiguration.f_time = {@(t) 0*t};
problemConfiguration.f_space = {@(x) 0*x};
problemConfiguration.u_0_x = @(x) (x > 0.25 && x < 0.75);
problemConfiguration.u_1_x = @(x) 0*x;
problemConfiguration.inital_conditions = true;
problemConfiguration.has_analytical_solution = false;



%% 2D examples
% problemConfiguration.d = 2;
%
% problemConfiguration.f_time = {@(t) 0*t};
% problemConfiguration.f_space_x = {@(x) 0*x};
% problemConfiguration.f_space_y = {@(y) 0*y};
% problemConfiguration.u_0_x = @(x) x* (x < 0.5) + (1-x).*(x > 0.5);
% problemConfiguration.u_0_y = @(y)  y* (y < 0.5) + (1-y).*(y > 0.5);
% problemConfiguration.u_1_x = @(x) 0*x;
% problemConfiguration.u_1_y = @(y) 0 *y;
% problemConfiguration.inital_conditions = true;
% problemConfiguration.has_analytical_solution = false;

% problemConfiguration.f_time = {@(t) 0*t};
% problemConfiguration.f_space_x = {@(x) 0*x};
% problemConfiguration.f_space_y = {@(y) 0*y};
% problemConfiguration.u_0_x = @(x) 1 * (x > 0.25) .* (x < 0.75);
% problemConfiguration.u_0_y = @(y) 1 * (y > 0.25) .* (y < 0.75);
% problemConfiguration.u_1_x = @(x) 0 * x;
% problemConfiguration.u_1_y = @(y) 0 * y;
% problemConfiguration.inital_conditions = true;
% problemConfiguration.has_analytical_solution = false;




problemConfiguration.refinementLevel_time = 5;
problemConfiguration.refinementLevel_space = 5;

resolution.x = 6;
resolution.y = 6;
resolution.t = 6;


problemConfiguration.mu = 1;

N_max = 20;

% N and M are used for the hierarchical error estimator which is currently
% not used. So just choose tolerance_M appropriately
tolerance_N = 1e-3;
tolerance_M = 1e-3;

Xi = linspace(.1,3,100);
XiTest = 0.5 * (Xi(1:end-1) + Xi(2:end));




%% Offline phase

% Calculate the fine solution for all Xis
% (to be replaced with a surrogate)

tic;

U = [];
sol = [];


problemConfiguration.mu = 1;

switch problemConfiguration.d
    case 1
        pOne = create1DWaveProblem(problemConfiguration);
    case 2
        pOne = create2DWaveProblem(problemConfiguration);
end


for i=1:length(Xi)
    mu = Xi(i);
    fprintf("Calculating the solution for mu = %f\n", mu);
    
    pC = problemConfiguration;
    pC.mu = mu;
    
    
    problem = pOne;
    problem.mu = mu;
    problem.A_space = mu * problem.A_space;
    problem.Q_space = mu^2 * problem.Q_space;
    
    switch pC.d
        case 1
            U = [U, solveProblem(problem)];
            sol{i} = get1Dsolution(problem, U(:,i),  resolution);
        case 2
            %             problem = create2DWaveProblem(problemConfiguration);
            
            
            U = [U, solveProblem(problem)];
            sol{i} = get2Dsolution(problem, U(:,i),  resolution);
    end
    
end


% pick the first mu before you enter the loop (how?)

% indices = randperm(length(Xi));
% offset = 2;
% 
% S = Xi(indices(1:offset+1));
% Y = U(:,indices(1:offset+1) 

S = Xi(ceil(length(Xi)/2));
Y = U(:,ceil(length(Xi)/2));
fprintf("N = %d (%d)\n", 1, N_max)
fprintf("Choosing %f as inital mu\n", S(1))
% 
% S = [S, Xi(1)];
% Y = [Y, U(:,1)];
% fprintf("N = %d (%d)\n", 2, N_max)
% fprintf("Choosing %f as second mu\n", S(2))
% 
% S = [S, Xi(end)];
% Y = [Y, U(:,end)];
% fprintf("N = %d (%d)\n", 3, N_max)
% fprintf("Choosing %f as second mu\n", S(end))

% X = sol(:,1);



%% Greedy algorithm
for i=1:N_max-2
    fprintf("N = %d (%d)\n", i+1, N_max)
    % TODO: break if a tolerance is reached
    
    %     % Calculate the reduced solution for all mu and calculate the error
    %     [A_Nq, f_Nq] = projectSystem(A_hq, f_hq, X, Y);
    parfor j=1:length(Xi)
        
        mu = Xi(j);
        if ismember(mu, S) 
            error(j) = 0;
            res(j) = 0;
            continue;
        end
        
        %         % Compute A(mu) and f(mu)
        %         A_N = zeros(size(A_Nq{1}));
        %         for q = 1:nummel(A_Nq)
        %             A_N = A_N + theta_Aq{q}(mu) * A_Nq{q};
        %         end
        %         f_N = zeros(size(f_Nq{1}));
        %         for q = 1:numel(f_Nq)
        %             f_N = f_N + theta_fq{q}(mu) * f_Nq{q};
        %         end
        
        p = pOne;
        p.mu = mu;
        p.A_space = mu * p.A_space;
        p.Q_space = mu^2 * p.Q_space;
        
        
        % Compute X
        
        % X = Y
        %          B_N = X' * (B) * Y
        
        
%         B_M = Y' * (kron(p.Q_time, p.M_space) ...
%             + kron(p.D_time, p.A_space') ...
%             + kron(p.D_time', p.A_space) ...
%             + kron(p.M_time, p.Q_space)) * Y;
%         f_M = Y' * p.rhs(:); % X or Y?
%         
%         u_M = B_M \ f_M;
%         
%         u_M_rec = Y * u_M;
        
        % N = M-1
        
        B_N = Y' * (kron(p.Q_time, p.M_space) ...
            + kron(p.D_time, p.A_space') ...
            + kron(p.D_time', p.A_space) ...
            + kron(p.M_time, p.Q_space)) * Y;
        f_N = Y' * p.rhs(:); % X or Y?
        
        u_N = B_N \ f_N;
        
        u_N_rec = Y * u_N;
        
        
        
        switch problemConfiguration.d
            case 1
                sol_N = get1Dsolution(p, u_N_rec,  resolution);
%                 sol_M = get1Dsolution(p, u_M_rec,  resolution);
            case 2
                sol_N = get2Dsolution(p, u_N_rec,  resolution);
%                 sol_M = get2Dsolution(p, u_M_rec,  resolution);
        end
        
        error(j) = sqrt(mean( (sol{j}-sol_N).^2, 'all'));
        
%         res(j) = sqrt(mean( (sol_M-sol_N).^2, 'all'));
   
%         res(j) = residuum(p, u_N_rec);
        
    end
    [maxError(i), index] = max(error);
%     [maxRes(i), indexRes] = max(res);
%     
    if maxError(i) < tolerance_M
        fprintf('Reached desired tolerance!\n')
        break;
    end
    
    %     fprintf("Choosing mu = %f as next mu for the basis (Error: %f, beta = %f)\n", ...
    %         Xi(indexRes), maxRes(i), maxRes(i) / maxError(i));
    
    fprintf("Choosing mu = %f as next mu for the basis (Error: %e)\n", ...
        Xi(index), maxError(i));
    
    subplot(1,2,1), hold off
    semilogy(maxError, '*-'), hold on, grid on
    xlabel('N')
    ylabel('Error')
    legend('Max Error')
    %     drawnow
    
    subplot(1,2,2), hold off
    semilogy(Xi, error, 'o-'), grid on
    xlabel('mu')
    ylabel('Error')
    legend('Error')
%     hold on
%     semilogy(Xi, res, '*--')
    drawnow

    
    S = [S, Xi(index)];
    Y = [Y, U(:,index)];
    
    
end
toc

% Hierachical error estimator
% index_N = find(maxError > tolerance_N, 1, 'last') + 1;
% S_N = S(1:index_N);
% Y_N = Y(:,1:index_N);
% 
% S_M = S;
% Y_M = Y;

Y_N = Y;



%% Online phase: Check error and hierarchical error estimator for all parameters

for j=1:length(XiTest)
    
    mu = XiTest(j);
    
    p = pOne;
    p.A_space = mu * p.A_space;
    p.Q_space = mu^2 * p.Q_space;
    p.mu = mu;
    
    % Compute X
    
    % X = Y
    %          B_N = X' * (B) * Y
    B_N = Y_N' * (kron(p.Q_time, p.M_space) ...
        + kron(p.D_time, p.A_space') ...
        + kron(p.D_time', p.A_space) ...
        + kron(p.M_time, p.Q_space)) * Y_N;
    f_N = Y_N' * p.rhs(:); % or Y?
    
    u_N = B_N \ f_N;
    
    u_N_rec = Y_N * u_N;
    
    
%     B_M = Y_M' * (kron(p.Q_time, p.M_space) ...
%         + kron(p.D_time, p.A_space') ...
%         + kron(p.D_time', p.A_space) ...
%         + kron(p.M_time, p.Q_space)) * Y_M;
%     f_M = Y_M' * p.rhs(:); % or Y?
%     
%     u_M = B_M \ f_M;
%     
%     u_M_rec = Y_M * u_M;
    
    
    
    
    switch problemConfiguration.d
        case 1
            sol_N = get1Dsolution(p, u_N_rec,  resolution);
%             sol_M = get1Dsolution(p, u_M_rec,  resolution);
        case 2
            sol_N = get2Dsolution(p, u_N_rec,  resolution);
%             sol_M = get2Dsolution(p, u_M_rec,  resolution);
    end
    
    errorTest(j) = sqrt(mean( (sol{j}-sol_N).^2, 'all'));
%     residualTest(j) = sqrt(mean( (sol_M-sol_N).^2, 'all'));
    
end

% Plotting
figure
semilogy(XiTest, errorTest, '*-'), hold on, grid on
% semilogy(XiTest, residualTest, 'o:')
semilogy(S, 0.5*min(errorTest)*ones(size(S)), 'k*')
xlabel('\mu')
ylabel('Error')
% legend('Error', 'Error estimator', 'Snapshot parameters')
legend('Error', 'Snapshot parameters')


%% Help functions

% Compute the residduum as an error estimator
% Input:
%   - problem
%   - parameter mu
%   - RB solution u_N_rec(mu) = Y * u_N(mu)
function res = residuum(p, u_N_rec)
res = 0;

if p.d ~= 1
    error('Not yet implemented')
end

rhs = p.rhs; % No dependency in mu


B = kron(p.Q_time, p.M_space) ...
    + kron(p.D_time, p.A_space') ...
    + kron(p.D_time', p.A_space) ...
    + kron(p.M_time, p.Q_space);

M_ = kron(p.M_time, p.M_space);


rhs = rhs(:);


% res = norm(B * u_N_rec - rhs(:));

for j = 1:size(B,2) %test every function of the test space
    
    r =  B(j,:) * u_N_rec;
    %   r = abs(r - rhs(j)) / sqrt(M_(j,j)); % todo: check this
        r = abs(r - rhs(j)) / sqrt(B(j,j)); % todo: check this
    res = max(r,  res);
end
end





function U = solveProblem(problem)


% Use Galerkin for solving
[uu,ss,vv]=svds(problem.rhs,1);
rhs1=uu(:,1)*sqrt(ss(1,1));
rhs2=vv(:,1)*sqrt(ss(1,1));
tolG=1e-8;
maxIt = 100;
info = 0;
[X1,X2]= ...
    Galerkin3(problem.M_space,2*problem.A_space,problem.Q_space,problem.Q_time,...
    (problem.D_time+problem.D_time')/2,problem.M_time,rhs1,rhs2,maxIt,tolG,1e-1,info);
U=X1*X2';
U=U(:);

% resolution.x = 8;
% resolution.t = 8;
% sol = get1Dsolution(problem, U, resolution);

end



