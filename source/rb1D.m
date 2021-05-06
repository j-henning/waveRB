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
problemConfiguration.u_0_x = @(x) x* (x < 0.5) + (1-x).*(x > 0.5);
problemConfiguration.u_1_x = @(x) 0*x;
problemConfiguration.inital_conditions = true;
problemConfiguration.has_analytical_solution = false;


% Example 3
% problemConfiguration.f_time = {@(t) 0*t};
% problemConfiguration.f_space = {@(x) 0*x};
% problemConfiguration.u_0_x = @(x) (x > 0.25 &&  x < 0.75);
% problemConfiguration.u_1_x = @(x) 0*x;
% problemConfiguration.inital_conditions = true;
% problemConfiguration.has_analytical_solution = false;

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




problemConfiguration.refinementLevel_time = 4;
problemConfiguration.refinementLevel_space = 4;

resolution.x = 5;
resolution.y = 5;
resolution.t = 5;


problemConfiguration.mu = 1;

N_max = 15;
tolerance = 1e-3;
Xi = linspace(.5,3,50);




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
    problemConfiguration.mu = mu;
    
    
    problem = pOne;
    problem.mu = mu;
    problem.A_space = mu * problem.A_space;
    problem.Q_space = mu^2 * problem.Q_space;
    
    switch problemConfiguration.d
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
S = Xi(ceil(length(Xi)/2));
Y = U(:,ceil(length(Xi)/2));

fprintf("N = %d (%d)\n", 1, N_max)
fprintf("Choosing %f as inital mu\n", S(1))

% X = sol(:,1);




for i=1:N_max-1
    fprintf("N = %d (%d)\n", i+1, N_max)
    % TODO: break if a tolerance is reached
    
    %     % Calculate the reduced solution for all mu and calculate the error
    %     [A_Nq, f_Nq] = projectSystem(A_hq, f_hq, X, Y);
    for j=1:length(Xi)
        
        mu = Xi(j);
        if ismember(mu, S)
            error(j) = 0;
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
        
        
        B_N = Y' * (kron(p.Q_time, p.M_space) ...
            + kron(p.D_time, p.A_space') ...
            + kron(p.D_time', p.A_space) ...
            + kron(p.M_time, p.Q_space)) * Y;
        f_N = Y' * p.rhs(:); % X or Y?
        
        u_N = B_N \ f_N;
        
        u_N_rec = Y * u_N;
        
        p.mu = Xi(j);
        
        switch problemConfiguration.d
            case 1
                sol_N = get1Dsolution(p, u_N_rec,  resolution);
            case 2
                sol_N = get2Dsolution(p, u_N_rec,  resolution);
        end
        
        error(j) = sqrt(mean( (sol{j}-sol_N).^2, 'all'));
        
        
        res(j) = residuum(problem, mu, u_N_rec);
        
    end
    [maxError(i), index] = max(error);
    maxRes(i) = max(res);
    if maxError(i) < tolerance
        fprintf('Reached desired tolerance!\n')
        break;
    end
    
    fprintf("Choosing %f as next mu for the basis (Error: %f)\n", ...
        Xi(index), maxError(i));
    S = [S, Xi(index)];
    Y = [Y, U(:,index)];
    
    
end
toc

%% Online phase


%% Check error for all parameters

for j=1:length(Xi)
    
    mu = Xi(j);
    
    p = pOne;
    p.A_space = mu * p.A_space;
    p.Q_space = mu^2 * p.Q_space;
    
    % Compute X
    
    % X = Y
    %          B_N = X' * (B) * Y
    B_N = Y' * (kron(p.Q_time, p.M_space) ...
        + kron(p.D_time, p.A_space') ...
        + kron(p.D_time', p.A_space) ...
        + kron(p.M_time, p.Q_space)) * Y;
    f_N = Y' * p.rhs(:); % or Y?
    
    u_N = B_N \ f_N;
    
    u_N_rec = Y * u_N;
    
    p.mu = mu;
    
    switch problemConfiguration.d
        case 1
            sol_N = get1Dsolution(p, u_N_rec,  resolution);
        case 2
            sol_N = get2Dsolution(p, u_N_rec,  resolution);
    end
    
    errorTest(j) = sqrt(mean( (sol{j}-sol_N).^2, 'all'));
    
end

subplot(1,2,1)
semilogy(Xi, errorTest, '*--'), hold on, grid on
semilogy(S, 0.5*min(errorTest)*ones(size(S)), '*')
xlabel('\mu')
ylabel('Error')

subplot(1,2,2)
semilogy(maxError, '*--'), hold on, grid on
semilogy(maxRes, 'o:')
semilogy([1 length(maxError)], [tolerance tolerance], 'k')
xlabel('N')
ylabel('Error')



% Compute the residduum as an error estimator
% Input:
%   - problem
%   - parameter mu
%   - RB solution u_N_rec(mu) = Y * u_N(mu)
function res = residuum(p, mu, u_N_rec)
res = 0;

if p.d ~= 1
    error('Not yet implemented')
end

rhs = p.rhs; % No dependency in mu


B = kron(p.Q_time, p.M_space) ...
        + kron(p.D_time, p.A_space') ...
        + kron(p.D_time', p.A_space) ...
        + kron(p.M_time, p.Q_space);
    
    
rhs = rhs(:);

for i =1:size(B,1)
    r = B(i,:) * u_N_rec - rhs(i);
    res = max(r,res);
end

% for i = 1:size(rhs,1)
%     for j = 1:size(rhs,2)
%     % Compute f(v_h, mu) = f(v_h)
%     
%     f = rhs(i,j);
%     % Compute a(u_N(mu), v; mu)
%     a = 0;
%     % Compute ||v_h||_{V_h}
%     normV = 0;
%     
%     res = max(res, abs(f - a) / normV);
%     end
% end


end





function U = solveProblem(problem)


% Use Galerkin for solving
[uu,ss,vv]=svds(problem.rhs,1);
rhs1=uu(:,1)*sqrt(ss(1,1));
rhs2=vv(:,1)*sqrt(ss(1,1));
tolG=1e-5;
maxIt = 100;
info = 0;
[X1,X2]= ...
    Galerkin3(problem.M_space,2*problem.A_space,problem.Q_space,problem.Q_time,...
    (problem.D_time+problem.D_time')/2,problem.M_time,rhs1,rhs2,maxIt,tolG,info);
U=X1*X2';
U=U(:);

% resolution.x = 8;
% resolution.t = 8;
% sol = get1Dsolution(problem, U, resolution);

end



