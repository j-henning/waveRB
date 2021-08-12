clear
close all
clc

format short

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
problemConfiguration.u_0_x = @(x) 100 * (x > 0.25 && x < 0.75) ...
    .* ((x-0.25).^2 .* (x -0.75).^2);
problemConfiguration.u_1_x = @(x) 0*x;
problemConfiguration.inital_conditions = true;
problemConfiguration.has_analytical_solution = false;

% Example 3
% problemConfiguration.d = 1; % Dimension
% problemConfiguration.f_time = {@(t) 0*t};
% problemConfiguration.f_space = {@(x) 0*x};
% problemConfiguration.u_0_x = @(x) (x > 0.25 && x < 0.75);
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




problemConfiguration.refinementLevel_time = 5;
problemConfiguration.refinementLevel_space = 5;

resolution.x = 7;
resolution.y = 7;
resolution.t = 7;


problemConfiguration.mu = 1;

maxIt = 100;
tolerance = 1e-4;
toleranceEst = 1e-10; % Error estimator

muMin = 1;
muMax = 1.2;

Xi = linspace(muMin,muMax,100);
XiTest = rand(1,100) * (muMax - muMin) + muMin;
XiTest = sort(XiTest, 'ascend');

Y_save = cell(maxIt,1); % save all Y_N (due to orthogonalization)



%% Offline phase

% Calculate the fine solution for all Xis
% (to be replaced with a surrogate)



problemConfiguration.mu = 1;

switch problemConfiguration.d
    case 1
        pOne = create1DWaveProblem(problemConfiguration);
    case 2
        pOne = create2DWaveProblem(problemConfiguration);
end


U = zeros(2^problemConfiguration.refinementLevel_time * ...
    2^problemConfiguration.refinementLevel_space, length(Xi));
sol = cell(length(Xi),1);


% Precompute the needed splines
[splines_time, splines_2_time, splines_space, splines_2_space] = ...
    precompute1DSplines(pOne,  resolution);

fprintf("Calculating the solutions for mu = %e, ..., %e", Xi(1), Xi(end));
parfor i=1:length(Xi)
    mu = Xi(i);
    
    pC = problemConfiguration;
    pC.mu = mu;
    
    
    problem = pOne;
    problem.mu = mu;
    problem.A_space = mu * problem.A_space;
    problem.Q_space = mu^2 * problem.Q_space;
    
    switch pC.d
        case 1
            U(:,i) = solveProblem(problem);
            sol{i} = get1DsolutionSmart(problem, U(:,i),  resolution,splines_time, splines_2_time, splines_space, splines_2_space);
        case 2
            U(:,i) = solveProblem(problem);
            sol{i} = get2Dsolution(problem, U(:,i),  resolution);
    end
    
end
fprintf(' Done!\n')

tDim = size(sol{1},1); % Maybe the dimensions are mixed up
xDim = size(sol{1},2); % Right now, the solution is quadratic
[X, T] = meshgrid(linspace(0,1,xDim), linspace(0,1,tDim));

% Plot some solutions
f = figure;
f.WindowState = 'maximized'; % Requires R2018a
subplot(3,3,1)
[~,h] = contourf(X, T, sol{1},50);
set(h,'LineColor','none')
xlabel('Time')
ylabel('Space')
title(['Solution for mu = ' num2str(Xi(1))])

subplot(3,3,2)
[~,h] = contourf(X, T, sol{ceil(length(Xi)/2)},50);
set(h,'LineColor','none')
xlabel('Time')
ylabel('Space')
title(['Solution for mu = ' num2str(Xi(ceil(length(Xi)/2)))])

subplot(3,3,3)
[~,h] = contourf(X, T, sol{end},50);
set(h,'LineColor','none')
xlabel('Time')
ylabel('Space')
title(['Solution for mu = ' num2str(Xi(end))])
view(2)
drawnow

% Pick a first mu
S = Xi(ceil(length(Xi)/2));
Y = U(:,ceil(length(Xi)/2)) / norm(U(:,ceil(length(Xi)/2)));
fprintf("N = %d (%d)\n", 1, maxIt)
fprintf("Choosing %f as inital mu\n", S(1))

Q_tM_h = kron(pOne.Q_time, pOne.M_space);

%% Greedy algorithm

for i=1:maxIt-2
    tic;
    fprintf("N = %d\n", i+1)
    
    %     % Create all four sub matrices for B
    %
    %     B11 = zeros(i,i); % dt^2, dt^2
    %     B12 = zeros(i,i); % dt^2, dx^2
    %     B21 = zeros(i,i); % dx^2, dt^2
    %     B22 = zeros(i,i); % dx^2, dx^2
    %
    %     % j and l run over all i^2 combinations
    %     % m and n run over all basis functions
    %     for j = 1:i
    %         for l = 1:i
    %             for m = 1:size(U,1)
    %                 for n = 1:size(U,1)
    %                    B11(j,l) = B11(j,l) + 0;
    %                    B12(j,l) = B12(j,l) + 0;
    %                    B21(j,l) = B21(j,l) + 0;
    %                    B21(j,l) = B21(j,l) + 0;
    %                 end
    %             end
    %         end
    %     end
    
    
    parfor j=1:length(Xi)
        
        
        
        mu = Xi(j);
        
        % Check if we already have mu in our set. In this case, continue
        % with the next parameter mu
        if ismember(mu, S)
            error(j) = 0;
            res(j) = 0;
            continue;
        end
        
        % Create a problem for the current paramter
        
        p = pOne;
        p.mu = mu;
        p.A_space = mu * p.A_space;
        p.Q_space = mu^2 * p.Q_space;
        
        
        % Compute the reduced system
        
        B_N = Y' * (kron(p.Q_time, p.M_space) ...
            + kron(p.D_time, p.A_space') ...
            + kron(p.D_time', p.A_space) ...
            + kron(p.M_time, p.Q_space)) * Y;
        f_N = Y' * p.rhs(:);
        
        % Compute the reduces solution u_N_rec
        u_N = B_N \ f_N;
        
        u_N_rec = Y * u_N;
        
        switch problemConfiguration.d
            case 1
                sol_rec = get1DsolutionSmart(p, u_N_rec,  resolution,splines_time, splines_2_time, splines_space, splines_2_space);
            case 2
                sol_rec = get2Dsolution(p, u_N_rec,  resolution);
        end
        
        error(j) = sqrt(mean( (sol{j}-sol_rec).^2, 'all'));
        
        % Compute an error estimator to compute the residuum
        % res(j) = residuum(p, u_N_rec);
        
    end
    % Check for which mu the error was the largest
    [maxError(i), index] = max(error);
    %     [maxRes(i), indexRes] = max(res);
    
    iterTime = toc;
    fprintf('This iteration took %fs\n', iterTime)
    
    
    
    
    
    
    %      Unew = U(:,index) - Y * Y' * U(:,index);
    %      Unew = Unew / norm(Unew);
    %      Y = [Y, Unew];
    
    
    % Draw the current progress
    subplot(3,3,4), hold off
    semilogy(maxError, '*-'), hold on, grid on
    semilogy([1 i], [tolerance tolerance], 'k--')
    semilogy([1 i], [toleranceEst toleranceEst], 'k-')
    ylim([0.1 * toleranceEst, max(maxError)])
    xlabel('N')
    ylabel('Error')
    legend('Max Error', 'Tolerance', 'Tolerance Error Estimator')
    title('Error over iterations')
    
    subplot(3,3,5), hold off
    semilogy(Xi, error, 'o-'), grid on
    xlabel('mu')
    ylabel('Error')
    legend('Error')
    title('Error over training set')
    
    subplot(3,3,6)
    histogram(S,10)
    xlabel('mu')
    ylabel('#Samples')
    title('Sample distribution')
    
    
    drawnow
    
    % add the found parameter to our basis
    fprintf("Choosing mu = %f as next mu for the basis (Error: %e)\n", ...
        Xi(index), maxError(i));
    
    
    
    % Modiefied Gram Schmidt (does not seem to work):
    q = U(:,index);
    %     for j=1:i
    %         q = q - (q'*Y(:,i)) * q;
    %     end
    %     q = q / norm(q);
    %
    S = [S, Xi(index)];
    Y = [Y, q];
    
    Y = orth(Y); % Could be done more efficient
    Y_save{i} = Y;
    
    
    % If the error is smaller than the tolerance, we are done
    if maxError(i) < toleranceEst
        fprintf('Reached desired tolerance!\n')
        break;
    end
end







fprintf('Starting online phase\n')
solTest = cell(length(XiTest),1);

% Hierachical error estimator
index_N = find(maxError > tolerance, 1, 'last') + 1;

if index_N == length(maxError) + 1 || index_N == 0
    warning('Did not reach desired tolerance. Using different splitting')
    
    tol_Alternative = 0.5 * (maxError(1) + maxError(end));
    index_N = find(maxError > tol_Alternative, 1, 'last') + 1;
end
S_N = S(1:index_N);
% Y_N = orth(Y(:,1:index_N));
Y_N = Y_save{index_N};

S_M = S;
Y_M = Y;

%% Calculation of the theta value for the online error estimator




minProb = @(x) minProblem(x, pOne, Y_N, Y_M,resolution,splines_time, splines_2_time, splines_space, splines_2_space);


% Split the interval from muMin to muMax into several subintervals and find
% for each one the best theta

theta = 0;

muInt = linspace(muMin, muMax, 20);
thetaInt = linspace(0, 2, 20);
minVal = Inf;

for i=1:length(muInt)
    minVal = Inf;
    for j = 1:length(thetaInt)
        val = minProb([muInt(i), thetaInt(j)]);
        if val < minVal
            minVal = val;
            thetaLocal = thetaInt(j);
        end    
    end
    thetaVec(i) = thetaLocal;
    
    if thetaLocal > 1
        warning('Theta is bigger than 1 for mu = %f!', muInt(i))
    end
end

theta = max(thetaVec);

% for i=1:length(muInt) - 1
%     [~,fval,exitflag,output] = fmincon(minProb,[(muInt(i) + muInt(i+1)) / 2, 0],[],[],[],[],[muInt(i)-eps, 0-eps], [muInt(i+1)+eps, 2+eps]);
%     x = output.bestfeasible.x;
%     theta = max(theta, x(2));
%
%     if theta > 1
%         warning('Theta is larger than 1! Stopping theta calculation!')
%         break;
%     end
% end


fprintf('Theta = %f\n', theta)




%% Online phase: Check error and hierarchical error estimator for different
% parameters


fprintf("Calculating the test solutions for mu = %e, ..., %e", XiTest(1), XiTest(end));
for i=1:length(XiTest)
    mu = XiTest(i);
    
    pC = problemConfiguration;
    pC.mu = mu;
    
    
    problem = pOne;
    problem.mu = mu;
    problem.A_space = mu * problem.A_space;
    problem.Q_space = mu^2 * problem.Q_space;
    
    switch pC.d
        case 1
            Utest = solveProblem(problem);
            solTest{i} = get1DsolutionSmart(problem, Utest,  resolution,splines_time, splines_2_time, splines_space, splines_2_space);
        case 2
            Utest = solveProblem(problem);
            solTest{i} = get2Dsolution(problem, Utest,  resolution);
    end
    
end
fprintf(' Done!\n')



fprintf('Testing all solutions')
for j=1:length(XiTest)
    
    mu = XiTest(j);
    
    p = pOne;
    p.A_space = mu * p.A_space;
    p.Q_space = mu^2 * p.Q_space;
    p.mu = mu;
    
    % Compute X
    
    % X = Y
    %          B_N = X' * (B) * Y
    
    % Calculte the reduced solution
    B_N = Y_N' * (kron(p.Q_time, p.M_space) ...
        + kron(p.D_time, p.A_space') ...
        + kron(p.D_time', p.A_space) ...
        + kron(p.M_time, p.Q_space)) * Y_N;
    f_N = Y_N' * p.rhs(:);
    
    u_N = B_N \ f_N;
    
    u_N_rec = Y_N * u_N;
    
    % Calculate the reduced solution of higher dimension (for the error
    % estimation)
    
    B_M = Y_M' * (kron(p.Q_time, p.M_space) ...
        + kron(p.D_time, p.A_space') ...
        + kron(p.D_time', p.A_space) ...
        + kron(p.M_time, p.Q_space)) * Y_M;
    f_M = Y_M' * p.rhs(:); % or Y?
    
    u_M = B_M \ f_M;
    
    u_M_rec = Y_M * u_M;
    
    switch problemConfiguration.d
        case 1
            sol_rec_N = get1DsolutionSmart(p, u_N_rec,  resolution,splines_time, splines_2_time, splines_space, splines_2_space);
            sol_rec_M = get1DsolutionSmart(p, u_M_rec,  resolution,splines_time, splines_2_time, splines_space, splines_2_space);
        case 2
            sol_rec_N = get2Dsolution(p, u_N_rec,  resolution);
            sol_rec_M = get2Dsolution(p, u_M_rec,  resolution);
    end
    
    
    
    errorTest(j) = sqrt(mean( (solTest{j}-sol_rec_N).^2, 'all'));
    errorEstTest(j) = sqrt(mean( (sol_rec_M-sol_rec_N).^2, 'all'));
end
fprintf(' Done!\n')
fprintf('Desired tolerance : %e\n', tolerance)
fprintf('Maximum test error: %e\n', max(errorTest))

errorEstTest = errorEstTest./(1-theta);

% Plotting

subplot(3,3,7:8)
semilogy(XiTest, errorTest, '*-'), hold on, grid on
semilogy(XiTest, errorEstTest, 'd:')
semilogy(S, 0.5*min(errorTest)*ones(size(S)), 'k*')
xlabel('\mu')
ylabel('Error')
legend('Error', 'Error estimator', 'Snapshot parameters')
title(['Error over test set (theta = ' num2str(theta) ')'])

subplot(3,3,9)
plot(XiTest,errorEstTest-errorTest), hold on
plot(XiTest(find(errorEstTest-errorTest < 0)), ...
    errorEstTest(find(errorEstTest-errorTest < 0)) ...
    -errorTest(find(errorEstTest-errorTest < 0)), 'ro')
title('Error Estimator - Error')
legend('Estimator - Error', 'Values smaller than 0')
grid on








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
    (problem.D_time+problem.D_time')/2,problem.M_time,rhs1,rhs2,maxIt,tolG,1e-3,info);
U=X1*X2';
U=U(:);

% resolution.x = 8;
% resolution.t = 8;
% sol = get1Dsolution(problem, U, resolution);

end

% Create a function for the minimization problem
function val = minProblem(x, pOne, Y_N, Y_M,resolution,splines_time, splines_2_time, splines_space, splines_2_space)
mu = x(1);
q = x(2);

% fprintf('Testing mu = %f, q = %f', mu, q);

p = pOne;
p.mu = mu;
p.A_space = mu * p.A_space;
p.Q_space = mu^2 * p.Q_space;


% Compute the reduced system

B_N = Y_N' * (kron(p.Q_time, p.M_space) ...
    + kron(p.D_time, p.A_space') ...
    + kron(p.D_time', p.A_space) ...
    + kron(p.M_time, p.Q_space)) * Y_N;
f_N = Y_N' * p.rhs(:);
u_N = B_N \ f_N;
u_N_rec = Y_N * u_N;

B_M = Y_M' * (kron(p.Q_time, p.M_space) ...
    + kron(p.D_time, p.A_space') ...
    + kron(p.D_time', p.A_space) ...
    + kron(p.M_time, p.Q_space)) * Y_M;
f_M = Y_M' * p.rhs(:);
u_M = B_M \ f_M;
u_M_rec = Y_M * u_M;


sol_rec_N = get1DsolutionSmart(p, u_N_rec,  resolution,splines_time, splines_2_time, splines_space, splines_2_space);
sol_rec_M = get1DsolutionSmart(p, u_M_rec,  resolution,splines_time, splines_2_time, splines_space, splines_2_space);


% Compute the reference solutions

problem = pOne;
problem.mu = mu;
problem.A_space = mu * problem.A_space;
problem.Q_space = mu^2 * problem.Q_space;


U = solveProblem(problem);
sol = get1DsolutionSmart(problem, U,  resolution,splines_time, splines_2_time, splines_space, splines_2_space);


fmu = sqrt(mean( (sol-sol_rec_M).^2, 'all'));
gmu = sqrt(mean( (sol-sol_rec_N).^2, 'all'));

val = abs(fmu - q * gmu);
% fprintf(' - val = %e\n', val)
end


