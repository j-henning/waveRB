function [sol, solvingTime] = timeStepping1D(problemConfiguration, K, space_refinement, resolution, splineOrder, tolerance)



%% Read all need values from the problem struct
f_time = problemConfiguration.f_time;
f_space = problemConfiguration.f_space;
u0 = problemConfiguration.u_0_x;
u1 = problemConfiguration.u_1_x;


laplace = 0; % 0 =  trial equals test space, 2 = trial equals laplace(test)
% Note: Use Laplace with causion!

x = linspace(0,1,2^resolution.x+1);
Kref = 2^resolution.t+1;

tau = 1 / (K-1);
sol = zeros(length(x), K);

T_space_test = InitUniNodes(space_refinement, splineOrder);
T_space_ansatz = InitUniNodes(space_refinement, splineOrder);

%% Calculation of the space matrices

M = StiffMat(T_space_ansatz, ... % Ansatz discretization
    T_space_test, ... % Test discretization
    splineOrder, ... % Ansatz spline order
    splineOrder, ... % Test spline order
    [1 1],... % Ansatz offset
    [1 1], ... % Test offset
    [laplace, 0] ... % diff
    );

A = StiffMat(T_space_ansatz, ... % Ansatz discretization
    T_space_test, ... % Test discretization
    splineOrder, ... % Ansatz spline order
    splineOrder, ... % Test spline order
    [1 1],... % Ansatz offset
    [1 1], ... % Test offset
    [1+laplace, 1] ... % diff
    );


%% Calculation of the right hand side
ntest = length(T_space_test) - splineOrder;
Fvals = zeros(ntest-2,K);


for k=1:K
    for i=1:length(f_time)
        temp = zeros(ntest-2,1);
        for j=2:ntest-1 % Use a offset of [1 1]
            g = @(x) f_space{i}(x) * Ndiff(T_space_test, splineOrder,0,j,x);
            for step = 0:splineOrder-1%bSplineOrder_space_test-1
                temp(j-1) = temp(j-1) ...
                    + Gaussq(T_space_test(j+step),T_space_test(j+step+1),g,3);
            end
        end
        Fvals(:, k) = Fvals(:,k) + temp .* f_time{i}((k-1)/(K-1));
    end
end


%% Calculation of u0 and u1
U0 = zeros(ntest-2,1);
U1 = zeros(ntest-2,1);


for j=2:ntest-1 % Use a offset of [1 1]
    g = @(x) u0(x) .* Ndiff(T_space_test, splineOrder,0,j,x);
    h = @(x) u1(x) .* Ndiff(T_space_test, splineOrder,0,j,x);
    for step = 0:splineOrder-1%bSplineOrder_space_test-1
        U0(j-1) = U0(j-1) ...
            + Gaussq(T_space_test(j+step),T_space_test(j+step+1),g,3);
        U1(j-1) = U1(j-1) ...
            + Gaussq(T_space_test(j+step),T_space_test(j+step+1),h,3);
    end
end

%% Time stepping
tic;
U = zeros(size(M,2), K);
% Save U0
U(:,1) = pcg(M,U0,tolerance);

% Calculate the first step:
U(:,2) = U(:,1) + tau* pcg(M, U1, 1e-7) + (tau*tau/2) * pcg(M,(Fvals(:,1)-A*U(:,1)), 1e-7);

ts4 = tau*tau/4;

% Calculate the other steps:
for k=3:K
    Fpp = Fvals(:,k-2); % penultimate time step
    Fp = Fvals(:,k-1); % last time step
    F = Fvals(:,k); % current time step
    
    % Crank-Nicolson
    [U(:, k), ~] = pcg(M + ts4 * A, ( ts4* (F + 2*Fp+ Fpp) ...
        - ts4 * A * (2*U(:, k-1) + U(:, k-2)) ...
        + M * (2*U(:,k-1) - U(:, k-2))), tolerance, 1e7, [], [], U(:, k-1));
    
    if mod(k, K/4) == 0
        fprintf('Time stepping progress: %.2f%%\n', 100*k/K)
    end
end
solvingTime = toc;

%% Get the solution
% Add zeros coefficients for the boundaries
u_full = zeros(size(U,1) + 2, size(U,2));
u_full(2:end-1,:) = U;


% Precalculate the splines (not really needed in 1D, but efficient
% in 2D and 3D)
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
    
    
    
    for k=1:K % Every time step
        if u_full(i,k) < eps && u_full(i,k) > - eps
            continue;
        end
        
        if laplace == 2
            spline = spline_x_2(space_x_indices);
        else
            spline = spline_x(space_x_indices);
        end
        
        sol(space_x_indices, k) = ...
            sol(space_x_indices, k) +  u_full(i,k) * spline';
    end
end

% Linearaly interpolate the solution to match Kref time steps
[X, Y] = meshgrid(x, linspace(0,1, K));
[Xq, Yq] = meshgrid(x, linspace(0,1, Kref));
sol = interp2(X,Y,sol',Xq,Yq)';
end
