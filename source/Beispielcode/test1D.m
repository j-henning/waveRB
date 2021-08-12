clear
close all
clc

addpath('splines');
addpath('splines/Utilities');

%% Specify all needed input data
problemConfiguration.bSplineOrder_time = 3;
problemConfiguration.bSplineOrder_space = 3;

problemConfiguration.offset_time_ansatz = [0, 2];
problemConfiguration.offset_time_test = [0, 2];

problemConfiguration.offset_space_ansatz = [1, 1];
problemConfiguration.offset_space_test = [1, 1];

problemConfiguration.d = 1; % Dimension


problemConfiguration.f_time = {@(t) 0*t};
problemConfiguration.f_space = {@(x) 0*x};
problemConfiguration.u_0_x = @(x)  1 * (x > 0.25 &&  x < 0.75);
problemConfiguration.u_1_x = @(x) 0*x;
problemConfiguration.inital_conditions = true;

%% Specify the resolution for plotting and the error calculation

ref_plotting = 11;
x = linspace(0,1,2^ref_plotting+1);
t = linspace(0,1,2^ref_plotting+1);

resolution.x = ref_plotting;
resolution.t = ref_plotting;


l2error = true; % Calculate the l2 error

tolerance = 1e-10;
maxIt = 100;
exactFlag = true;


%% Compute the analytical solution
sol_ana = dAlembert1D(problemConfiguration.u_0_x,...
    problemConfiguration.u_1_x, 1, length(x), 1, length(t));

%% Compute the numerical solution
problemConfiguration.refinementLevel_space = 6;
problemConfiguration.refinementLevel_time = 6;

fprintf('Creating problem...')
problem = create1DWaveProblem(problemConfiguration);
fprintf(' Done!\n')

% Galerkin
[uu,ss,vv]=svds(problem.rhs,1); % Problem.rhs is a matrix
rhs1=uu(:,1)*sqrt(ss(1,1));
rhs2=vv(:,1)*sqrt(ss(1,1));
tolG=tolerance;
info=1;
[X1,X2,restot,~]= Galerkin3(problem.M_space,2*problem.A_space,...
    problem.Q_space,problem.Q_time,(problem.D_time+problem.D_time')/2,...
    problem.M_time,rhs1,rhs2,maxIt,tolG,1e-2,info);
U=X1*X2';
U_galerkin=U(:);


% Get the solution
fprintf('Getting the galerkin solution\n');
solGalerkin = get1Dsolution(problem, U_galerkin, resolution);

% Plotting
figure
subplot(1,2,1)
s = surf(sol_ana); s.EdgeAlpha = 0;
title('Analytical solution')

subplot(1,2,2)
s = surf(solGalerkin); s.EdgeAlpha = 0;
title('Numerical solution')
