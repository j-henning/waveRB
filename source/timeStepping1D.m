clear
close all
clc

addpath('../Toolkit-Splines');
addpath('../Toolkit-Splines/Utilities');
addpath('../Solvers');





%% Define the example
% Example 1
% f_time = {@(t) 2+0*t, @(t) t.^2};
% f_space = {@(x) sin(2*pi*x), @(x) 4*pi^2*sin(2*pi*x)};
% u0 = @(x) 0.*x;
% u1 = @(x) 0.*x;
% u_analytical = @(x,t) t.^2 .* sin(2*pi*x);
% dAlemebert = false;

% Example 2
f_time = {@(t) 0*t};
f_space = {@(x) 0*x};
u0 = @(x) x .* (x < 0.5) + (1-x) .* (x >= 0.5);
u1 = @(x) 0.*x;
u_analytical = @(x,t) 0;
dAlemebert = true;

% Example 3
% f_time = {@(t) 0*t};
% f_space = {@(x) 0*x};
% u0 = @(x) (x > 0.25 &&  x < 0.75);
% u1 = @(x) 0.*x;
% u_analytical = @(x,t) 0;
% dAlemebert = true;


% Example 3
% ctime = 5;
% cspace = 5;
% f_time = {@(t) (2-ctime^2*t^2).*sin(ctime*t)+4*ctime*t.*cos(ctime*t), ...
%     @(t) t.^2.*sin(ctime*t)};
% f_space = {@(x) sin(cspace*pi*x), @(x) cspace^2*pi^2*sin(cspace*pi*x)};
% u0 = @(x) 0.*x;
% u1 = @(x) 0.*x;
% u_analytical = @(x,t) t.^2 .* sin(ctime*t).* sin(cspace*pi*x);
% dAlemebert = false;

% Example 3
% u0 = @(x) (x > 0.25) .* (x < 0.75);% sin(2*pi*x);
% u1 = @(x) 0.*x;
% dAlemebert = true;


%% Define the ansatz and trial splines
splineOrder = 3;
laplace = 0; % 0 =  trial equals test space, 2 = trial equals laplace(test)

%% Define the resolutions to test
space_refinements = 1:6;
%time_refinements = 1:7;

% Define the resolution for the L2 error calculation (and plotting)

ref_plotting = 11;
x = linspace(0,1,2^ref_plotting+1);
Kref = 2^ref_plotting+1;

if dAlemebert
    sol_ref = dAlembert1D (u0, u1, 1, length(x), 1, Kref);
else
    [X, T] = meshgrid(x, linspace(0,1,Kref));
    sol_ref = u_analytical(X,T)';
end
% sol_ref = solRef;

for refinement_space = space_refinements
    for refinement_time = refinement_space% time_refinements
        fprintf('Space refinement: %d, Time refinement: %d\n', ...
            refinement_space, refinement_time);
        
        K = 2^ refinement_time; % Number of time steps
        tau = 1 / (K-1);
        sol = zeros(length(x), K);
        
        % If the rhs is zero, we can use dAlemebert to get an analytical
        % solution
        
        
        % T_space_ansatz= InitUniNodes(refinementLevel_space, ...
        %                 (bSplineOrder_space_ansatz==0)*bSplineOrder_space_test ...
        %                 + (bSplineOrder_space_ansatz~=0)*bSplineOrder_space_ansatz);
        
        
        T_space_test = InitUniNodes(refinement_space, splineOrder); %refinement, splineOrder);
        T_space_ansatz = InitUniNodes(refinement_space, splineOrder);
        
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
        
        funM = @(x) M * x;
        funA = @(x) A * x;
        
        
        
        %% Calculation of the right hand side
        ntest = length(T_space_test) - splineOrder;
        Fvals = zeros(ntest-2,K);
        
        if ~dAlemebert
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
        U(:,1) = pcg(M,U0,1e-7);
        
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
                + M * (2*U(:,k-1) - U(:, k-2))), 1e-7, 1e7, [], [], U(:, k-1));
            
            if mod(k, K/4) == 0
                fprintf('Time stepping progress: %.2f%%\n', 100*k/K)
            end
        end
        times(refinement_space, refinement_time) =  toc;
        
        
        %% Get the solution
        tic;
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
                
                %                     spline = reshape(spline, ...
                %                         [length(space_x_indices), length(space_y_indices)]);
                sol(space_x_indices, k) = ...
                    sol(space_x_indices, k) +  u_full(i,k) * spline';
            end
        end
        toc
        % Linearaly interpolate the solution to match Kref time steps
        [X, Y] = meshgrid(x, linspace(0,1, K));
        [Xq, Yq] = meshgrid(x, linspace(0,1, Kref));
        solInt = interp2(X,Y,sol',Xq,Yq)';
        %                 s = surf(solInt);
        %                 s.EdgeAlpha = 0;
        %                 return
        %
        
        %         surf(sol, 'FaceAlpha',0), hold on
        %         s = surf(sol_dAlembert, 'FaceAlpha',0.7);
        %         s.EdgeAlpha = 0;
        %         return
        
        
        l2error(refinement_space, refinement_time) = sqrt(mean((solInt-sol_ref).^2, 'all'));
    end
end

% semilogy(space_refinements, l2error(space_refinements,:), 'o--', 'LineWidth', 2)
% xlabel('Space refinement', 'Interpreter', 'Latex', 'FontSize', 20)
% ylabel('$L_2$ error', 'Interpreter', 'Latex', 'FontSize', 20)
% title('1D Laplace case',  'Interpreter', 'Latex', 'FontSize', 20)
% grid on
% legend('Time refinement = 1' , 'Time refinement = 2', 'Time refinement = 3',...
%     'Time refinement = 4','Time refinement = 5','Time refinement = 6', ...
%     'Time refinement = 7','Time refinement = 8','Time refinement = 9')


% figure
% s = surf(solInt, 'FaceAlpha',.5); hold on
% s.FaceColor = 'red';
% s.EdgeAlpha = 0;
% 
% s = surf(sol_ref, 'FaceAlpha',0.5);
% s.FaceColor = 'blue';
% s.EdgeAlpha = 0;

%% Write the whole data into files
% space_refinement = space_refinements';
% writetable(table(space_refinement,times, l2error), ...
%            '../data/1Dexample2-timestepping','Delimiter',' ')

% refinement = space_refinements';
% time = diag(times);
% error = diag(l2error);
% writetable(table(refinement,time, error), ...
%            '../data/1Dexample3-timestepping','Delimiter',' ')



