clear
close all
clc

format shortE


addpath('../splines');
addpath('../splines/Utilities');
addpath('../solvers');

% load('3DreferenceSol')


%% Define the example
% Example 1
% mm = 3; nn = 2; oo = 1;
% f_time = {@(t) 2+ 0.*t; @(t) (t.^2)*(mm^2 + nn^2 + oo^2)*(pi^2)};
% f_space_x = {@(x) sin(mm*pi*x); @(x) sin(mm*pi*x)};
% f_space_y = {@(y) sin(nn*pi*y); @(y) sin(nn*pi*y)};
% f_space_z = {@(z) sin(oo*pi*z); @(z) sin(oo*pi*z)};
% u_analytical = @(x,y,z,t) t.^2 .* sin(mm*pi*x).*sin(nn*pi*y).*sin(oo*pi*z);
% inital_conditions = false;

% Example 2
% f_time = {@(t) 0*t};
% f_space_x = {@(x) 0*x};
% f_space_y = {@(y) 0*y};
% f_space_z = {@(z) 0*z};
% u0_x = @(x) x * (x < 0.5) + (1-x).*(x > 0.5);
% u0_y = @(y) y * (y < 0.5) + (1-y).*(y > 0.5);
% u0_z = @(z) z * (z < 0.5) + (1-z).*(z > 0.5);
% u1_x = @(x) 0*x;
% u1_y = @(y) 0*y;
% u1_z = @(z) 0*z;
% u_analytical = @(x,y,z,t) 0.*t.*x.*y.*z; % Not known
% inital_conditions = true;
% load('..\3D-example2.mat')
% sol_ref = solRef;



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



% Example 4
% f_time = {@(t) 0*t};
% f_space_x = {@(x) 0*x};
% f_space_y = {@(y) 0*y};
% u0_x = @(x) sin(pi*x); u0_y = @(y) sin(pi*y);
% u1_x = @(x) 0*x; u1_y = @(y) 0 *y;
% u_analytical = @(x,y,t) 0.*t.*x.*y; % Not known
% inital_conditions = true;

%% Example 5 (radial)
% u0r = @(r) 1.* (r < 0.2);
f_time = {@(t) 0*t};
f_space_x = {@(x) 0*x};
f_space_y = {@(y) 0*y};

funR = @(x,y,z) sqrt((x-0.5).^2 + (y-0.5).^2 + (z-0.5).^2);  % Returns the radius
u0= @(x,y,z)  (1-5 * funR(x,y,z)) .* ( funR(x,y,z) <= 0.2);
% u0= @(x,y) 1 + 0.* x .*y;
u1 = @(x,y,z)  0.*x.*y.*z;
% u_analytical = @(x,y,t) 0.*t.*x.*y; % Not known
inital_conditions = true;
c = 0.2;



%% Define the ansatz and trial splines
laplace = false;
splineOrder = 3; % at least four if laplace is true, otherwise at least 3

%% Define the resolutions to test
% space_refinements = 1:4;
% time_refinements = 1:6;



% Define the resolution for the L2 error calculation (and plotting)
x = linspace(0,1, 2^6+1);
y = linspace(0,1, 2^6+1);
z = linspace(0,1, 2^6+1);
t = linspace(0,1, 2^6+1);
Kref = 2^6+1;

% Compute the analytical solution
[X, Y, Z,  T] = ndgrid(x, y, z, t);
u0r = @(r) (1-5*abs(r)) .*  (abs(r) <= 0.2);
for i=1:2^6+1
    for j=1:2^6+1
        for k=1:2^6+1
            r = sqrt((0.5 -X(i,j,1))^2 + (0.5 - Y(i,j,1))^2 + (0.5 - Z(i,j,1))^2);
            sol_ref(i,j,k,:) = dAlemenbertSphere(r,t,c, u0r);
        end
    end
end


for refinement_space = 6%space_refinements
    for refinement_time = refinement_space%time_refinements
        fprintf('#######################################\n')
        fprintf('Space refinement: %d, Time refinement: %d\n', ...
            refinement_space, refinement_time);
        fprintf('#######################################\n')
        
        K = 2^ refinement_time; % Number of time steps
        tau = 1 / (K-1);
        sol = zeros(length(x), length(y), length(z), K);
        
        T_space_test = InitUniNodes(refinement_space, splineOrder);
        T_space_ansatz = InitUniNodes(refinement_space, splineOrder);
        
        %% Calculation of the space matrices
        fprintf('Calculating matrices...')
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
            % Provide functions handles for the matrix vector
            % multiplication instead of calculating the large matrices
            funM = @(x) reshape(kron(M1D, M1D) * reshape(x, [dimM^2, dimM]) * O'...
                + kron(O, M1D) * reshape(x, [dimM^2, dimM]) * M1D'...
                + kron(M1D, O) * reshape(x, [dimM^2, dimM]) * M1D', [], 1);
            
            funA = @(x) reshape(kron(M1D, M1D) * reshape(x, [dimM^2, dimM]) * S'... %
                + kron(O, M1D) * reshape(x, [dimM^2, dimM]) * Q'...
                + kron(M1D, O) * reshape(x, [dimM^2, dimM]) * Q'...
                + kron(Q, M1D) * reshape(x, [dimM^2, dimM]) * O'...
                + kron(S, M1D) * reshape(x, [dimM^2, dimM]) * M1D'...
                + kron(Q, O) * reshape(x, [dimM^2, dimM]) * M1D'...
                + kron(M1D, Q) * reshape(x, [dimM^2, dimM]) * O'...
                + kron(O, Q) * reshape(x, [dimM^2, dimM]) * M1D'...
                + kron(M1D, S) * reshape(x, [dimM^2, dimM]) * M1D',...
                [], 1);
            
            
            M = kron(kron(O, M1D),M1D) + kron(kron(M1D, O),M1D) + kron(kron(M1D, M1D),O);
            A = kron(S, kron(M1D, M1D)) + kron(Q, kron(O, M1D)) + kron(Q, kron(M1D, O)) ...
                + kron(O, kron(Q, M1D)) + kron(M1D, kron(S, M1D)) + kron(M1D, kron(Q, O)) ...
                + kron(O, kron(M1D, Q)) + kron(M1D, kron(O, Q)) + kron(M1D, kron(M1D, S));
            
        else
            
            % Provide functions handles for the matrix vector
            % multiplication instead of calculating the large matrices
            funM = @(x) reshape(kron(M1D, M1D) * reshape(x, [dimM^2, dimM]) * M1D', [], 1);
            funA = @(x) reshape(kron(M1D, M1D) * reshape(x, [dimM^2, dimM]) * Q' ...
                + kron(Q, M1D) * reshape(x, [dimM^2, dimM]) * M1D' ...
                + kron(M1D, Q) * reshape(x, [dimM^2, dimM]) * M1D', [], 1);
            
            M = kron(M1D, kron(M1D, M1D));
            A = kron(Q, kron(M1D, M1D)) + kron(M1D, kron(Q, M1D)) ...
                + kron(M1D, kron(M1D, Q));
        end
        
        fprintf(' Done.\n')
        
        %   Uncomment if you want to check the eigenvalues
        %                         figure
        %                         subplot(1,3,1)
        %                         E = eig(full(M));
        %                         plot(real(E), imag(E), '*'), grid on
        %                         title('Eigenvalues of M')
        %
        %                         subplot(1,3,2)
        %                         E = eig(full(A));
        %                         plot(real(E), imag(E), '*'), grid on
        %                         title('Eigenvalues of A')
        %
        %
        %                         subplot(1,3,3)
        %                         E = eig(full(M + tau * tau / 4 * A));
        %                         plot(real(E), imag(E), '*'), grid on
        %                         title('Eigenwerte von M + tau*tau / 4 A')
        %                         drawnow
        %
        %                         return
        %
        
        
        
        
        
        %% Calculation of the right hand side
        fprintf('Calculating rhs...')
        ntest = length(T_space_test) - splineOrder;
        Fvals = zeros((ntest-sum(test_offset)).^3,K);
        
        % Todo: Adapt this to the new form of the right hand side
        %         for k=1:K
        %             for i=1:length(f_time)
        %                 % First space dimension
        %                 tempX = zeros(ntest-sum(test_offset),1);
        %                 for j=1+test_offset(1):ntest-test_offset(2)
        %                     g = @(x) f_space_x{i}(x) * Ndiff(T_space_test, splineOrder,0,j,x);
        %                     for step = 0:splineOrder-test_offset(1)%bSplineOrder_space_test-1
        %                         tempX(j-test_offset(1)) = tempX(j-test_offset(1)) ...
        %                             + Gaussq(T_space_test(j+step),T_space_test(j+step+1),g,3);
        %                     end
        %                 end
        %                 % Second space dimension
        %                 tempY = zeros(ntest-sum(test_offset),1);
        %                 for j=1+test_offset(1):ntest-test_offset(2)
        %                     g = @(y) f_space_y{i}(y) * Ndiff(T_space_test, splineOrder,0,j,y);
        %                     for step = 0:splineOrder-test_offset(1)%bSplineOrder_space_test-1
        %                         tempY(j-test_offset(1)) = tempY(j-test_offset(1)) ...
        %                             + Gaussq(T_space_test(j+step),T_space_test(j+step+1),g,3);
        %                     end
        %                 end
        %                 % Third space dimension
        %                 tempZ = zeros(ntest-sum(test_offset),1);
        %                 for j=1+test_offset(1):ntest-test_offset(2)
        %                     g = @(z) f_space_z{i}(z) * Ndiff(T_space_test, splineOrder,0,j,z);
        %                     for step = 0:splineOrder-test_offset(1)%bSplineOrder_space_test-1
        %                         tempZ(j-test_offset(1)) = tempZ(j-test_offset(1)) ...
        %                             + Gaussq(T_space_test(j+step),T_space_test(j+step+1),g,3);
        %                     end
        %                 end
        %                 %                 Fvals(:, k) = Fvals(:,k) + kron(kron(tempX, tempY), tempZ) .* f_time{i}((k-1)/(K-1));
        %                 Fvals(:, k) = Fvals(:,k) + kron(kron(tempZ, tempY), tempX) .* f_time{i}((k-1)/(K-1));
        %             end
        %         end
        %         fprintf(' Done.\n')
        
        %% Calculation of u0 and u1
        fprintf('Calculating u0 and u1...')
        U0 = zeros(ntest-sum(test_offset),ntest-sum(test_offset),ntest-sum(test_offset));
        U1 = zeros(ntest-sum(test_offset),ntest-sum(test_offset),ntest-sum(test_offset));
        if inital_conditions
            for j=1+test_offset(1):ntest-test_offset(2)

                for k=1+test_offset(1):ntest-test_offset(2)
                    for l=1+test_offset(1):ntest-test_offset(2)
                        g = @(z)@(y)@(x) u0(x,y,z) .* ...
                            Ndiff(T_space_test, splineOrder,0,j,x) .* ...
                            Ndiff(T_space_test, splineOrder,0,k,y) .*...
                            Ndiff(T_space_test, splineOrder,0,l,z); %*u0
                        h = @(z)@(y)@(x) u1(x,y,z) .* ...
                            Ndiff(T_space_test, splineOrder,0,j,x) .* ...
                            Ndiff(T_space_test, splineOrder,0,k,y) .* ...
                            Ndiff(T_space_test, splineOrder,0,l,z); %*u1
                        
                        for stepX = 0:splineOrder-1
                            for stepY = 0:splineOrder-1
                                for stepZ = 0:splineOrder-1
                                    % Integrate from a to b
                                    a = [T_space_test(j+stepX); T_space_test(k+stepY); T_space_test(l+stepZ)];
                                    b = [T_space_test(j+stepX+1); T_space_test(k+stepY+1);T_space_test(l+stepZ+1) ];
                                    center = 0.5 * (a+b);
                                    if abs(u0(center(1),center(2),center(3))) > eps || ...
                                            abs(u0(a(1),a(2),a(3))) > eps || ...
                                            abs(u0(a(1),a(2),b(3))) > eps || ...
                                            abs(u0(a(1),b(2),a(3))) > eps || ...
                                            abs(u0(a(1),b(2),b(3))) > eps || ...
                                            abs(u0(b(1),b(2),b(3))) > eps || ...
                                            abs(u0(b(1),b(2),a(3))) > eps || ...
                                            abs(u0(b(1),a(2),b(3))) > eps || ...
                                            abs(u0(b(1),a(2),a(3))) > eps
                                        U0(j-test_offset(1), ...
                                        k-test_offset(1), ...
                                        l-test_offset(1)) = ...
                                        U0(j-test_offset(1), ...
                                        k-test_offset(1), ...
                                        l-test_offset(1)) + Gaussq(a,b,g,5);
                                    end
                                    
                                    if abs(u1(center(1),center(2),center(3))) > eps || ...
                                            abs(u1(a(1),a(2),a(3))) > eps || ...
                                            abs(u1(a(1),a(2),b(3))) > eps || ...
                                            abs(u1(a(1),b(2),a(3))) > eps || ...
                                            abs(u1(a(1),b(2),b(3))) > eps || ...
                                            abs(u1(b(1),b(2),b(3))) > eps || ...
                                            abs(u1(b(1),b(2),a(3))) > eps || ...
                                            abs(u1(b(1),a(2),b(3))) > eps || ...
                                            abs(u1(b(1),a(2),a(3))) > eps
                                        U1(j-test_offset(1), ...
                                        k-test_offset(1), ...
                                        l-test_offset(1)) = ...
                                        U1(j-test_offset(1), ...
                                        k-test_offset(1), ...
                                        l-test_offset(1)) + Gaussq(a,b,h,5);
                                    end
                                    
                                    
                                end
                            end
                        end
                    end
                end
            end
        end
        U0 = U0(:);
        U1 = U1(:);
        
        
        %         if inital_conditions
        %             tempXu0 = zeros(ntest-sum(test_offset),1);
        %             tempXu1 = zeros(ntest-sum(test_offset),1);
        %             for j=1+test_offset(1):ntest-test_offset(2)
        %                 g = @(x) u0_x(x) * Ndiff(T_space_test, splineOrder,0,j,x);
        %                 h = @(x) u1_x(x) * Ndiff(T_space_test, splineOrder,0,j,x);
        %
        %                 for step = 0:splineOrder-1
        %                     tempXu0(j-test_offset(1)) = tempXu0(j-test_offset(1)) ...
        %                         + Gaussq(T_space_test(j+step),T_space_test(j+step+1),g,3);
        %                     tempXu1(j-test_offset(1)) = tempXu1(j-test_offset(1)) ...
        %                         + Gaussq(T_space_test(j+step),T_space_test(j+step+1),h,3);
        %                 end
        %             end
        %
        %             tempYu0 = zeros(ntest-sum(test_offset),1);
        %             tempYu1 = zeros(ntest-sum(test_offset),1);
        %             for j=1+test_offset(1):ntest-test_offset(2)
        %                 g = @(y) u0_y(y) * Ndiff(T_space_test, splineOrder,0,j,y);
        %                 h = @(y) u1_y(y) * Ndiff(T_space_test, splineOrder,0,j,y);
        %
        %                 for step = 0:splineOrder-1
        %                     tempYu0(j-test_offset(1)) = tempYu0(j-test_offset(1)) ...
        %                         + Gaussq(T_space_test(j+step),T_space_test(j+step+1),g,3);
        %                     tempYu1(j-test_offset(1)) = tempYu1(j-test_offset(1)) ...
        %                         + Gaussq(T_space_test(j+step),T_space_test(j+step+1),h,3);
        %                 end
        %             end
        %
        %
        %             tempZu0 = zeros(ntest-sum(test_offset),1);
        %             tempZu1 = zeros(ntest-sum(test_offset),1);
        %             for j=1+test_offset(1):ntest-test_offset(2)
        %                 g = @(y) u0_y(y) * Ndiff(T_space_test, splineOrder,0,j,y);
        %                 h = @(y) u1_y(y) * Ndiff(T_space_test, splineOrder,0,j,y);
        %
        %                 for step = 0:splineOrder-1
        %                     tempZu0(j-test_offset(1)) = tempZu0(j-test_offset(1)) ...
        %                         + Gaussq(T_space_test(j+step),T_space_test(j+step+1),g,3);
        %                     tempZu1(j-test_offset(1)) = tempZu1(j-test_offset(1)) ...
        %                         + Gaussq(T_space_test(j+step),T_space_test(j+step+1),h,3);
        %                 end
        %             end
        %             U0 = kron(kron(tempXu0, tempYu0), tempZu0);
        %             U1 = kron(kron(tempXu1, tempYu1), tempZu1);
        %         else
        %             U0 = zeros((ntest-sum(test_offset)).^3,1);
        %             U1 = zeros((ntest-sum(test_offset)).^3,1);
        %         end
        fprintf(' Done.\n')
        
        
        %% Time stepping
        fprintf('Time stepping...')
        tic
        U = zeros(size(U0,1), K);
        % Save U0
        U(:,1) = pcg(M, U0);
        
        % Calculate the first step:
        [temp, flag] = pcg(M, U1, 1e-6, 1e5);
        
        if flag ~= 0
            warning('CG did not converge')
        end
        
        [temp2, flag] = pcg(M, Fvals(:,1)-funA(U(:,1)), 1e-6, 1e5);
        if flag ~= 0
            warning('CG did not converge')
        end
        
        U(:,2) = U(:,1) + tau*temp  + (tau*tau/2) * temp2;
        ts4 = tau*tau/4;
        funMA = @(x) funM(x) + ts4 * funA(x);
        
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
            
            [U(:, k), flag] = pcg(M + ts4 * A,  ts4* (F + 2*Fp+ Fpp) ...
                - ts4 * A * (2*U(:, k-1) + U(:, k-2)) ...
                + M * (2*U(:,k-1) - U(:, k-2)),...
                1e-6, ... % Tolerance
                1e5, ... % maxIt
                [], [], ... %Preconditioner
                U(:,k-1)); % Inital value
            
            if flag ~= 0
                warning('GMRES did not converge')
            end
            if mod(k, K/4) == 0
                fprintf('Time stepping progress: %.2f%%\n', 100*k/K)
            end
        end
        fprintf(' Done.\n')
        times(refinement_space, refinement_time) = toc;
        
        %% Get the solution
        fprintf('Get the solution...')
        
        % Add zeros coefficients for the boundaries
        space_size = nthroot(size(U,1),3) + sum(test_offset);
        u_full = zeros(space_size,space_size,space_size,size(U,2));
        u_full(1+test_offset(1):end-test_offset(2), ...
            1+test_offset(1):end-test_offset(2), ...
            1+test_offset(1):end-test_offset(2),:) = ...
            reshape(U, [nthroot(size(U,1),3),nthroot(size(U,1),3), ...
            nthroot(size(U,1),3), size(U,2)]);
        
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
        
        
        
        
        tic
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
                
                for k=1:size(u_full,3) % every space node in z
                    spline_z = splines(k,:);
                    
                    if laplace
                        spline_z_2 = splines_2(k,:);
                        
                        space_z_indices = min(find(spline_z, 1, 'first'), find(spline_z_2, 1, 'first')):...
                            max(find(spline_z, 1, 'last'), find(spline_z_2, 1, 'last'));
                    else
                        space_z_indices = find(spline_z, 1, 'first'):find(spline_z, 1, 'last');
                    end
                    
                    for l=1:K % Every time step
                        if u_full(i,j,k,l) < eps && u_full(i,j,k,l) > - eps
                            continue;
                        end
                        
                        if laplace
                            spline = kron(spline_z_2(space_z_indices), kron(spline_y(space_y_indices), spline_x(space_x_indices))) ...
                                + kron(spline_z(space_z_indices), kron(spline_y_2(space_y_indices), spline_x(space_x_indices))) ...
                                + kron(spline_z(space_z_indices), kron(spline_y(space_y_indices), spline_x_2(space_x_indices)));
                        else
                            
                            spline = kron(spline_z(space_z_indices), kron(spline_y(space_y_indices), spline_x(space_x_indices)));
                            
                        end
                        
                        spline = reshape(spline, ...
                            [length(space_x_indices), length(space_y_indices),...
                            length(space_z_indices)]);
                        sol(space_x_indices, space_y_indices, space_z_indices, l) = ...
                            sol(space_x_indices, space_y_indices, space_z_indices, l) + ...
                            u_full(i,j,k,l) * spline;
                    end
                end
            end
        end
        
        toc
        
        fprintf(' Done.\n')
        
        
        % Linearaly interpolate the solution to match Kref time steps
        [X, Y, Z, T] = ndgrid(x, y, z, linspace(0,1, K));
        [Xq, Yq, Zq, Tq] = ndgrid(x, y, z, linspace(0,1, Kref));
        solInt = interpn(X, Y, Z, T, sol,...
            Xq, Yq, Zq, Tq);
        %     solInt = sol(:,:,:,end);
        
        l2error(refinement_space, refinement_time) = sqrt(mean((solInt-sol_ref).^2, 'all'));
        
        diag(l2error)
    end
end
% figure
% semilogy(space_refinements, l2error(space_refinements,:), 'o--', 'LineWidth', 2)
% xlabel('Space refinement', 'Interpreter', 'Latex', 'FontSize', 20)
% ylabel('$L_2$ error', 'Interpreter', 'Latex', 'FontSize', 20)
% title('3D wave equation - Laplace case',  'Interpreter', 'Latex', 'FontSize', 20)
% grid on
% legend('Time refinement = 1','Time refinement = 2','Time refinement = 3',...
%     'Time refinement = 4','Time refinement = 5','Time refinement = 6',...
%     'Time refinement = 7')
%
% figure
% semilogy(space_refinements, times_gmres(space_refinements,:), 'o--', 'LineWidth', 2)
% xlabel('Space refinement', 'Interpreter', 'Latex', 'FontSize', 20)
% ylabel('Solving time $[s]$', 'Interpreter', 'Latex', 'FontSize', 20)
% title('3D wave equation - Laplace case',  'Interpreter', 'Latex', 'FontSize', 20)
% grid on
% legend('Time refinement = 1','Time refinement = 2','Time refinement = 3',...
%     'Time refinement = 4','Time refinement = 5','Time refinement = 6',...
%     'Time refinement = 7')
%
% figure
% semilogy(time_refinements, times_gmres(:,time_refinements), 'o--', 'LineWidth', 2)
% xlabel('Time refinement', 'Interpreter', 'Latex', 'FontSize', 20)
% ylabel('Solving time $[s]$', 'Interpreter', 'Latex', 'FontSize', 20)
% title('3D wave equation - Laplace case',  'Interpreter', 'Latex', 'FontSize', 20)
% grid on
% legend('Space refinement = 1' , 'Space refinement = 2','Space refinement = 3',...
%     'Space refinement = 4','Space refinement = 5')
%
% figure
% for i=1:100
% s = surf(squeeze(solInt(25,:,:,i)));
% s.EdgeAlpha = 0;
% title('Numerical solution')
% drawnow
% pause(0.5)
% end


% figure
% loglog(reshape(l2error, [], 1), reshape(times_gmres, [], 1), '*'), grid on
% xlabel('L2 Error')
% ylabel('CPU Time')

%% Write the whole data into files

% refinement = space_refinements';
% error = diag(l2error);
% time = diag(times);
% writetable(table(refinement,time, error), ...
%     '../data/3Dexample1-timestepping','Delimiter',' ')

semilogy(diag(l2error), '*--', 'LineWidth', 2), grid on

% Produces NaN for r = 0
function sol = dAlemenbertSphere(r,t,c, u0)
% if r == 0
%     sol = NaN;
% elseif r > c * t
    sol = 1./(2.*r) .*( (r+c.*t) .* u0(r + c.*t) + (r - c.*t).*u0(r-c.*t));
% else
%     sol = 1./(2.*r) .*( (r+c.*t) .* u0(r + c.*t) - (c.*t - r).*u0(c.*t-r));
% end
end

