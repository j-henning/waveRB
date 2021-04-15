clear
close all
clc

addpath('../Toolkit-Splines');
addpath('../Toolkit-Splines/Utilities');
addpath('../Solvers');





%% Specify all needed input data

% Example 1
f_time = {@(t) 2+0*t, @(t) t.^2};
f_space = {@(x) sin(2*pi*x), @(x) 4*pi^2*sin(2*pi*x)};
u_0_x = @(x) 0.*x;
u_1_x = @(x) 0.*x;
u_analytical = @(x,t) t.^2 .* sin(2*pi*x);
dAlemebert = false;
has_analytical_solution = true;


% Example 2
% f_time = {@(t) 0*t};
% f_space = {@(x) 0*x};
% u_0_x = @(x) x* (x < 0.5) + (1-x).*(x > 0.5);
% u_1_x = @(x) 0*x;
% inital_conditions = true;
% has_analytical_solution = false;
% load('1D-example2')

% Example 3
% f_time = {@(t) 0*t};
% f_space = {@(x) 0*x};
% u_0_x = @(x) (x > 0.25 &&  x < 0.75);
% u_1_x = @(x) 0*x;
% inital_conditions = true;
% has_analytical_solution = false;
% load('1D-example3')

bSplineOrder_time = 3;
bSplineOrder_space = 3;

offset_time_ansatz = [0, 2];
offset_time_test = [0, 2];
offset_space_ansatz = [1, 1];
offset_space_test = [1, 1];

plotting = false;

% Specify the solution for plotting and the error calculation

ref_plotting = 11;
x = linspace(0,1,2^ref_plotting+1);
t = linspace(0,1,2^ref_plotting+1);


l2error = true; % Calculate the l2 error

tolerance = 1e-10;
maxIt = 200;
exactFlag = true;


space_refinement = 1:10;
for refinementLevel_space = space_refinement
    for refinementLevel_time = refinementLevel_space
        
        
        
        
        
        %% Calculate the matrices
        disp('Compute matrices')
        
        [Q_time, M_time, N_time, Q_space, M_space, N_space, rhs] ...
            = computeMatrices1D(refinementLevel_time, bSplineOrder_time, ...
            refinementLevel_space, bSplineOrder_space, f_time, f_space, ...
            u_0_x, u_1_x);
        
        disp('done')
        for ii=1:3 % Test all three solutions
            
            %% Sove the linear equation system
            
            % tic;
%                         B = kron(Q_time, M_space) ...
%                             + kron(N_time,N_space')...
%                             + kron(N_time',N_space)...
%                             + kron(M_time, Q_space);
%                         return
            %             B=sparse(B);
            %             tic
            %             U = B \ reshape(rhs, [numel(rhs) 1]);
            %             toc
            %             norm(B*U-reshape(rhs, [numel(rhs) 1]))
            % time = toc;
            % disp(['Solving time with mldivide: ' num2str(time)])
            
            funA=@(X)( M_space*X*Q_time'+N_space'*X*N_time'+N_space*X*N_time+Q_space*X*M_time');
            funB=@(x) reshape(funA(reshape(x, [sqrt(length(x)), sqrt(length(x))])), ...
                [length(x) 1]); % only works if space ref = time ref
            
            
            
            rhsfull=reshape(rhs, [numel(rhs) 1]);
            rhsfull=rhsfull/norm(rhsfull);
            params.Q_time=Q_time; params.M_time=M_time;
            params.M_space=M_space; params.Q_space=Q_space;
            params.A_space=N_space; params.D_time=N_time;
            params.d=1;
       
            switch ii
                case 1
                    % Product operator preconditioner
                    
                    params.precond='optimal';
                    tt=tic;
                    [U, iterCGopt(refinementLevel_space, refinementLevel_time)]= ...
                        pcg_fun4(funA,rhsfull,0*rhsfull,params,maxIt,tolerance,...
                        size(M_space,1),size(M_time,2),exactFlag);
                    U_cg_optimal=U(:)*norm(rhs(:));
                    timeCGopt(refinementLevel_space, refinementLevel_time) = toc(tt);
                    solvingErrorCGopt(refinementLevel_space, refinementLevel_time) = ...
                        norm(funB(U_cg_optimal) - reshape(rhs, [numel(rhs), 1])) ...
                        / norm(reshape(rhs, [numel(rhs), 1]));
                    
                    fprintf('Time refinement: %d, Space refinement: %d, Time to solve: %f\n', ...
                        refinementLevel_time, refinementLevel_space,  timeCGopt(refinementLevel_space, refinementLevel_time))
                    
                case 2
                    % Lyap operator preconditioner
                    params.precond='lyap';
                    tt=tic;
                    [U, iterCGlyap(refinementLevel_space, refinementLevel_time)]=...
                        pcg_fun4(funA,rhsfull,0*rhsfull,params,maxIt,tolerance,size(M_space,1),size(M_time,2),exactFlag);
                    timeCGlyap(refinementLevel_space, refinementLevel_time) = toc(tt);
                    U_cg_lyap=U(:)*norm(rhs(:));
                    solvingErrorCGlyap(refinementLevel_space, refinementLevel_time) = ...
                        norm(funB(U_cg_lyap) - reshape(rhs, [numel(rhs), 1])) ...
                        / norm(reshape(rhs, [numel(rhs), 1]));
                    
                    fprintf('Time refinement: %d, Space refinement: %d, Time to solve: %f\n', ...
                        refinementLevel_time, refinementLevel_space,  timeCGlyap(refinementLevel_space, refinementLevel_time))
                    
                case 3
                    % Galerkin
                    [uu,ss,vv]=svds(rhs,1);
                    rhs1=uu(:,1)*sqrt(ss(1,1));
                    rhs2=vv(:,1)*sqrt(ss(1,1));
                    path(path,'./valeria/')
                    tolG=tolerance;
                    info=1;
                    tt=tic;
                    [X1,X2,restot,iterGalerkin(refinementLevel_space, refinementLevel_time)]= ...
                        Galerkin3(M_space,2*N_space,Q_space,Q_time,(N_time+N_time')/2,M_time,rhs1,rhs2,maxIt,tolG,info);
                    timeGalerkin(refinementLevel_space, refinementLevel_time) = toc(tt);
                    U=X1*X2';
                    U_galerkin=U(:);
                    solvingErrorGalerkin(refinementLevel_space, refinementLevel_time) = ...
                        norm(funB(U_galerkin) - reshape(rhs, [numel(rhs), 1]))...
                        / norm(reshape(rhs, [numel(rhs), 1]));
                    
                    
                    
                    
                    
                    fprintf('Galerkin: Time refinement: %d, Space refinement: %d, Time to solve: %f\n', ...
                        refinementLevel_time, refinementLevel_space,  timeGalerkin(refinementLevel_space, refinementLevel_time))
            end
            
            
            switch ii
                case 1
                    U = U_cg_optimal;
                    fprintf('Testing the cg_opt solution\n');
                case 2
                    U = U_cg_lyap;
                    fprintf('Testing the cg_lyap solution\n');
                case 3
                    U = U_galerkin;
                    fprintf('Testing the galerkin solution\n');
            end
            
            %% Plot the solution
            if ~plotting && ~l2error
                continue
            end
            
            tic;
            % Recalculate some needed variables
            T_time = InitUniNodes(refinementLevel_time, bSplineOrder_time);
            T_space = InitUniNodes(refinementLevel_space, bSplineOrder_space);
            
            
            nsol_time = length(T_time) - bSplineOrder_time;
            nsol_space = length(T_space) - bSplineOrder_space;
            
            
            dim_time =  nsol_time-(offset_time_ansatz(1)+offset_time_ansatz(2));
            
            u = reshape(full(U)', numel(U) / dim_time, dim_time);
            
            % Preparing boundary for plot
            u_full = zeros(nsol_space, nsol_time);
            
            
            u_full(1+offset_space_ansatz(1):nsol_space-offset_space_ansatz(2),...
                1+offset_time_ansatz(1):nsol_time-offset_time_ansatz(2)) = u;
            
            
            sol = zeros(length(x), length(t));
            
            total_iterations = numel(u_full);
            iteration = 0;
            
            
            % Precalculate the splines
            
            splines_space = zeros(size(u_full,1), length(x));
            splines_2_space = zeros(size(u_full,1), length(x));
            
            for i=1:size(splines_space,1)
                if i < 5 || i > size(splines_space,1) - 5
                    for j=1:size(splines_space,2)
                        g = @(s) Ndiff(T_space, bSplineOrder_space,0,i,(s==0)*(s+eps) + (s==1)*(s-eps) + (s > 0 && s < 1)*s);
                        h = @(s) Ndiff(T_space, bSplineOrder_space,2,i,(s==0)*(s+eps) + (s==1)*(s-eps) + (s > 0 && s < 1)*s);
                        
                        splines_space(i,j) = g(x(j));
                        splines_2_space(i,j) = h(x(j));
                    end
                else
                    offset = 2^(ref_plotting - refinementLevel_space);
                    splines_space(i, offset+1:end) = splines_space(i-1, 1:end-offset);
                    splines_2_space(i, offset+1:end) = splines_2_space(i-1, 1:end-offset);
                end
                
                
                
            end
            
            splines_time = zeros(size(u_full,2), length(t));
            splines_2_time = zeros(size(u_full,2), length(t));
            
            for i=1:size(splines_time,1)
                if i < 5 || i > size(splines_time,1) - 5
                    for j=1:size(splines_time,2)
                        g = @(s) Ndiff(T_time, bSplineOrder_time,0,i,(s==0)*(s+eps) + (s==1)*(s-eps) + (s > 0 && s < 1)*s);
                        h = @(s) Ndiff(T_time, bSplineOrder_time,2,i,(s==0)*(s+eps) + (s==1)*(s-eps) + (s > 0 && s < 1)*s);
                        
                        splines_time(i,j) = g(x(j));
                        splines_2_time(i,j) = h(x(j));
                        
                    end
                else
                    offset = 2^(ref_plotting - refinementLevel_time);
                    splines_time(i, offset+1:end) = splines_time(i-1, 1:end-offset);
                    splines_2_time(i, offset+1:end) = splines_2_time(i-1, 1:end-offset);
                end
            end
            
            
            for i=1:size(u_full,1) % every space node in x
                if mod(i, size(u_full,1)/8) == 0
                    disp(['Plotting preparation: ' num2str(i / size(u_full,1) * 100) '%'])
                end
                spline_x = splines_space(i,:);
                spline_x_2 = splines_2_space(i,:);
                space_x_indices = min(find(spline_x, 1, 'first'), find(spline_x_2, 1, 'first')):...
                    max(find(spline_x, 1, 'last'), find(spline_x_2, 1, 'last'));
                spline_x = spline_x(space_x_indices);
                spline_x_2 = spline_x_2(space_x_indices);
                
                %                 for j=1:size(u_full, 2) % every space node in y
                %
                %                     spline_y = splines_space(j,:);
                %                     spline_y_2 = splines_2_space(j,:);
                %                     space_y_indices = min(find(spline_y, 1, 'first'), find(spline_y_2, 1, 'first')):...
                %                         max(find(spline_y, 1, 'last'), find(spline_y_2, 1, 'last'));
                %                     spline_y = spline_y(space_y_indices);
                %                     spline_y_2 = spline_y_2(space_y_indices);
                
                for k=1:size(u_full,2) % every time node
                    if abs(u_full(i,k)) < eps
                        continue;
                    end
                    
                    spline_t = splines_time(k,:);
                    spline_t_2 = splines_2_time(k,:);
                    time_indices = min(find(spline_t, 1, 'first'), find(spline_t_2, 1, 'first')):...
                        max(find(spline_t, 1, 'last'), find(spline_t_2, 1, 'last'));
                    spline_t = spline_t(time_indices);
                    spline_t_2 = spline_t_2(time_indices);
                    
                    % Build the whole spline
                    spline = kron(spline_t_2, spline_x) ...
                        - kron(spline_t, spline_x_2);
                    
                    spline = reshape(spline, ...
                        [length(space_x_indices),  length(time_indices)]);
                    sol(space_x_indices, time_indices) = ...
                        sol(space_x_indices, time_indices) + ...
                        u_full(i,k) * spline;
                end
                %  end
            end
            
            
            
            if plotting
                %   ylimits = [-1 1];
                ylimits = [min(sol, [], 'all'), max(sol, [], 'all')];
                [X, Y] = meshgrid(y, x);
                fh = figure();
                fh.WindowState = 'maximized';
                for i=1:2:length(t)
                    if has_analytical_solution
                        subplot(1,2,1)
                    end
                    plot(sol(:,i))
                    xlabel('x')
                    ylabel('y')
                    title(['Numerical solution at t = ' num2str(t(i))])
                    ylim(ylimits)
                    
                    if has_analytical_solution
                        subplot(1,2,2)
                        xval = linspace(0,1,100);
                        plot(xval, u_analytical(xval, t(i)))
                        xlabel('x')
                        ylabel('y')
                        title(['Analytical solution at t = ' num2str(t(i))])
                        ylim(ylimits)
                    end
                    drawnow
                    pause(1)
                end
                
              
            end
            toc
            
            
            if l2error
                if has_analytical_solution
                    [X, T] = ndgrid(x, t);
                    sol_ana = u_analytical(X,T);
                else
                    sol_ana = sol_ref;
                end
              
                switch ii
                    case 1
                        errorCGopt(refinementLevel_space, refinementLevel_time) = sqrt(mean( (sol-sol_ana).^2, 'all'));
                    case 2
                        errorCGlyap(refinementLevel_space, refinementLevel_time) = sqrt(mean( (sol-sol_ana).^2, 'all'));
                    case 3
                        errorGalerkin(refinementLevel_space, refinementLevel_time) = sqrt(mean( (sol-sol_ana).^2, 'all'));
                end
            end
            
        end
    end
end



%% Write the whole data into files

refinement = space_refinement';

errorCGopt = diag(errorCGopt);
timeCGopt = diag(timeCGopt);
iterCGopt = diag(iterCGopt);

errorCGlyap = diag(errorCGlyap);
timeCGlyap = diag(timeCGlyap);
iterCGlyap = diag(iterCGlyap);

errorGalerkin = diag(errorGalerkin);
timeGalerkin = diag(timeGalerkin);
iterGalerkin = diag(iterGalerkin);

solvingErrorCGopt = diag(solvingErrorCGopt);
solvingErrorCGlyap = diag(solvingErrorCGlyap);
solvingErrorGalerkin = diag(solvingErrorGalerkin);

writetable(table(refinement,iterCGopt, timeCGopt, errorCGopt), ...
    'data/1Dexample1-opt-exact-1-1e-10','Delimiter',' ')
writetable(table(refinement,iterCGlyap, timeCGlyap, errorCGlyap), ...
    'data/1Dexample1-lyap-exact-1-1e-10','Delimiter',' ')
writetable(table(refinement,iterGalerkin, timeGalerkin, errorGalerkin), ...
    'data/1Dexample1-galerkin-1e-10','Delimiter',' ')
