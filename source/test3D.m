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

problemConfiguration.d = 3; % Dimension

% Example 1 (very smooth)
mm = 3; nn = 2; oo = 1;
problemConfiguration.f_time = {@(t) 2+ 0.*t; ...
    @(t) (t.^2)*(mm^2 + nn^2 + oo^2)*(pi^2)};
problemConfiguration.f_space_x = {@(x) sin(mm*pi*x); @(x) sin(mm*pi*x)};
problemConfiguration.f_space_y = {@(y) sin(nn*pi*y); @(y) sin(nn*pi*y)};
problemConfiguration.f_space_z = {@(z) sin(oo*pi*z); @(z) sin(oo*pi*z)};
problemConfiguration.u_0_x = @(x) 0*x;
problemConfiguration.u_0_y = @(y) 0*y;
problemConfiguration.u_0_z = @(z) 0*z;
problemConfiguration.u_1_x = @(x) 0*x;
problemConfiguration.u_1_y = @(y) 0*y;
problemConfiguration.u_1_z = @(z) 0*z;

problemConfiguration.u_analytical = @(x,y,z,t) t.^2 .* ...
    sin(mm*pi*x).*sin(nn*pi*y).*sin(oo*pi*z);
problemConfiguration.has_analytical_solution = true;

% Example 2
% f_time = {@(t) 0*t};
% f_space_x = {@(x) 0*x};
% f_space_y = {@(y) 0*y};
% f_space_z = {@(z) 0*z};
% u_0_x = @(x) x* (x < 0.5) + (1-x).*(x > 0.5);
% u_0_y = @(y) y* (y < 0.5) + (1-y).*(y > 0.5);
% u_0_z = @(z) z* (z < 0.5) + (1-z).*(z > 0.5);
% u_1_x = @(x) 0*x;
% u_1_y = @(y) 0*y;
% u_1_z = @(z) 0*z;
% inital_conditions = true;
% has_analytical_solution = false;
% load('3D-example2.mat')

tolerance = 1e-5; % For the CG or Galerkin method
maxIt = 50;

resolution.x = 5;
resolution.y = 5;
resolution.z = 5;
resolution.t = 5;

plotting = false;
% Specify the resolution for plotting and the error calculation

% x = linspace(0,1,2^6+1); % Increase for large space refinement levels
% y = x;
% z = x;
% t = linspace(0,1,2^6+1);

l2error = true; % Calculate the l2 error'

exactFlag = false;


for refinementLevel_space = 1:6
    for refinementLevel_time = refinementLevel_space
        
        problemConfiguration.refinementLevel_space = refinementLevel_space;
        problemConfiguration.refinementLevel_time = refinementLevel_time;
        
        fprintf('Creating problem...')
        problem = create3DWaveProblem(problemConfiguration);
        fprintf(' Done!\n')
        
        
        fprintf('\n')
        fprintf('\n')
        disp(['Size of the time matrices: ' ...
            num2str(length(problem.Q_time)) ' x ' ...
            num2str(length(problem.Q_time))])
        disp(['Size of the space matrices: ' ...
            num2str(length(problem.Q_space)) ' x ' ...
            num2str(length(problem.Q_space))])
        
        %% Sove the linear equation system
        
        
        % B = kron(Q_time, M_space) ...
        %     + kron(N_time,N_space')...
        %     + kron(N_time',N_space)...
        %     + kron(M_time, Q_space);
        % tic;
        % U = B \ rhs(:);
        % time = toc;
        % disp(['Solving time: ' num2str(time)])
        
        
        % Specify a function handle for the pcg methods
        funA=@(X)( problem.M_space * X * problem.Q_time' ...
            + problem.A_space' * X * problem.D_time' ...
            + problem.A_space * X * problem.D_time ...
            + problem.Q_space * X * problem.M_time');
        %  funA=@(X)( M_space*X*Q_time'+N_space'*X*N_time'+N_space*X*N_time+Q_space*X*M_time');
        
        
        funB=@(x) reshape(funA(reshape(x, [nthroot(length(x),4)^3, nthroot(length(x),4)])), ...
            [length(x) 1]); % only works if space ref = time ref
        
        rhsfull=reshape(problem.rhs, [numel(problem.rhs) 1]);
        rhsfull=rhsfull/norm(rhsfull);
        %         params.Q_time=Q_time; params.M_time=M_time;
        %         params.M_space=M_space; params.Q_space=Q_space;
        %         params.A_space=N_space; params.D_time=N_time;
        %         params.M_spacelocal=M_space_local;
        %         params.Q_spacelocal=Q_space_local;
        %         params.N_spacelocal=N_space_local;
        %         params.d=3;
        %         exactflag=0;
        
        for ii=1:3 % Test all three solutions
            switch ii
                case 1
                    % Product operator preconditioner
                    fprintf('Starting CG optimal')
                    problem.precond='optimal';
                    tt=tic;
                    [U, iter_opt(refinementLevel_space, ...
                        refinementLevel_time)]=pcg_fun4(funA,rhsfull, ...
                        0*rhsfull,problem,30,tolerance,...
                        size(problem.M_space,1), ...
                        size(problem.M_time,2),exactFlag);
                    U_cg_optimal=U(:)*norm(problem.rhs(:));
                    solvingTimes_cg_opt(refinementLevel_space, refinementLevel_time) = toc(tt);
                    
                    solvingErrorCGopt(refinementLevel_space, refinementLevel_time) = ...
                        norm(funB(U_cg_optimal) - reshape(problem.rhs, [numel(problem.rhs), 1])) / norm(reshape(problem.rhs, [numel(problem.rhs), 1]));
                    
                    fprintf('Time refinement: %d, Space refinement: %d, Time to solve: %f\n', ...
                        refinementLevel_time, refinementLevel_space,  solvingTimes_cg_opt(refinementLevel_space, refinementLevel_time))
                    
                case 2
                    % Lyap operator preconditioner
                    fprintf('Starting CG lyaponov')
                    problem.precond='lyap';
                    tt=tic;
                    [U, iter_lyap(refinementLevel_space, ...
                        refinementLevel_time)]=pcg_fun4(funA,rhsfull,...
                        0*rhsfull,problem,30,tolerance, ...
                        size(problem.M_space,1), ...
                        size(problem.M_time,2),exactFlag);
                    solvingTimes_cg_lyap(refinementLevel_space, refinementLevel_time) = toc(tt);
                    U_cg_lyap=U(:)*norm(problem.rhs(:));
                    
                    solvingErrorCGlyap(refinementLevel_space, refinementLevel_time) = ...
                        norm(funB(U_cg_lyap) - reshape(problem.rhs, [numel(problem.rhs), 1])) / norm(reshape(problem.rhs, [numel(problem.rhs), 1]));
                    
                    fprintf('Time refinement: %d, Space refinement: %d, Time to solve: %f\n', ...
                        refinementLevel_time, refinementLevel_space,  solvingTimes_cg_lyap(refinementLevel_space, refinementLevel_time))
                case 3
                    % Galerkin
                    fprintf('Starting Galerkin')
                    [uu,ss,vv]=svds(problem.rhs,1);
                    rhs1=uu(:,1)*sqrt(ss(1,1));
                    rhs2=vv(:,1)*sqrt(ss(1,1));
                    path(path,'./valeria/')
                    maxitG=50;
                    tolG=tolerance;
                    info=1;
                    tt=tic;
                    [X1,X2,restot,iter_galerkin(refinementLevel_space, ...
                        refinementLevel_time)]=Galerkin3(problem.M_space,...
                        2*problem.A_space,problem.Q_space,problem.Q_time, ...
                        (problem.D_time+problem.D_time')/2,...
                        problem.M_time,rhs1,rhs2,maxitG,tolG,info);
                    solvingTimes_Galerkin(refinementLevel_space, refinementLevel_time) = toc(tt);
                    U=X1*X2';
                    U_galerkin=U(:);
                    
                    solvingErrorGalerkin(refinementLevel_space, refinementLevel_time) = ...
                        norm(funB(U_galerkin) - reshape(problem.rhs, [numel(problem.rhs), 1])) / norm(reshape(problem.rhs, [numel(problem.rhs), 1]));
                    
                    
                    fprintf('Galerkin: Time refinement: %d, Space refinement: %d, Time to solve: %f\n', ...
                        refinementLevel_time, refinementLevel_space,  solvingTimes_Galerkin(refinementLevel_space, refinementLevel_time))
                    
            end
            
            
            %% Plot the solution (only recommended for small problems)
            if ~plotting && ~l2error
                continue;
            else
                
                
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
                
                sol = get3Dsolution(problem, U, resolution);
            
                if l2error
                    
                    switch ii
                        case 1
                            errorCGopt(refinementLevel_space, refinementLevel_time) = calculate3DL2Error(problem, sol);
                        case 2
                            errorCGlyap(refinementLevel_space, refinementLevel_time) = calculate3DL2Error(problem, sol);
                        case 3
                            errorGalerkin(refinementLevel_space, refinementLevel_time) = calculate3DL2Error(problem, sol);
                    end
                    
                end
                
            end
        end
    end
end


errorCGopt = diag(errorCGopt);
errorCGlyap = diag(errorCGlyap);
errorGalerkin = diag(errorGalerkin);

solvingErrorCGopt = diag(solvingErrorCGopt);
solvingErrorCGlyap = diag(solvingErrorCGlyap);
solvingErrorGalerkin = diag(solvingErrorGalerkin);


solvingTimes_cg_opt = diag(solvingTimes_cg_opt);
solvingTimes_cg_lyap = diag(solvingTimes_cg_lyap);
solvingTimes_Galerkin = diag(solvingTimes_Galerkin);

figure
subplot(1,2,1)
semilogy(errorCGopt, '*--', 'LineWidth', 3), hold on, grid on
semilogy(errorCGlyap, 'd--', 'LineWidth', 3)
semilogy(errorGalerkin, 's--', 'LineWidth', 3)
title('L2 Error')
legend('CG opt', 'CG lyap', 'Galerkin')
xlabel('Space and time refinements')

subplot(1,2,2)
semilogy(solvingErrorCGopt, '*--', 'LineWidth', 3), hold on, grid on
semilogy(solvingErrorCGlyap, 'd--', 'LineWidth', 3)
semilogy(solvingErrorGalerkin, 's--', 'LineWidth', 3)
title('||Ax-b||/||b||')
legend('CG opt', 'CG lyap', 'Galerkin')
xlabel('Space and time refinements')


figure
semilogy(solvingTimes_cg_opt, '*--', 'LineWidth', 3), hold on, grid on
semilogy(solvingTimes_cg_lyap, 'd--', 'LineWidth', 3)
semilogy(solvingTimes_Galerkin, 's--', 'LineWidth', 3)
title('Walltime')
legend('CG opt', 'CG lyap', 'Galerkin')
xlabel('Space and time refinements')


% figure
% subplot(1,2,1)
% surf(l2err_galekin)
% xlabel('space refinement')
% ylabel('time refinement')
% title('L2 error')
% 
% subplot(1,2,2)
% surf(solvingTimes_Galerkin)
% xlabel('space refinement')
% ylabel('time refinement')
% title('Galerkin CPU Time')
%
% figure
% subplot(1,2,1)
% semilogy(l2err_galekin', '*--', 'LineWidth', 2)
% grid on
% xlabel('Space refinements')
% legend('Time refinement = 1' , 'Time refinement = 2', 'Time refinement = 3',...
%     'Time refinement = 4','Time refinement = 5','Time refinement = 6', ...
%     'Time refinement = 7','Time refinement = 8','Time refinement = 9', ...
%     'Time refinement = 10')
% title('L2 error')
%
% subplot(1,2,2)
% semilogy(solvingTimes_Galerkin', '*--', 'LineWidth', 2)
% grid on
% xlabel('Space refinements')
% legend('Time refinement = 1' , 'Time refinement = 2', 'Time refinement = 3',...
%     'Time refinement = 4','Time refinement = 5','Time refinement = 6', ...
%     'Time refinement = 7','Time refinement = 8','Time refinement = 9', ...
%     'Time refinement = 10')
% title('Time to solve')


%
% loglog(reshape(l2errCG_opt, [], 1)', reshape(solvingTimes_cg_optimal,[], 1)', '*'), hold on
% %loglog(reshape(l2errCG_lyap, [], 1), reshape(solvingTimes_cg_lyap,[], 1), 's:')
% loglog(reshape(l2err_galekin, [], 1), reshape(solvingTimes_Galerkin,[], 1), 'd')
% legend('CG opt', 'Galerkin')
% xlabel('L2 Error')
% ylabel('CPU time')
% grid on

%% Write the whole data into files

timeCGopt = diag(solvingTimes_cg_opt);
% timeCGlyap = diag(solvingTimes_cg_lyap);
% timeGalerkin = diag(solvingTimes_Galerkin);

iterCGopt = diag(iter_opt);
% iterCGlyap = diag(iter_lyap);
% iterGalerkin = diag(iter_galerkin);


refinement = (1:4)';
% writetable(table(refinement,iterCGopt, timeCGopt, errorCGopt), ...
%     'data/3Dexample2-opt-exact-0','Delimiter',' ')
% writetable(table(refinement,iterCGlyap, timeCGlyap, errorCGlyap), ...
%     'data/3Dexample1-lyap-exact-0','Delimiter',' ')
% writetable(table(refinement,iterGalerkin, timeGalerkin, errorGalerkin), ...
%     'data/3Dexample1-galerkin','Delimiter',' ')