function [] = test3D(example, method, maxRefinement, tolerance, maxIt, exactFlag)

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

switch example
    % Example 1 (very smooth)
    case 1
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
        
        
    case 2
        problemConfiguration.f_time = {@(t) 0*t};
        problemConfiguration.f_space_x = {@(x) 0*x};
        problemConfiguration.f_space_y = {@(y) 0*y};
        problemConfiguration.f_space_z = {@(z) 0*z};
        problemConfiguration.u_0_x = @(x) x* (x < 0.5) + (1-x).*(x > 0.5);
        problemConfiguration.u_0_y = @(y) y* (y < 0.5) + (1-y).*(y > 0.5);
        problemConfiguration.u_0_z = @(z) z* (z < 0.5) + (1-z).*(z > 0.5);
        problemConfiguration.u_1_x = @(x) 0*x;
        problemConfiguration.u_1_y = @(y) 0*y;
        problemConfiguration.u_1_z = @(z) 0*z;
        problemConfiguration.has_analytical_solution = false;
        problemConfiguration.referenceSolutionPath = '3D-example2';
        
    case 3
        problemConfiguration.f_time = {@(t) 0*t};
        problemConfiguration.f_space_x = {@(x) 0*x};
        problemConfiguration.f_space_y = {@(y) 0*y};
        problemConfiguration.f_space_z = {@(z) 0*z};
        problemConfiguration.u_0_x = @(x) 1 * (x > 0.25) .* (x < 0.75);
        problemConfiguration.u_0_y = @(y) 1 * (y > 0.25) .* (y < 0.75);
        problemConfiguration.u_0_z = @(z) 1 * (z > 0.25) .* (z < 0.75);
        problemConfiguration.u_1_x = @(x) 0*x;
        problemConfiguration.u_1_y = @(y) 0*y;
        problemConfiguration.u_1_z = @(z) 0*z;
        problemConfiguration.has_analytical_solution = false;
        problemConfiguration.referenceSolutionPath = '3D-example3';
      
end

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

resolution.x = 7;
resolution.y = 7;
resolution.z = 7;
resolution.t = 7;

% Create all the solution vectors (needed for writing the tables)
refinements = (1:maxRefinement)';
iterCGopt = zeros(maxRefinement,1);
timeCGopt = zeros(maxRefinement,1);
errorCGopt = zeros(maxRefinement,1);
solvingErrorCGopt = zeros(maxRefinement,1);

iterCGlyap = zeros(maxRefinement,1);
timeCGlyap = zeros(maxRefinement,1);
errorCGlyap = zeros(maxRefinement,1);
solvingErrorCGlyap = zeros(maxRefinement,1);

iterGalerkin = zeros(maxRefinement,1);
timeGalerkin = zeros(maxRefinement,1);
errorGalerkin = zeros(maxRefinement,1);
solvingErrorGalerkin = zeros(maxRefinement,1);


for refinement = 1:maxRefinement
    refinementLevel_space = refinement;
    refinementLevel_time = refinement;
    
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
    
    switch method
        case 1
            % Product operator preconditioner
            fprintf('Starting CG optimal\n')
            problem.precond='optimal';
            tt=tic;
            [U, iterCGopt(refinement)]=pcg_fun4(funA,rhsfull, ...
                0*rhsfull,problem,maxIt,tolerance,...
                size(problem.M_space,1), ...
                size(problem.M_time,2),exactFlag);
            U_cg_optimal=U(:)*norm(problem.rhs(:));
            timeCGopt(refinement) = toc(tt);
            
            solvingErrorCGopt(refinement) = ...
                norm(funB(U_cg_optimal) - reshape(problem.rhs, [numel(problem.rhs), 1])) / norm(reshape(problem.rhs, [numel(problem.rhs), 1]));
            
            fprintf('Refinement: %d, Time to solve: %f\n', ...
                refinement,  timeCGopt(refinement))
            
            fprintf('Testing the cg_opt solution\n');
            sol = get3Dsolution(problem, U_cg_optimal, resolution);
            
            errorCGopt(refinement) = calculate3DL2Error(problem, sol);
            
            % Write everything in a file in case the whole script does
            % not finish
            t = table(refinements, iterCGopt, timeCGopt, errorCGopt, solvingErrorCGopt)
            writetable(t, ...
                strcat('data/3Dexample', num2str(example), '-CG-opt-exact-', ...
                num2str(exactFlag), '-maxIt-', num2str(maxIt), ...
                '-tolerance-', num2str(tolerance)),'Delimiter',' ')
            
        case 2
            % Lyap operator preconditioner
            fprintf('Starting CG lyaponov\n')
            problem.precond='lyap';
            tt=tic;
            [U, iterCGlyap(refinement)]=pcg_fun4(funA,rhsfull,...
                0*rhsfull,problem,maxIt,tolerance, ...
                size(problem.M_space,1), ...
                size(problem.M_time,2),exactFlag);
            timeCGlyap(refinement) = toc(tt);
            U_cg_lyap=U(:)*norm(problem.rhs(:));
            
            solvingErrorCGlyap(refinement) = ...
                norm(funB(U_cg_lyap) - reshape(problem.rhs, [numel(problem.rhs), 1])) / norm(reshape(problem.rhs, [numel(problem.rhs), 1]));
            
            fprintf('Refinement: %d, Time to solve: %f\n', ...
                refinement,  timeCGlyap(refinement))
            
            fprintf('Testing the cg_lyap solution\n');
            sol = get3Dsolution(problem, U_cg_lyap, resolution);
            
            errorCGlyap(refinement) = calculate3DL2Error(problem, sol);
            
            % Write everything in a file in case the whole script does
            % not finish
            t = table(refinements,iterCGlyap, timeCGlyap, errorCGlyap, solvingErrorCGlyap)
            writetable(t, ...
                strcat('data/3Dexample', num2str(example), '-CG-lyap-exact-', ...
                num2str(exactFlag), '-maxIt-', num2str(maxIt), ...
                '-tolerance-', num2str(tolerance)),'Delimiter',' ')
        case 3
            % Galerkin
            fprintf('Starting Galerkin\n')
            [uu,ss,vv]=svds(problem.rhs,1);
            rhs1=uu(:,1)*sqrt(ss(1,1));
            rhs2=vv(:,1)*sqrt(ss(1,1));
            % path(path,'./valeria/')
            info=1;
            tt=tic;
            [X1,X2,restot,iterGalerkin(refinement)]=Galerkin3(problem.M_space,...
                2*problem.A_space,problem.Q_space,problem.Q_time, ...
                (problem.D_time+problem.D_time')/2,...
                problem.M_time,rhs1,rhs2,maxIt,tolerance,info);
            timeGalerkin(refinement) = toc(tt);
            U=X1*X2';
            U_galerkin=U(:);
            
            solvingErrorGalerkin(refinement) = ...
                norm(funB(U_galerkin) - reshape(problem.rhs, [numel(problem.rhs), 1])) / norm(reshape(problem.rhs, [numel(problem.rhs), 1]));
            
            
            fprintf('Galerkin: Refinement: %d, Time to solve: %f\n', ...
                refinement,  timeGalerkin(refinement))
            
            fprintf('Testing the galerkin solution\n');
            sol = get3Dsolution(problem, U_galerkin, resolution);
            
            errorGalerkin(refinement) = calculate3DL2Error(problem, sol);
            
            % Write everything in a file in case the whole script does
            % not finish
            t = table(refinements,iterGalerkin, timeGalerkin, errorGalerkin, solvingErrorGalerkin)
            writetable(t, ...
                strcat('data/3Dexample', num2str(example), '-Galerkin-', ...
                num2str(exactFlag), '-maxIt-', num2str(maxIt), ...
                '-tolerance-', num2str(tolerance)),'Delimiter',' ')
            
         case 4
            % GalerkinNEW
            fprintf('Starting Galerkin\n')
            [uu,ss,vv]=svds(problem.rhs,1);
            rhs1=uu(:,1)*sqrt(ss(1,1));
            rhs2=vv(:,1)*sqrt(ss(1,1));
            % path(path,'./valeria/')
            info=1;
            tt=tic;
            [X1,X2,restot,iterGalerkin(refinement)]=Galerkin3NEW(problem.M_space,...
                2*problem.A_space,problem.Q_space,problem.Q_time, ...
                (problem.D_time+problem.D_time')/2,...
                problem.M_time,rhs1,rhs2,maxIt,tolerance,tolerance,info);
            timeGalerkin(refinement) = toc(tt);
            U=X1*X2';
            U_galerkin=U(:);
            
            solvingErrorGalerkin(refinement) = ...
                norm(funB(U_galerkin) - reshape(problem.rhs, [numel(problem.rhs), 1])) / norm(reshape(problem.rhs, [numel(problem.rhs), 1]));
            
            
            fprintf('Galerkin: Refinement: %d, Time to solve: %f\n', ...
                refinement,  timeGalerkin(refinement))
            
            fprintf('Testing the galerkin solution\n');
            sol = get3Dsolution(problem, U_galerkin, resolution);
            
            errorGalerkin(refinement) = calculate3DL2Error(problem, sol);
            
            % Write everything in a file in case the whole script does
            % not finish
            t = table(refinements,iterGalerkin, timeGalerkin, errorGalerkin, solvingErrorGalerkin)
            writetable(t, ...
                strcat('data/3Dexample', num2str(example), '-GalerkinNEW-', ...
                num2str(exactFlag), '-maxIt-', num2str(maxIt), ...
                '-tolerance-', num2str(tolerance)),'Delimiter',' ')
                
    end
end
end