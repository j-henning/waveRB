clear
close all
clc

addpath('../splines');
addpath('../splines/Utilities');
addpath('../solvers');



%% Specify all needed input data
problemConfiguration.bSplineOrder_time = 3; %3
problemConfiguration.bSplineOrder_space = 3; %3

problemConfiguration.offset_time_ansatz = [0, 2]; %[0, 2]
problemConfiguration.offset_time_test = [0, 2]; %[0, 2]
problemConfiguration.offset_space_ansatz = [1, 1]; %[1, 1];
problemConfiguration.offset_space_test = [1, 1]; %[1, 1];

problemConfiguration.d = 3; % Dimension

funR = @(x,y,z) sqrt((x-0.5).^2 + (y-0.5).^2 + (z-0.5).^2);  % Returns the radius

smoothing = false;

refinements = 1:4;
method=3; % Specify which methods should be tested (CG Opt, CG Lyap, Galerkin)

%% Example 1 continuous rhs
% problemConfiguration.f_time = {@(t) 0*t};
% problemConfiguration.u0 = @(x,y,z) (1-5*funR(x,y,z)).*( funR(x,y,z) <= 0.2);
% u0r = @(r) (1-5*abs(r)) .* (abs(r) <= 0.2) .* (abs(r) >= 0);
% problemConfiguration.u1 = @(x,y,z)  0.*x.*y.*z;
% problemConfiguration.inital_conditions = true;
% problemConfiguration.c = .2;

%% Example 2 discontinuous rhs
problemConfiguration.f_time = {@(t) 0*t};
problemConfiguration.u0 = @(x,y,z) ( funR(x,y,z) <= 0.2);
u0r = @(r) (abs(r) <= 0.2);
problemConfiguration.u1 = @(x,y,z)  0.*x.*y.*z;
problemConfiguration.inital_conditions = true;
problemConfiguration.c = .1;





% Specify the solution for plotting and the error calculation
% Should be higher than the number of refinements

resolution.x = max(6, refinements(end)+1);
resolution.y = max(6, refinements(end)+1);
resolution.z = max(6, refinements(end)+1);
resolution.t = max(6, refinements(end)+1);

% Define the resolution for the L2 error calculation (and plotting)
x = linspace(0,1, 2^resolution.x+1);
y = linspace(0,1, 2^resolution.y+1);
z = linspace(0,1, 2^resolution.z+1);

% Compute the analytical solution
t = linspace(0,1,2^resolution.t+1);
[X, Y, Z, T] = ndgrid(x, y, z, t);


for i=1:2^resolution.x+1
    for j=1:2^resolution.y+1
        for k=1:2^resolution.z+1
            r = funR(x(i), y(j), z(k));
            sol_ref(i,j,k,:) = dAlemenbertSphere(r,t,problemConfiguration.c, u0r);
        end
    end
end

tolerance = 1e-12;
maxIt = 1000;
exactFlag = false;


for refinementLevel_space = refinements
    for refinementLevel_time = refinementLevel_space
        
        problemConfiguration.refinementLevel_space = refinementLevel_space;
        problemConfiguration.refinementLevel_time = refinementLevel_time;
        
        fprintf('Creating problem...')
        problem = create3DWaveProblemImproved(problemConfiguration);
        fprintf(' Done!\n')
        
        % Specify a function handle for the pcg methods
        funA=@(X)( problem.M_space * X * problem.Q_time' ...
            + problem.A_space' * X * problem.D_time' ...
            + problem.A_space * X * problem.D_time ...
            + problem.Q_space * X * problem.M_time');
        
        
        
        %% Sove the linear equation system
        
        rhsfull=reshape(problem.rhs, [numel(problem.rhs) 1]);
        rhsfull=rhsfull/norm(rhsfull);
        
        switch method
            case 1
                % Product operator preconditioner
                
                problem.precond='optimal';
                tt=tic;
                [U, iterCGopt(refinementLevel_space, refinementLevel_time)]= ...
                    pcg_fun4(funA,rhsfull,0*rhsfull,problem,maxIt,tolerance, 1e-2,...
                    size(problem.M_space,1),size(problem.M_time,2),exactFlag);
                U_cg_optimal=U(:)*norm(problem.rhs(:));
                timeCGopt(refinementLevel_space, refinementLevel_time) = toc(tt);
                
                fprintf('Time refinement: %d, Space refinement: %d, Time to solve: %f\n', ...
                    refinementLevel_time, refinementLevel_space,  timeCGopt(refinementLevel_space, refinementLevel_time))
                
            case 2
                % Lyap operator preconditioner
                problem.precond='lyap';
                tt=tic;
                [U, iterCGlyap(refinementLevel_space, refinementLevel_time)]=...
                    pcg_fun4(funA,rhsfull,0*rhsfull,problem,maxIt,tolerance,1e-2,size(problem.M_space,1),size(problem.M_time,2),exactFlag);
                timeCGlyap(refinementLevel_space, refinementLevel_time) = toc(tt);
                U_cg_lyap=U(:)*norm(problem.rhs(:));
                
                fprintf('Time refinement: %d, Space refinement: %d, Time to solve: %f\n', ...
                    refinementLevel_time, refinementLevel_space,  timeCGlyap(refinementLevel_space, refinementLevel_time))
                
            case 3
                % Galerkin
                [uu,ss,vv]=svds(problem.rhs,1);
                rhs1=uu(:,1)*sqrt(ss(1,1));
                rhs2=vv(:,1)*sqrt(ss(1,1));
                
                
                %      path(path,'./valeria/')
                
                info=1;
                tt=tic;
                [X1,X2,restot,iterGalerkin(refinementLevel_space, refinementLevel_time)]= ...
                    Galerkin3(problem.M_space,2*problem.A_space, ...
                    problem.Q_space,problem.Q_time, ...
                    (problem.D_time+problem.D_time')/2, ...
                    problem.M_time,rhs1,rhs2,maxIt,tolerance,1e-2,info);
                timeGalerkin(refinementLevel_space, refinementLevel_time) = toc(tt);
                U=X1*X2';
                U_galerkin = U(:);
                
                fprintf('Galerkin: Time refinement: %d, Space refinement: %d, Time to solve: %f\n', ...
                    refinementLevel_time, refinementLevel_space,  timeGalerkin(refinementLevel_space, refinementLevel_time))
                
            case 4
                B = kron(problem.Q_time, problem.M_space) ...
                    + kron(problem.D_time, problem.A_space')...
                    + kron(problem.D_time', problem.A_space)...
                    + kron(problem.M_time, problem.Q_space);
                
                U = B\problem.rhs(:);
        end
        
        
        switch method
            case 1
                U = U_cg_optimal;
                fprintf('Testing the cg_opt solution\n');
            case 2
                U = U_cg_lyap;
                fprintf('Testing the cg_lyap solution\n');
            case 3
                U = U_galerkin;
                fprintf('Testing the galerkin solution\n');
            case 4
                fprintf('Testing the backslash solution\n');
        end
        
        %% Get the solution
        
        
        % Smoothing
        
        if smoothing
            res.x = refinementLevel_space;
            res.y = refinementLevel_space;
            res.z = refinementLevel_space;
            res.t = refinementLevel_time;
            
            sol = get3Dsolution(problem, U, res, true);
            
            
            
            xCoarse = linspace(0,1,2^res.x+1);
            yCoarse = linspace(0,1,2^res.y+1);
            zCoarse = linspace(0,1,2^res.z+1);
            tCoarse = linspace(0,1,2^res.t+1);
            
            solSmooth = zeros(size(sol,1)+1, ...
                size(sol,2)+1, ...
                size(sol,3)+1, ...
                size(sol,4)+1);
            
            for i = 2:size(solSmooth,1)-1
                for j = 2:size(solSmooth,2)-1
                    for k = 2:size(solSmooth,3)-1
                        for l = 2:size(solSmooth,4)-1
                            solSmooth(i,j,k,l) = 1/16 * ...
                                (sol(i-1,j-1,k-1,l-1) ...
                                + sol(i-1,j-1,k-1,l) ...
                                + sol(i-1,j-1,k,l-1) ...
                                + sol(i-1,j-1,k,l) ...
                                + sol(i-1,j,k-1,l-1) ...
                                + sol(i-1,j,k-1,l) ...
                                + sol(i-1,j,k,l-1) ...
                                + sol(i-1,j,k,l) ...
                                + sol(i,j-1,k-1,l-1) ...
                                + sol(i,j-1,k-1,l) ...
                                + sol(i,j-1,k,l-1) ...
                                + sol(i,j-1,k,l) ...
                                + sol(i,j,k-1,l-1) ...
                                + sol(i,j,k-1,l) ...
                                + sol(i,j,k,l-1) ...
                                + sol(i,j,k,l));
                        end
                    end
                end
            end
            
            % Interp2 does not extrapolate, so we have to do it
            % by hand before we feed the data into the method
            solSmooth(:,:,:,1) = solSmooth(:,:,:,2) + (solSmooth(:,:,:,2) - solSmooth(:,:,:,3));
            %
            solSmooth(:,:,:,end) = solSmooth(:,:,:,end-1) + (solSmooth(:,:,:,end-1) - solSmooth(:,:,:,end-2));
            
            [X, Y, Z, T] = ndgrid(xCoarse, yCoarse, zCoarse, tCoarse);
            [Xfine, Yfine, Zfine, Tfine] = ndgrid(x,y,z,t);
            sol = interpn(X,Y,Z,T,solSmooth,Xfine,Yfine,Zfine,Tfine, 'linear');
        else
            sol = get3Dsolution(problem, U, resolution, false);
        end
        
        
        
        
        errorL2(refinementLevel_space) = sqrt(mean((sol-sol_ref).^2, 'all','omitnan'));
        
        
        
    end
end

err = errorL2;

conv_rate = log(err(end)/ err(end-1)) / log(1/2)



%% Plot the solution

% f = figure('WindowState','maximized');
% for i=1:size(sol_ref,4)
%     subplot(2,2,1)
%     s = surf(squeeze(sol(:,:,floor(0.5 * 2^resolution.z), i))); s.EdgeColor = 'none';
%     title(num2str(i/size(sol_ref,4)));
%
%     subplot(2,2,2)
%     s = surf(squeeze(sol_ref(:,:,floor(0.5 * 2^resolution.z), i))); s.EdgeColor = 'none';
%     title(num2str(i/size(sol_ref,4)));
%
%     subplot(2,2,[3 4])
%     s = surf(squeeze(sol(:,:,floor(0.5 * 2^resolution.z), i)) - squeeze(sol_ref(:,:,floor(0.5 * 2^resolution.z), i))); s.EdgeColor = 'none';
%     title(num2str(i/size(sol_ref,4)));
%     drawnow
%     pause(0.1)
%
%
%
%
% end


% Produces NaN for r = 0
function sol = dAlemenbertSphere(r,t,c, u0)
if r == 0
    sol = NaN;
else
    sol = ( (r+c*t) .* u0(r + c*t) + (r - c*t).*u0(r-c*t))./(2.*r) ;
end
end


