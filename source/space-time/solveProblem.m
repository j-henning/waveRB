% Solves the linear system of a space time problem
% Input:
% problem    - Problem struct
% method     - 'Galerkin' (default), 'cg-lyap', or 'cg-opt'
% maxIt      - Maximum number of iterations (default: 1000)
% tolerance1 - Backward error tolerance (default 1e-5)
% tolerance2 - Residual error tolerance (default 1e-2)
% info       - Display info during the solving procedure (default 0)
% Output:
% U          - Solution vector
% iterations - Number of iterations
function [U, iterations] = solveProblem(problem, method, maxIt, tolerance1, tolerance2, info)
if nargin < 2
    method = 'galerkin';
end
if nargin < 3
    maxIt = 1000;
end
if nargin < 4
    tolerance1 = 1e-5; 
end
if nargin < 5
    tolerance2 = 1e-2;
end
if nargin < 6
    info = 0;
end

funA=@(X)( problem.M_space * X * problem.Q_time' ...
    + problem.A_space' * X * problem.D_time' ...
    + problem.A_space * X * problem.D_time ...
    + problem.Q_space * X * problem.M_time');

switch method
    case 'galerkin'
        % Use Galerkin for solving
        [uu,ss,vv]=svds(problem.rhs,1);
        rhs1=uu(:,1)*sqrt(ss(1,1));
        rhs2=vv(:,1)*sqrt(ss(1,1));
        
        [X1,X2, ~, iterations]= ...
            Galerkin(problem.M_space,2*problem.A_space,problem.Q_space,problem.Q_time,...
            (problem.D_time+problem.D_time')/2,problem.M_time,rhs1,rhs2,maxIt,tolerance1,tolerance2,info);
        U=X1*X2';
        U=U(:);
        
    case 'cg-lyap'
        problem.precond='lyap';
        rhsfull=reshape(problem.rhs, [numel(problem.rhs) 1]);
        rhsfull=rhsfull/norm(rhsfull);
        [U, iterations] = pcg_fun(funA,rhsfull,0*rhsfull,problem,maxIt,tolerance1, tolerance2,...
            size(problem.M_space,1),size(problem.M_time,2),0);
        U = U(:)*norm(problem.rhs(:));
    case 'cg-opt'
        problem.precond='optimal';
        rhsfull=reshape(problem.rhs, [numel(problem.rhs) 1]);
        rhsfull=rhsfull/norm(rhsfull);
        [U, iterations] = pcg_fun(funA,rhsfull,0*rhsfull,problem,maxIt,tolerance1, tolerance2,...
            size(problem.M_space,1),size(problem.M_time,2),0);
        U = U(:)*norm(problem.rhs(:));
        
    case 'backslash'
        B = kron(problem.Q_time, problem.M_space) ...
            + kron(problem.D_time, problem.A_space') ...
            + kron(problem.D_time', problem.A_space) ...
            + kron(problem.M_time, problem.Q_space);
        
        U = B \ problem.rhs(:);
        
end
end
