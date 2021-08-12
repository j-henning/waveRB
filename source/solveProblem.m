function [U] = solveProblem(problem, method, maxIt, tolerance1, tolerance2, info)
if nargin < 2
    method = 'galerkin';
end
if nargin < 3
    maxIt = 100;
end
if nargin < 4
    tolerance1 = 1e-10;
end
if nargin < 5
    tolerance2 = 1e-4;
end
if nargin < 6
    info = 0;
end

switch method
    case 'galerkin'
        % Use Galerkin for solving
        [uu,ss,vv]=svds(problem.rhs,1);
        rhs1=uu(:,1)*sqrt(ss(1,1));
        rhs2=vv(:,1)*sqrt(ss(1,1));

        [X1,X2]= ...
            Galerkin3(problem.M_space,2*problem.A_space,problem.Q_space,problem.Q_time,...
            (problem.D_time+problem.D_time')/2,problem.M_time,rhs1,rhs2,maxIt,tolerance1,tolerance2,info);
        U=X1*X2';
        U=U(:);
        
        
end



end