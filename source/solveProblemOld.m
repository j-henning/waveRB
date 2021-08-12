
function U = solveProblemOld(problem)


% Use Galerkin for solving
[uu,ss,vv]=svds(problem.rhs,1);
rhs1=uu(:,1)*sqrt(ss(1,1));
rhs2=vv(:,1)*sqrt(ss(1,1));
tolG=1e-8;
maxIt = 200;
info = 0;
[X1,X2]= ...
    Galerkin3(problem.M_space,2*problem.A_space,problem.Q_space,problem.Q_time,...
    (problem.D_time+problem.D_time')/2,problem.M_time,rhs1,rhs2,maxIt,tolG,1e-6,info);
U=X1*X2';
U=U(:);

% resolution.x = 8;
% resolution.t = 8;
% sol = get1Dsolution(problem, U, resolution);

end