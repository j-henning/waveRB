% By C.Mollet
%
% !! This is only a very first attempt and not testet or optimized !!
% Multigrid Method to solve Petrov-Galerkin solutions
% using B-splines

function [ uk ] = MGM( A, f, u, l, ksol, ktest, offsetsol, offsettest , nu1, nu2, mu )

% Todo: 
% *) Implement a version for minimal residual Petrov-Galerkin
% *) Optimize and other smoother

% A-priori smoothing:
D = sparse(diag(diag(A)));
L = sparse(tril(A,-1));
R = sparse(triu(A,1));
%I = speye(length(A));

uk = u;
for i=1:nu1
    %uk = 0.5 * D^(-1)*f + 0.5 * (I - D^(-1)*(L+R))*uk;
    %uk = -(D+L)^(-1) * R * uk + (D+L)^(-1) * f;
    [L,R] = ilu(A,struct('type','ilutp','droptol',1e-6 ));
    uk = gmres(A , f , 5 , 1e-1 * 2^( -max(ksol) * max(l) ) , size(A,1), L , R,uk);
end

% Corse grod correction:
M1 = RefineMat(l(1)-1,ksol(1),offsetsol(1,:));
M2 = RefineMat(l(2)-1,ksol(2),offsetsol(2,:));
p = kron(M1,M2);

M1T = RefineMat(l(1)-1,ktest(1),offsettest(1,:));
M2T = RefineMat(l(2)-1,ktest(2),offsettest(2,:));
r = kron(M1T',M2T');


B = r * A * p;
g = r * (f - A * uk);

v = zeros(size(M1,2)*size(M2,2),1);

if( min(l)==2 )
    %R = chol(A);
    %v = inv(R') * inv(R) * g;
    v = B\g;
else
    for i=1:mu
        %v = MGM(l-1, B, g, v, nu1, nu2, mu);
        v = MGM(B,g,v,l-1,ksol,ktest,offsetsol,offsettest,nu1,nu2,mu);
    end
end
uk = uk + r' * v;

% A-Posteriori Smoothing:
for i=1:nu2
    %uk = 0.5 * D^(-1)*f + 0.5 * (I - D^(-1)*(L+R))*uk;
    uk = -(D+L)^(-1) * R * uk + (D+L)^(-1) * f;
end


end

