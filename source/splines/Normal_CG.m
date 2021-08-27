% By C. Mollet
%
% Function to solve normal equations with intermediate Riesz mapping, i.e.,
% B'*RY^(-1)*B*u = B'*RY^(-1)*f
% by using an outer cg-loop and avoiding the explicit
% assembling of the matrix B'*RY^(-1)*B and therefore saves
% massively memory usage or/and computation time.
% Because: Since RY is sparse, its inverse do not necessarily
% need to be sparse, but even if so, the explicit invertation would be
% much to expansive for large matrices and an approximation (RY\B) would be
% fully populated.
% To this end, this functions only acts on the input vecor and right hand
% side and circumvent the explicit assembling of the matrix RY^(-1)*B
% in this way.
%
% input:
% ------
% 'B'       : System matrix
% 'RY'      : Gram matrix w.r.t. test space
% 'f'       : Right hand side (not in "minimal residuum form")
% 'kmax'    : Maximum number of iterations
% 'tol'     : Required accuracy of relative residual
% 'x'       : Initial Vector
%
% output:
% -------
% 'x'       : Solution vector
% 'counter' : Number of iterations
% 'res'     : Relative residuum

function [ x , counter , res ] = Normal_CG( B , RY , f , kmax, tol, x )

% Transfer right hand side to "minimal residuum form"
% and calculate initial residuum
g = B'* (RY\f);
r = B'*(RY\ (f-B*x) );
d = r;
ng = norm(g,2);

counter = 0;
res = norm(r,2);
for k=1:kmax
    
    % Check if required accuracy is achieved and also
    % capture the trivial case of right hand side zero
    if (norm(r,2)/ng < tol) || (norm(r,2) < eps*eps)
        break
    end
    % CG-iteration; instead of applying it directly
    % to B'*RY^(-1)*B, we avoid the calculation of the
    % inverse of RY here
    z = B'*(RY\(B*d));
    alpha = (r'*r)/(d'*z);
    x = x + alpha*d;
    r_old = r;
    r = r - alpha*z;
    beta = (r'*r)/(r_old'*r_old);
    d = r + beta*d;
    counter = counter + 1;
    res = norm(r,2)/ng;
end
end