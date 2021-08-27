% By C. Mollet

% Function to set up system matrices of Petrov-Galerkin discretizations
% with B-splines on unit interval [0,1]. Not restricted to uniform
% discretiazion

% input:
% ------
% 'Tsol, Ttest'           : Set of extended knots for B-spline basis wrt to
%                           Solution and test space
% 'ksol, ktest'           : Order of B-splines in solution and test space
% 'diff'                  : Classification; diff=(d1,d2) yields
%                           <phi^(d1),phi^(d2)> , where phi^(d1) is the
%                           d1's derivative. E.g. diff=(1,1) -> Laplace and
%                           diff=(0,0) -> Massmatrix
% 'offsetsol, offsettest' : Offset due to boundary condition, e.g. value 1
%                           for zero boundary condition (first and/or last B-spline
%                           will be excluded)
%                           offsetsol(1) left boundary solution space
%                           offsetsol(2) right boundary solution space
%                           offsettest(1) left boundary test space
%                           offsettest(2) right boundary test space
%
% output:
% -------
% 'A'      : System matrix according to input (cf. 'diff')

function A = StiffMat(Tsol,Ttest,ksol,ktest,offsetsol,offsettest,diff)

% Determine number of (full) B-splines/dimension
ntest = length(Ttest)-ktest;
nsol = length(Tsol)-ksol;

% Determine required Order for quadrature to obtain exact solution
% Gauss-Quadrature is exact up to order 2k+1
order = ceil((ktest+ksol - 1)/2);

% Setting up system matrix
A = sparse(ntest-(offsettest(1)+offsettest(2)),nsol-(offsetsol(1)+offsetsol(2)));

% Preallocate relevant indices which overlap
valsol = 1+offsetsol(1):nsol-offsetsol(2);
valtest = 1+offsettest(1):ntest-offsettest(2);
[Vsol,Vtest] = meshgrid(valsol,valtest);
[s,t] = find ((Tsol(Vsol + ksol) > Ttest(Vtest)) & (Tsol(Vsol) < Ttest(Vtest+ktest)) );
for it = [s,t]'
    i = Vtest(it(1),it(2));
    l = Vsol(it(1),it(2));
    

    % Integrand at nodes 'l' and 'i'
        f = @(x) Ndiff(Tsol,ksol,diff(1),l,x)...
            * Ndiff(Ttest,ktest,diff(2),i,x);
    % Quadrature on each polynomial part of the B-splines/Integrand
    
    % Only iteration over the intersection of the supports
    % needed
    left = max(Ttest(i),Tsol(l));
    right = min(Ttest(i+ktest),Tsol(l+ksol));
    T = [Ttest(i:i+ktest);Tsol(l:l+ksol)];
    ind = (left<=T) & (right>=T) ;
    T = unique(T(ind));
    for step = 1:length(T)-1
        A(i-offsettest(1),l-offsetsol(1)) = A(i-offsettest(1),l-offsetsol(1)) + Gaussq(T(step),T(step+1),f,order);
    end
end

end