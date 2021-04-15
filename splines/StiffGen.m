% By C. Mollet
%
% Function to generate sum of tensor product stiffness matrices
%
% input:
% ------
% 'Tsol, Ttest'           : Cell array of extended knot sequences for B-splines
%                           wrt to solution space and test space
% 'ksol, ktest'           : Order of B-Splines in solution- and testspace
% 'Maind'                 : Classification of matrices arranged via StiffMat
% 'offsetsol, offsettest' : 2D-Vector Offset due to boundary condition, e.g. value 1
%                           for zero boundary condition (first and last B-spline
%                           will be excluded)
% 'secForm'               : Optional Factor to handle 2nd form of parabolic full weak form
%                           (Helpconstant)
%
% output:
% -------
% 'B'      : System Matrix according to input

function [ B ] = StiffGen( Tsol,Ttest,ksol,ktest,offsetsol,offsettest,Matind, secForm )

% Setting sign according to first or second weak formulation of parabolic
% PDEs
if nargin < 8
    c = 1.0;
elseif secForm == 1
    c = -1.0;
elseif secForm == 0
    c = 1.0;
else
    error('Wrong Specification of Form');
end
[dim,terms] = size(Matind);

% First term init
A = StiffMat(Tsol{dim},Ttest{dim},ksol(dim),ktest(dim),offsetsol(dim,:),offsettest(dim,:),Matind{dim,1}');
for i=dim-1:-1:1
    A = kron( StiffMat(Tsol{i},Ttest{i},ksol(i),ktest(i),offsetsol(i,:),offsettest(i,:),Matind{i,1}') , A );
end
B = c * A;

% loop over remaining terms
for n=2:terms
    A = StiffMat(Tsol{dim},Ttest{dim},ksol(dim),ktest(dim),offsetsol(dim,:),offsettest(dim,:),Matind{dim,n}');
    for i=dim-1:-1:1
        A = kron( StiffMat(Tsol{i},Ttest{i},ksol(i),ktest(i),offsetsol(i,:),offsettest(i,:),Matind{i,n}') , A );
    end
    B = B + A;
end


end

