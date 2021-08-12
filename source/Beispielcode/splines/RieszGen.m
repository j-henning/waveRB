% By C. Mollet
%
% Helpfunction to Assemble Riesz Matrices wrt B-Splines (using StiffMat)
%
% input:
% ------
% 'T'       : Cell array of extended knot sequences for B-splines
% 'k'       : Order of B-Splines
% 'offset'  : 2D-Vector Offset due to boundary conditions, e.g. value 1
%             for zero boundary condition (first and last B-spline
%             will be excluded)
% 'Matind'  : (Cell) Indices correpsonding to appearing derivatives;
%             Matind(d,t) encodes the order of Sobolev space w.r.t. the d-th.
%             coordinate and t-th. intersection
%
% output:
% -------
% 'R'      : Gram Matrix (discretized Riesz mapping)

function [ R ] = RieszGen( T , k , offset , Matind )

[dim,terms] = size(Matind);

% Calculating R
% First term init
A = RieszMat(T{dim}, k(dim) , offset(dim,:) , Matind(dim,1));
for i=dim-1:-1:1
    A = kron( RieszMat(T{i}, k(i) , offset(i,:) , Matind(i,1)) , A );
end
R = A;

% loop
for n=2:terms
    A = RieszMat(T{dim}, k(dim) , offset(dim,:) , Matind(dim,n));
    for i=dim-1:-1:1
        A = kron( RieszMat(T{i}, k(i) , offset(i,:) , Matind(i,n)) , A );
    end
    R = R + A;
end

end

