% By C. Mollet
%
% Helpfunction to Assemble Riesz Matrices wrt B-Splines (using StiffMat)
%
% input:
% ------
% 'T'      : Set of extended knots for B-spline basis
% 'k'      : Order of B-splines
% 'offset' : Offset due to boundary conditions
% 'order'  : Order of differntiability of space H^order
%
% output:
% -------
% 'R'      : Riesz Matrix


function [ R ] = RieszMat( T,k,offset, order )

% Assembling Gram matrices on primal Sobolev spaces (positive order)
if order >= 0
    R = StiffMat(T,T,k,k,offset,offset,[0,0]);
    for i = 1:order
        R = R + StiffMat(T,T,k,k,offset,offset,[i,i]);
    end
else
    % Assembling discrete Riesz discretizations on dual Sobolev spaces 
    %(positive order), i.e., R_L2*R_V^(-1)*R_L2 for ansatz functions in L2
    P = StiffMat(T,T,k,k,offset,offset,[0,0]);
    R = P;
    for i = 1:abs(order)
        R = R + StiffMat(T,T,k,k,offset,offset,[i,i]);
    end
    R = P * (R\P);
end


end

