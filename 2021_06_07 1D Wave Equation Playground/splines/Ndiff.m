% By C. Mollet
%
% Helpfunction to evaluate B-Spline derivatives N^(diff)_{i,k,T}(x)
%
% Input:
% ------
% 'T'     : Set of extended knots
% 'k'     : B-Spline Oder
% 'diff'  : Differntial order
% 'i'     : Position index
% 'x'     : Evaluation point
%
% Output:
% -------
% 'Nx'    : Value of N^(diff)_{i,k,T}(x)

function Nx = Ndiff(T,k,diff,i,x)
% Recursively defined via a general rule of derivatives of
% B-Splines, combined with a Neville scheme for the point evaluation

% Smoothness: N_{i,k,T} \subset H^{k-1}
if k <= diff
  error('Spline not sufficiently smooth!');
end

if diff == 0
    % Case without derivative
    % Neville scheme with a delta sequence as coefficient vector
    Nx = Nev(sparse(1,i,1,1,length(T)-k),{T},k,x);
    
else
    % Formular for derivatives
    Nx = (k-1)*( Ndiff(T,k-1,diff-1,i,x) * Fac(T(i+k-1),T(i)) - Ndiff(T,k-1,diff-1,i+1,x) * Fac(T(i+k),T(i+1)) );
    
end



end

% Simple helpfunction for coincident nodes
function value = Fac(a,b)

if a~=b
    value = 1/(a-b);
else
    value = 0;
end

end