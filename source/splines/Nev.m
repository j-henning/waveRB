% By C. Mollet
%
% Function to Perform a tensor product Neville-like de Boor scheme
%
% Input:
% ------
% 'c'     : n-D vectorized expansion coefficients
% 'T'     : Set of extended knots of B-splines given as a cell element {T1,T2,...,Tn}
% 'k'     : Order of B-splines in each direction
% 'x'     : Evaluation point (n-D)
%
% Output:
% -------
% 'Val' : Value of tensorized spline with expanion coefficients 'c' at point 'x'

function Val = Nev(c,T,k,x)

% Set number of elements
dim = length(T);

n = zeros(dim,1);
for i = 1:dim
    n(i) = length(T{i})-k(i);
end

if dim == 1
    
    %% Implementation of the 1D Neville scheme %%
    
    T1d = T{1};
    
    % initialize helpvector detemined by 'c'
    d = zeros(length(c),k);
    d(:,1)=c;
    
    % capture the case where 'x' does not belong to the global interval
    sizeT = numel(T1d);
    if (x<T1d(1)) || (x>=T1d(sizeT))
        Val = 0;
        return;
    end
    
    % Determine the subinterval in which 'x' lies
    for j=k:sizeT-k
        if (x < T1d(j+1))
            i=j;
            break;
        end
    end
    
    % Applying Neville scheme
    for r=1:(k-1)
        for s=i-(k-1)+r:i
            % Determination of c_s^[r]=d(s,r+1):
            d(s,r+1)=((x-T1d(s))/(T1d(s+k-r)-T1d(s)))*d(s,r) + ((T1d(s+k-r)-x)/(T1d(s+k-r)-T1d(s)))*d(s-1,r);
            % Note: Case T1d(s+k-r)==T2d(s) cannot appear in a reasonable setting
        end
    end
    Val = d(i,k);
    
else
    % initialize 'outer' expansion coefficients
    d = zeros(1,n(1));
    l = prod(n(2:dim));
    for i=1:n(1)
        % recursive inner Neville scheme to detemine 'outer'
        % coefficients/functions
        d(i)= Nev( c(1+(i-1)*l:i*l) ,T(2:dim),k(2:dim),x(2:dim));
    end
    % outer 1D-Neville scheme
    Val = Nev(d,T(1),k(1),x(1));
    
end


end