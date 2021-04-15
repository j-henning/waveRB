% By C. Mollet
%
% Assembles a Refinement Matrix of hierachical B-spline bases.
% Only works for unirom grids.
%

function [ M ] = RefineMat( lev , k , offset )

if nargin == 2
    offset = [0,0];
end

% Refined Mesh
T = InitUniNodes(lev+1,k);

% Coarse Grid
T_coarse = InitUniNodes(lev,k);

% Refinement Matrix
n = numel(T)-k;
n_coarse = numel(T_coarse)-k;
Mv = sparse(n,n_coarse);
% boundary adaptation (left)
% Notice: Becomes unstable for high orders (Matrix badly scaled)
for i=1:k-1
    N = zeros(i+1);
    v = zeros(i+1,1);
    for a=1:i+1
        for b=1:i+1
            N(a,b) = Ndiff(T,k,0,i+b-1,T(a+k-1+i-1));
        end
        v(a) = Ndiff(T_coarse,k,0,i,T(a+k-1+i-1));
    end
    Mv(i:2*i,i) = N \ v;
end

% Inner Refinements
a = zeros(k+1,1);
for l = 0:k
    a(l+1) = 2^(1-k) * nchoosek(k,l);
end

ind = k;
for j = k:n_coarse-k+1
    Mv(ind:ind+k,j) = a;
    ind = ind + 2;
end

count = 1;
% Mirrored to right boundary
for j = n_coarse:-1:n_coarse-k+2
    Mv(:,j) = Mv(n:-1:1,count);
    count = count + 1;
end

%M = Mv;
M = Mv(1+offset(1):n-offset(2),1+offset(1):n_coarse-offset(2));


end

