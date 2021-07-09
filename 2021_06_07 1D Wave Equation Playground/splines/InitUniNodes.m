% By C. Mollet
%
% Function to initialize set of nodes for B-Splines of order k
% on a uniform grid
%
% Input:
% ------
% 'j' : Level of refinement
% 'k' : Order of B-Splines
%
% Output:
% -------
% T : Corresponding nodes

function T = InitUniNodes(j,k)

% Setting number of Nodes (for test- and solutionspace)
n = 2^j+k-1;

% Init Nodes ("Full" w/o boundary offset )
T = zeros(n+k,1);
for i=k+1:n
    T(i)=(i-k) * 2^(-j);
end
T(n+1:n+k)=ones(k,1);

end