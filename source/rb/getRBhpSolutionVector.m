% Return the coefficient vector for a given binary tree and a parameter mu
% and a problem p
% Input:
% tree         - Reduced Basis model
% mu           - Query parameter
% Output:
% U            - Solution vector
function U = getRBhpSolutionVector(tree, mu)
index = findInterval(tree, mu);
% p = changeWaveSpeed(p, mu);

% TODO: Calculate the matrix products in the offline phase and just scale
% and and it up through in the online phase

% B = tree{index}.Y' * (coefficients{1}(mu)  * matrices{1}) * tree{index}.Y;
% for i = 2:length(matrices)
% B = B +  tree{index}.Y' * (coefficients{i}(mu)  * matrices{i}) * tree{index}.Y; 
% end
B = tree{index}.coefficients{1}(mu)  * tree{index}.matrices{1};
for i = 2:length(tree{index}.matrices)
B = B + tree{index}.coefficients{i}(mu) *  tree{index}.matrices{i}; 
end

f = tree{index}.f;

U = B \ f;
U = tree{index}.Y * U;
end