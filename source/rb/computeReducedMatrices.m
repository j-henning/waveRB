% Compute and store the reduced affine matrices as well as the coefficient
% functions for each subinterval
function [tree] = computeReducedMatrices(tree, problem, matrices, coefficients)
tree = computeReducedMatricesRecursive(tree, 1, problem, matrices, coefficients);
end


% Recurisive version of computeReducedMatrices
function [tree] = computeReducedMatricesRecursive(tree, index, problem, matrices, coefficients)
leftChild = 2 * index;
rightChild = 2 * index + 1;
if leftChild > length(tree) || isnan(tree{leftChild}.muMin)
    % There are no more children, so we work on this node
    for i=1:length(matrices)
        tree{index}.matrices{i} = tree{index}.Y' * matrices{i} * tree{index}.Y;
        tree{index}.coefficients{i} = coefficients{i};
        tree{index}.f = tree{index}.Y' * problem.rhs(:);
    end
else
    tree = computeReducedMatricesRecursive(tree, leftChild, problem, matrices, coefficients);
    tree = computeReducedMatricesRecursive(tree, rightChild, problem, matrices, coefficients);
end
end
