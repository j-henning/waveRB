% Returns a vector containing all basis sizes
% Input:
% tree  - Reduced basis moded
% Output
% sizes - [Vector containing all basis sizes (unsorted)
function [sizes] = returnBasisSize(tree)
sizes = returnBasisSizeRecursive(tree, 1, []);
end

% Recursive version of findInterval
function [sizes] = returnBasisSizeRecursive(tree, index, sizes)
leftChild = 2 * index;
rightChild = 2 * index + 1;
if leftChild > length(tree)
    % There are no more children, so we are done
    sizes = [sizes, size(tree{index}.Y,2)];
elseif isnan(tree{leftChild}.muMin)
    % The interval is no further refinded
    sizes = [sizes, size(tree{index}.Y,2)];
else
    sizes = returnBasisSizeRecursive(tree, leftChild, sizes);
    sizes = returnBasisSizeRecursive(tree, rightChild, sizes);
end
end