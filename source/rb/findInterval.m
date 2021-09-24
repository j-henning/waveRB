% Returns the index of the interval which contains mu
% Input:
% tree  - Reduced basis model
% mu    - parameter to search for
% Output
% index - Index of the interval containing mu (-1 is mu is lies not inside
%         the tree)
function [index] = findInterval(tree, mu)
if mu < tree{1}.muMin || mu > tree{1}.muMax
    index =  -1;
else
    index = findIntervalRecursive(tree, 1, mu);
end
end

% Recursive version of findInterval
function [index] = findIntervalRecursive(tree, index, mu)
leftChild = 2 * index;
rightChild = 2 * index + 1;
if leftChild > length(tree)
    % There are no more children, so we are done
    return;
elseif isnan(tree{leftChild}.muMin)
    % The interval is no further refinded
    return;
elseif tree{index}.center > mu
    index = findIntervalRecursive(tree, leftChild, mu);
else
    index = findIntervalRecursive(tree, rightChild, mu);
end
end