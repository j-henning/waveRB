function [index] = findInterval(tree, mu)
    index = findIntervalRecursive(tree, 1, mu);
end

% Recursively search
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