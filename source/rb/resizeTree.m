% Changes the size of a tree by either throwing some data away or by adding
% new empty data
% Input:
% tree      - Old tree
% newHeight - New desired height
% Output
% tree      - New tree with the old data and the height newHeight
function [tree] = resizeTree(tree, newHeight)
currentHeight = floor( log(length(tree)) / log(2))+1;

if currentHeight > newHeight
    tree = tree(1:newHeight);
elseif currentHeight < newHeight
    i = 2^newHeight - 1;
    while i > 2^currentHeight - 1
        tree{i} = struct('muMin', NaN, 'muMax', NaN, 'U', NaN, 'Y', NaN, ...
            'Xi', NaN,'S', NaN, 'solutions', NaN, 'maxError', Inf, 'center', NaN);
        i = i - 1;
    end
end
end