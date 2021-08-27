function [] = plotError(tree, color)
plotErrorRecursive(tree, color, 1);
end

function [] = plotErrorRecursive(tree, color, i)
% Check if we have children
if 2*i > length(tree) || isnan(tree{2*i}.muMin)
    semilogy([tree{i}.muMin, tree{i}.muMax], [tree{i}.maxError, tree{i}.maxError], 'LineWidth', 2, 'Color', color), hold on
else
    plotErrorRecursive(tree, color, 2*i);
    plotErrorRecursive(tree, color, 2*i+1);
end
end