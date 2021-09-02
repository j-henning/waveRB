% Visualize the hierarchical interval splitting of a Reduced Basis model
% Input:
% tree - Reduced basis model
function [] = plotIntervals(tree)
colors = {'b', 'r', 'g', 'm', 'y'};
hold on
height = floor( log(length(tree)) / log(2))+1;
for i=1:length(tree)
    if isnan(tree{i}.muMin)
        continue;
    end
    
    level = floor( log(i) / log(2));
    plot([tree{i}.muMin, tree{i}.muMax], [level, level], ...
        strcat(colors{mod((level), length(colors))+1 }, '+-'), 'LineWidth', 3)
    plot(tree{i}.Xi, level*ones(length(tree{i}.Xi),1), 'k|', 'LineWidth', 2)
    
    if level > 0
        plot([tree{i}.muMin, tree{i}.muMin], [level, level-1], 'k-', 'LineWidth', 2)
        plot([tree{i}.muMax, tree{i}.muMax], [level, level-1], 'k-', 'LineWidth', 2)
    end
    
    if level < height - 1 && ~isnan(tree{2*i}.muMin)
        for j = 1:length(tree{i}.Xi)
            plot([tree{i}.Xi(j), tree{i}.Xi(j)], [level, level+1], 'k:', 'LineWidth', 1)
        end
    end
end
xlabel('\mu', 'FontSize', 20)
ylabel('Refinement level', 'FontSize',13)
curtick = get(gca, 'yTick');
yticks(unique(round(curtick)));
drawnow
end