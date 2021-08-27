function [tree] = applyHGreedy(tree, tolerance, N, pOne, resolution, splines, solver, maxIt, tolerance1, tolerance2)
[tree] = appylyHGreedyRecursive(tree, tolerance, N, pOne, resolution, splines, 1, solver, maxIt, tolerance1, tolerance2);
end

function [tree] = appylyHGreedyRecursive(tree, tolerance, N, pOne, resolution, splines, i, solver, maxIt, tolerance1, tolerance2)
% Check if we are on a leaf node
if isnan(tree{i}.muMin)
    return;
end

% Check if we already did refine this node (only happens if we call the
% applyHGreedy method on a already refined tree)

if 2*i < length(tree) && ~isnan(tree{2*i}.muMin)
    % In this case, directly process the child nodes and skip anything else
    tree = appylyHGreedyRecursive(tree, tolerance, N, pOne, resolution, splines, 2*i, solver, maxIt, tolerance1, tolerance2);
    tree = appylyHGreedyRecursive(tree, tolerance, N, pOne, resolution, splines, 2*i+1, solver, maxIt, tolerance1, tolerance2);
    return
end

[muNew, muNewIndex, tree{i}.maxError, tree{i}.Xi, tree{i}.solutions, tree{i}.U] = ...
    hGreedy(tree{i}.Xi, tree{i}.solutions, tree{i}.U, tree{i}.Y, N, tree{i}.muMin, tree{i}.muMax, pOne, resolution, splines, solver, maxIt, tolerance1, tolerance2);
% hGreedy(Xi, solutions, U, Y, N, muMin, muMax, pOne, resolution, splines, solver, tolerance1, tolerance2)

% Check if we reached the tolerance
if tree{i}.maxError < tolerance
    return;
end

tree{i}.center = 0.5 * (tree{i}.S + muNew);

% Check if we have space for children
if 2*i > length(tree)
%     fprintf('Can not refine further because there is not enough space in the tree. Error = %e\n', tree{i}.maxError)
    return;
end


% Save the left subinterval
leftIndices = find(tree{i}.Xi < tree{i}.center);
tree{2*i}.muMin = tree{i}.muMin;
tree{2*i}.muMax = tree{i}.center;
tree{2*i}.U = tree{i}.U(:,leftIndices);

tree{2*i}.Xi = tree{i}.Xi(leftIndices);
tree{2*i}.solutions = tree{i}.solutions(leftIndices);

if tree{i}.S < tree{i}.center
    tree{2*i}.Y = tree{i}.Y;
    tree{2*i}.S =  tree{i}.S;
else
    tree{2*i}.Y =  tree{i}.U(:,muNewIndex); % TODO
    tree{2*i}.S =  muNew;
end

% Save the right subinterval
rightIndices = find(tree{i}.Xi >= tree{i}.center);
tree{2*i+1}.muMin = tree{i}.center;
tree{2*i+1}.muMax = tree{i}.muMax;
tree{2*i+1}.U = tree{i}.U(:,rightIndices);
tree{2*i+1}.Xi = tree{i}.Xi(rightIndices);
tree{2*i+1}.solutions = tree{i}.solutions(rightIndices);
if tree{i}.S >= tree{i}.center
    tree{2*i+1}.Y = tree{i}.Y;
    tree{2*i+1}.S = tree{i}.S;
else
    tree{2*i+1}.Y = tree{i}.U(:,muNewIndex); % TODO
    tree{2*i+1}.S = muNew;
end

% Recursion
tree = appylyHGreedyRecursive(tree, tolerance, N, pOne, resolution, splines, 2*i, solver, maxIt, tolerance1, tolerance2);
tree = appylyHGreedyRecursive(tree, tolerance, N, pOne, resolution, splines, 2*i+1, solver, maxIt, tolerance1, tolerance2);


end