function [tree] = enrich(tree, Np, pOne, splines, resolution, solver, maxIt, tolerance1, tolerance2)
tree = enrichRecursive(tree, Np, 1, pOne, splines, resolution, solver, maxIt, tolerance1, tolerance2);

end

% Recursively search
function [tree] = enrichRecursive(tree, Np, index, pOne, splines, resolution, solver, maxIt, tolerance1, tolerance2)
leftChild = 2 * index;
rightChild = 2 * index + 1;
if leftChild > length(tree) || isnan(tree{leftChild}.muMin)
    % There are no more children, so we refine this node
    newN = Np  - length(tree{index}.Xi);
    XiNew = rand(newN, 1) * (tree{index}.muMax - tree{index}.muMin) ...
        + tree{index}.muMin;
    snapshotsNew = cell(newN,1);
    UNew = zeros(size(tree{index}.U,1), newN);
    
    % Compute the new snapshots
    parfor k = 1:newN
        mu = XiNew(k);
        problem = changeWaveSpeed(pOne, mu);

        UNew(:,k) = solveProblem(problem, solver, maxIt, tolerance1, tolerance2);
        snapshotsNew{k} = getSolution(problem, UNew(:,k), ...
            resolution, splines);
    end
    % Concatinate the parameters and the snapshots
    tree{index}.Xi = [tree{index}.Xi; XiNew];
    tree{index}.solutions = [tree{index}.solutions; snapshotsNew];
    tree{index}.U = horzcat(tree{index}.U, UNew);
else
    tree = enrichRecursive(tree, Np, leftChild, pOne, splines, resolution, solver, maxIt, tolerance1, tolerance2);
    tree = enrichRecursive(tree, Np, rightChild, pOne, splines, resolution, solver, maxIt, tolerance1, tolerance2);
end
end