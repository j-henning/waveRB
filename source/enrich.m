function [tree] = enrich(tree, Np, pOne, splines, resolution)
tree = enrichRecursive(tree, Np, 1, pOne, splines, resolution);

end

% Recursively search
function [tree] = enrichRecursive(tree, Np, index, pOne, splines, resolution)
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
        problem = pOne;
        
        problem.mu = mu;
        problem.A_space = mu * problem.A_space;
        problem.Q_space = mu^2 * problem.Q_space;
        
        
        UNew(:,k) = solveProblem(problem);
        snapshotsNew{k} = get1DsolutionSmart(problem, UNew(:,k), ...
            resolution, splines.splines_time, ...
            splines.splines_2_time, ...
            splines.splines_space, ...
            splines.splines_2_space);
    end
    % Concatinate the parameters and the snapshots
    tree{index}.Xi = [tree{index}.Xi; XiNew];
    tree{index}.solutions = [tree{index}.solutions; snapshotsNew];
    tree{index}.U = horzcat(tree{index}.U, UNew);
else
    tree = enrichRecursive(tree, Np, leftChild, pOne, splines, resolution);
    tree = enrichRecursive(tree, Np, rightChild, pOne, splines, resolution);
end
end