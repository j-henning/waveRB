% Apply the pGreedy algorithm on all leaf nodes of the tree
% Input:
% tree       - Reduced Basis model
% Np         - Maximum number of snapshots
% tolerance  - Desired tolerance
% pOne       - Problem to copy from
% splines    - Spline data
% resolution - Resolution for the evaluation
% Ouput:
% tree       - 'Fitted' RB model
function [tree] = pGreedy(tree, Np, tolerance, pOne, splines, resolution)
tree = pGreedyRecursive(tree, Np, tolerance, 1, pOne, splines, resolution);

end

% Recursive version of pGreedy
function [tree] = pGreedyRecursive(tree, Np, tolerance, index, pOne, splines, resolution)
leftChild = 2 * index;
rightChild = 2 * index + 1;
if leftChild > length(tree) || isnan(tree{leftChild}.muMin)
    % There are no more children, so we can execute the greedy algorithm
    
    maxError = zeros(Np,1);
    err = zeros(Np,1);
    for k = 1:Np
        parfor l=1:length(tree{index}.Xi)
            
            mu = tree{index}.Xi(l);
            
            % Check if we already have mu in our set. In this case, continue
            % with the next parameter mu
            if ismember(mu, tree{index}.S)
                err(l) = 0;
                continue;
            end
            
            % Create a problem for the current paramter
            p = pOne;
            p.mu = mu;
            p.A_space = mu * p.A_space;
            p.Q_space = mu^2 * p.Q_space;
            
            % Compute the reduced system
            B_N = tree{index}.Y' * (kron(p.Q_time, p.M_space) ...
                + kron(p.D_time, p.A_space') ...
                + kron(p.D_time', p.A_space) ...
                + kron(p.M_time, p.Q_space)) * tree{index}.Y;
            f_N = tree{index}.Y' * p.rhs(:);
            % Compute the reduces solution u_N_rec
            u_N = B_N \ f_N;
            u_N_rec = tree{index}.Y * u_N;
            sol_rec = getSolution(p, u_N_rec,  resolution, ...
                splines);
            
            % Compute the error
            err(l) = sqrt(mean( (tree{index}.solutions{l}-sol_rec).^2, 'all'));
        end
        % Check for which mu the error was the largest
        [maxError(k), indMax] = max(err);
        if maxError(k) < tolerance || k == Np
            tree{index}.maxError = maxError(k);
            %             maxErrors = [maxErrors, maxError];
%                         semilogy(maxError), hold on
            break;
        end
        
        tree{index}.S = [tree{index}.S, tree{index}.Xi(indMax)];
        
        
        snapshot = tree{index}.U(:,indMax);
        tree{index}.Y = [tree{index}.Y, snapshot];
        tree{index}.Y = orth(tree{index}.Y);
    end
else
    tree = pGreedyRecursive(tree, Np, tolerance, leftChild, pOne, splines, resolution);
    tree = pGreedyRecursive(tree, Np, tolerance, rightChild, pOne, splines, resolution);
end
end