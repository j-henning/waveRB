function [err] = calculate2DL2Error(p, sol)
if p.has_analytical_solution
    % TODO: Check meshgrid / ndgrid!
    [X, Y, T] = meshgrid(linspace(0,1, size(sol,1)), ...
        linspace(0,1, size(sol,2)), linspace(0,1, size(sol,3)));
    sol_ana = p.u_analytical(X,Y,T);
else % Load the desired solution
    load(p.referenceSolutionPath)
end
err = sqrt(mean( (sol-sol_ana).^2, 'all'));
end
