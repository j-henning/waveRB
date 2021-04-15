function [err] = calculate1DL2Error(p, sol)
if p.has_analytical_solution
    [X, T] = ndgrid(linspace(0,1, size(sol,1)), linspace(0,1, size(sol,2)));
    sol_ana = p.u_analytical(X,T);
elseif p.dAlemebert
    error('Not yet implemented!')
else % Load the desired solution
    error('Not yet implemented!')
end
err = sqrt(mean( (sol-sol_ana).^2, 'all'));
end

