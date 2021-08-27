function [solutionTS] = get1DTSSolution(problemTS, UTS, resolution, splines)
% Add zeros coefficients for the boundaries
u_full = zeros(size(UTS,1) + 2, size(UTS,2));
u_full(2:end-1,:) = UTS;


solutionTS = zeros(2^resolution.x + 1, problemTS.K);

for i=1:size(u_full,1) % every space node in x
    spline_x = splines.space(i,:);
    
    
    space_x_indices = find(spline_x, 1, 'first'):find(spline_x, 1, 'last');
    
    
    
    
    for k=1:problemTS.K % Every time step
        if u_full(i,k) < eps && u_full(i,k) > - eps
            continue;
        end
        
        spline = spline_x(space_x_indices);
        solutionTS(space_x_indices, k) = solutionTS(space_x_indices, k) +  u_full(i,k) * spline';
    end
end

% Linearaly interpolate the solution to match Kref time steps

x = linspace(0,1,2^resolution.x + 1);
[X, Y] = meshgrid(x, linspace(0,1, problemTS.K));
[Xq, Yq] = meshgrid(x, linspace(0,1, 2^resolution.t + 1));
solutionTS = interp2(X,Y,solutionTS',Xq,Yq)';



end