function [solutionTS] = get3DTSSolution(problemTS, UTS, resolution, splines)
% Add zeros coefficients for the boundaries

space_size = nthroot(size(UTS,1),3) + sum(problemTS.offset_space);
u_full = zeros(space_size,space_size,space_size,size(UTS,2));
u_full(1+problemTS.offset_space(1):end-problemTS.offset_space(2), ...
    1+problemTS.offset_space(1):end-problemTS.offset_space(2), ...
    1+problemTS.offset_space(1):end-problemTS.offset_space(2),:) = ...
    reshape(UTS, [nthroot(size(UTS,1),3),nthroot(size(UTS,1),3), ...
    nthroot(size(UTS,1),3), size(UTS,2)]);


solutionTS = zeros(2^resolution.x + 1, ...
    2^resolution.y + 1, ...
    2^resolution.z + 1, problemTS.K);

for i=1:size(u_full,1) % every space node in x
    spline_x = splines.space(i,:);
    
    
    space_x_indices = find(spline_x, 1, 'first'):find(spline_x, 1, 'last');
    
    
    for j=1:size(u_full,2) % every space node in y
        spline_y = splines.space(j,:);
        space_y_indices = find(spline_y, 1, 'first'):find(spline_y, 1, 'last');
        
        for k=1:size(u_full,3) % every space node in z
            spline_z = splines.space(k,:);
            space_z_indices = find(spline_z, 1, 'first'):find(spline_z, 1, 'last');
            
            for l=1:problemTS.K % Every time step
                if u_full(i,j,k,l) < eps && u_full(i,j,k,l) > - eps
                    continue;
                end

                spline = kron(spline_z(space_z_indices), kron(spline_y(space_y_indices), spline_x(space_x_indices)));
                spline = reshape(spline, ...
                    [length(space_x_indices), length(space_y_indices),...
                    length(space_z_indices)]);
                
                solutionTS(space_x_indices, space_y_indices, space_z_indices, l) = ...
                    solutionTS(space_x_indices, space_y_indices, space_z_indices, l) + ...
                    u_full(i,j,k,l) * spline;
            end
        end
    end
end



% Linearaly interpolate the solution to match Kref time steps

% x = linspace(0,1, 2^resolution.x+1);
% y = linspace(0,1, 2^resolution.y+1);
% z = linspace(0,1, 2^resolution.z+1);
% [X, Y, Z, T] = ndgrid(x, y, z, linspace(0,1, problemTS.K));
% [Xq, Yq, Zq, Tq] = ndgrid(x, y, z, linspace(0,1,2^resolution.t+1));
% solutionTS = interpn(X, Y, Z, T, solutionTS, Xq, Yq, Zq, Tq);

t_problem = linspace(0,1, problemTS.K);
t_evaluation = linspace(0,1,2^resolution.t+1);

solutionTSFull = solutionTS;
solutionTS = zeros(2^resolution.x+1, 2^resolution.y+1, 2^resolution.z+1, ...
    2^resolution.t+1);

for i=1:2^resolution.t+1
    [~, index_closest] = min(abs(t_problem-t_evaluation(i)));
    solutionTS(:,:,:,i) = solutionTSFull(:,:,:,index_closest);

end





end