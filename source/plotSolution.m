function [] = plotSolution(solution, resolution)

switch ndims(solution)
    case 2 % 1D case, we can just use a simple surf plot
        x = linspace(0, 1, 2^resolution.x + 1);
        t = linspace(0, 1, 2^resolution.t + 1);
        
        [X, T] = meshgrid(x,t);
        surf(X, T, solution)
        xlabel('t')
        ylabel('x')
        zlabel('u(x,t)')
    case 3
        error('Not yet implemented')
    case 4 % 3D case: We slice through the solution at various times
        f = figure('WindowState','maximized');
        for i=1:size(solution,4)
            s = surf(squeeze(solution(:,:,floor(0.5 * 2^resolution.z), i))); s.EdgeColor = 'none';
            title(num2str(i/size(solution,4))); 
            pause(0.1)
            drawnow   
        end
        
end

end