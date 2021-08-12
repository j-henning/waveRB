function [] = plotSolution(solution, resolution)

switch ndims(solution)
    case 2
        x = linspace(0, 1, 2^resolution.x + 1);
        t = linspace(0, 1, 2^resolution.t + 1);
        
        [X, T] = meshgrid(x,t);
        surf(X, T, solution)
        xlabel('t')
        ylabel('x')
        zlabel('u(x,t)')
    case 3
        error('Not yet implemented')
    case 4
        error('Not yet implemented')
        
end

end