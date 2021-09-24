% Displays the space-time solution in an easy manner, depending of the
% dimension:
% 1D: Surf plot
% 2D: TBD
% 3D: Animation of a surf plot for which the third dimension gets only
%     shown in the center
% Input:
% solution   - Solution vector
% resolution - Resolution used
% export     - If a video should be exported (default = false)
% name       - Name of the video file (default = 'wave.avi')
function [] = plotSolution(solution, resolution, export, name)

if nargin < 3
   export = false;
end

if nargin < 4
    name = 'wave.avi';
end


if export && ismatrix(solution)
    warning('Can not create video for 1D example')
    export = false;
elseif export
    v = VideoWriter(name);
    v.FrameRate = ceil(size(solution,ndims(solution))/10);
    open(v);
end

switch ndims(solution)
    case 2 % 1D case, we can just use a simple surf plot
        x = linspace(0, 1, 2^resolution.x + 1);
        t = linspace(0, 1, 2^resolution.t + 1);
        
        [X, T] = meshgrid(x,t);
        s = surf(X, T, solution); 
        s.EdgeColor = 'none';
        xlabel('t')
        ylabel('x')
        zlabel('u(x,t)')
    case 3
        minVal = min(solution, [], 'all');
        maxVal = max(solution, [], 'all');
        f = figure('WindowState','maximized');
        for i=1:size(solution,3)
            s = surf(solution(:,:, i)); 
            s.EdgeColor = 'none';
            title(num2str(i/size(solution,4)));
            zlim([minVal, maxVal])
            drawnow
            if export
               frame = getframe(gcf);
               writeVideo(v,frame); 
            end
%              pause(0.25)
        end
    case 4 % 3D case: We slice through the solution at various times
        minVal = min(solution, [], 'all');
        maxVal = max(solution, [], 'all');
        f = figure('WindowState','maximized');
        
        % The whole animation should took around 10 seconds to display (if
        % possible)
        durating = 10 / size(solution,4); 
        
        for i=1:size(solution,4)
            tic;
            s = surf(squeeze(solution(:,:,floor(0.5 * 2^resolution.z), i))); s.EdgeColor = 'none';
            title(num2str(i/size(solution,4)));
             zlim([minVal, maxVal])
            drawnow
            t = toc;
            if t < durating
               pause(durating - t); 
            end
            if export
               frame = getframe(gcf);
               writeVideo(v,frame); 
            end
%             pause(0.25)
        end
        
end

if export
    close(v);
end
end