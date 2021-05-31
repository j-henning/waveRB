clear
close all
clc

% u0 = @(r) sin(2*pi*r);
u0 = @(r) 1.* (r < 0.25);
% u0 = @(r) r;

c = .2;

t = linspace(0,1,100);
r = linspace(0,1);

% for k=1:length(t)
%     for j=1:length(r)
%         sol(j) =  dAlemenbertSphere(r(j),t(k),c, u0);
%     end
%     plot(r,sol)
%     drawnow
%     title(num2str(t(k)))
%     pause(0.5)
%     
% end

[X, Y] = meshgrid(linspace(0,1,100), linspace(0,1,100));
% 
for i=1:100
    for j=1:100
        r = sqrt((0.5 -X(i,j))^2 + (0.5 - Y(i,j))^2);
        sol(i,j,:) = dAlemenbertSphere(r,t,c, u0);
    end
end


zlimits = [min(sol,[], 'all'), max(sol,[],'all')];
for k=1:length(t)
    s = surf(sol(:,:,k)); s.LineStyle = 'none';
    title(num2str(t(k)))
 %  zlim(zlimits)
    drawnow
    
    pause(0.1)
end
        


function sol = dAlemenbertSphere(r,t,c, u0)
if r > c * t
    sol = 1./(2.*r) .*( (r+c.*t) .* u0(r + c.*t) + (r - c.*t).*u0(r-c.*t));
else
    sol = 1./(2.*r) .*( (r+c.*t) .* u0(r + c.*t) - (c.*t - r).*u0(r-c.*t));
end
end