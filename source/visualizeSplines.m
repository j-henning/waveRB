clear
close all
clc

addpath('../splines');
addpath('../splines/Utilities');
addpath('../solvers');


% x = linspace(0,1,500);
% refinement = 3;
% 
% for index=1:2^refinement
% f = timeSplines(refinement,index,1);
% 
% plot(x, f(x)), hold on
% end
% 


% T = linspace(0,1,9);
% T = [-T(2), T, T(end) + T(2)]
% x = linspace(0,1,500);
% 
% 
% for i=1:8
% f = @(x) -(x - T(i)).^3 .* (x - T(i+2)).^3 .* (x >= T(i)) .* (x <= T(i+2));
% fxx = @(x) 6 * (T(i) -x) .* (x - T(i+2)) .* ...
%     (T(i)^2 + 3 * T(i) * T(i+2) - 5 * T(i) * x + T(i+2)^2 - 5 * T(i+2) * x ...
%     + 5 * x.^2).* (x >= T(i)) .* (x <= T(i+2));
% subplot(1,2,1)
% plot(x, f(x)), hold on
% subplot(1,2,2)
% plot(x, fxx(x)), hold on
% 
% end



T = InitUniNodes(3, 3);

splines = [];

% T = linspace(-1,1,9);
x = linspace(0,1,32);

% for i=2:length(T)-4
%     
%     val = [];
%     val2 = [];
%     for j=1:length(x)
%         val(j) = Ndiff(T,3,0,i,x(j));
%         val2(j) = Ndiff(T,3,1,i,x(j));
%     end
%     
%     splines = [splines; val];
%     plot(x, val), hold on
% %     plot(val2, '--')
% end
% 
% splines = [x; splines]'

i = 5;

val = [];
val2 = [];
for j=1:length(x)
        val(j) = Ndiff(T,3,0,i,x(j));
        val2(j) = Ndiff(T,3,2,i,x(j));
end

spline = kron(val2, val') - kron(val, val2');
surf(spline)


[X, Y] = meshgrid(x,x);

data = [X(:), Y(:), spline(:)];
