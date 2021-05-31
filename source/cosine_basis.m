% for i = 1:10
% 
% plot(cos((i*2+1)/2 * pi * linspace(0,1))), hold on
% drawnow
% pause(2)
% end

num = 32;
Q = zeros(num,num);
M = zeros(num,num);
N = zeros(num,num);
for i = 1:num
%        for j = 1:num
j = i;
       funM = @(t) cos((i*2+1)/2 * pi * t) .* cos((j*2+1)/2 * pi * t);
       funQ = @(t) cos((i*2+1)/2 * pi * t) .* ((i*2+1)/2 * pi)^2 ...
           .* cos((j*2+1)/2 * pi * t) .* ((j*2+1)/2 * pi)^2;
       funN = @(t) -cos((i*2+1)/2 * pi * t) .* ((i*2+1)/2 * pi)^2 ...
           .* cos((j*2+1)/2 * pi * t);
        Q(i,j) = integral(funQ, 0, 1);
        M(i,j) = integral(funM, 0, 1);
        N(i,j) = integral(funN, 0, 1);
%        end
end



cond(Q)
cond(M)
cond(N)