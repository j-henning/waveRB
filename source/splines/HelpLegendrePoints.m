% By C. Mollet
%
% Helpfunction to allocate Gauss-Legendre Points of order 4
% on subintervals according to refinement
%
% input:
% ------
% 'ref'    : Level of refinement
% 'l'      : Left boundary
% 'r'      : right boundary
%
% output:
% -------
% 'points' : Set of Gauss-Legendre points

function points = HelpLegendrePoints(ref,l,r)

num_Interval = 2^ref;

order = 4;
a = l:2^(-ref):r-2^(-ref);
b = l+2^(-ref):2^(-ref):r;
points = zeros(order*num_Interval,1);
for i = 1:num_Interval
    
    points(1+(i-1)*order:i*order,1) = (b(i)-a(i))/2 * ...
        [ -0.8611363115940525752239...
        -0.3399810435848562648027...
        0.3399810435848562648027...
        0.8611363115940525752239] ...
        +(a(i)+b(i))/2;
    
end