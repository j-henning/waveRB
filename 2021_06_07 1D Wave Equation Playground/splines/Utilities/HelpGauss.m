% By C. Mollet
%
% Helpfunction to calculate Gauss-Legendre quadrature of order 4
% on subintervals according to refinement.
%
% input:
% ------
% 'u'      : Function values at Gauss-Legendre points
% 'ref'    : Level of refinement
% 'l'      : Left boundary
% 'r'      : right boundary
%
% output:
% -------
% 'val' : \int_a^b f dx (resp. approx.)

function val = HelpGauss(u,ref,l,r)

num_Interval = 2^ref;
val = 0;
a = l:2^(-ref):r-2^(-ref);
b = l+2^(-ref):2^(-ref):r;

order = 4;

for i=1:num_Interval
    val = val + GaussLegendre(u(1+(i-1)*order:i*order),a(i),b(i));
end