% By C. Mollet
% 1D Gauss-Legende Quadrature of Ordnung 4
%
% Input:
% ------
% 'f'    : Integrand evaluated as Gauss-Legendre points
% 'a'    : left integration limit
% 'b'    : right integration limit
%
% Output:
% -------
% 'val'  :\int_a^b f (resp. approx.)
%


function val = GaussLegendre(f,a,b)

alpha = [ ...
0.3478548451374538573731 ...
0.65214515486254614263...
0.65214515486254614263...
0.3478548451374538573731];

% Algo
val = 0;
for i=1:length(alpha)
    val = val + alpha(i)*f(i);
end
val =  (b-a)/2 * val;


end