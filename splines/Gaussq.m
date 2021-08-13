% By C. Mollet
%
% Gauss-Legende Quadrature of Ordnung k (exactness up to degree 2k-1)
% Works for arbitrary dimensions given via the integration limits
%
% Input:
% ------
% 'a'    : left integration limits (array)
% 'b'    : right integration limits (array)
% 'f'    : Integrand
% 'k'    : Order
%
% Output:
% -------
% 'val'  :\int_a^b f (resp. approx.)
%
% NOTE: For technical reasons, f needs to be in a nested format, i.e.,
% instead of f(x,y) use f(x)(y). In Matlab: Instead of @(x,y) f, use
% @(y)@(x) f.
% There is a workaround by Matlab but since it is not that important, I
% haven't implemented it yet

function val = Gaussq(a,b,f,k)

% Case to capture, resp. give a meaningful error, "unnested format"
if nargin(f) ~= 1
    error('The integrand f needs to be in a nested format, i.e., instead of e.g. @(x,y) f, use @(x)@(y) f.');
end

dim = length(a);

% switch over order 'k'
switch k
    
    case 1
        % setting weights und nodes
        alpha = 2;
        x = 0;
        
        if dim == 1
            % Algo
            val = 0;
            val = val + alpha*f((b(1)-a(1))/2*x + (a(1)+b(1))/2);
            val = val*(b(1)-a(1))/2;
            
        else
            val = 0;
            val = val + alpha*Gaussq(a(1:dim-1),b(1:dim-1),f((b(dim)-a(dim))/2*x + (a(dim)+b(dim))/2),k);
            val = val*(b(dim)-a(dim))/2;
        end
        
        
    case 2
        
        % setting weights und nodes
        alpha = [1 1];
        x = [-5.773502691896258e-01 5.773502691896258e-01];
        
        if dim == 1
            % Algo
            %val = 0;
            h = (b(1)-a(1))/2;
            m = (a(1)+b(1))/2;
            val = h*sum( alpha.*arrayfun(f , h*x+m ) );
            
        else
            val = 0;
            for i=1:2
                val = val + alpha(i)*Gaussq(a(1:dim-1),b(1:dim-1),f((b(dim)-a(dim))/2*x(i) + (a(dim)+b(dim))/2),k);
            end
            val = val*(b(dim)-a(dim))/2;
        end
        
    case 3
        
        % setting weights und nodes
        alpha = [5.555555555555556e-01 , 8.888888888888889e-01 , 5.555555555555556e-01];
        x = [-7.745966692414834e-01 , 0 , 7.745966692414834e-01];
        
        if dim == 1
            % Algo
            val = 0;
            for i=1:3
                val = val + alpha(i)*f((b(1)-a(1))/2*x(i) + (a(1)+b(1))/2);
            end
            val = val*(b(1)-a(1))/2;
            
        else
            val = 0;
            for i=1:3
                val = val + alpha(i)*Gaussq(a(1:dim-1),b(1:dim-1),f((b(dim)-a(dim))/2*x(i) + (a(dim)+b(dim))/2),k);
            end
            val = val*(b(dim)-a(dim))/2;
        end
        
    case 4
        % setting weights und nodes
        alpha = [3.478548451374539e-01 , 6.521451548625462e-01 , 6.521451548625462e-01 , 3.478548451374539e-01 ];
        x = [ -8.611363115940526e-01 , -3.399810435848563e-01 , 3.399810435848563e-01 , 8.611363115940526e-01 ];
        
        if dim == 1
            % Algo
            val = 0;
            for i=1:4
                val = val + alpha(i)*f((b(1)-a(1))/2*x(i) + (a(1)+b(1))/2);
            end
            val = val*(b(1)-a(1))/2;
            
        else
            val = 0;
            for i=1:4
                val = val + alpha(i)*Gaussq(a(1:dim-1),b(1:dim-1),f((b(dim)-a(dim))/2*x(i) + (a(dim)+b(dim))/2),k);
            end
            val = val*(b(dim)-a(dim))/2;
        end
        
        
    otherwise
        % setting weights und nodes
        alpha = [2.369268850561891e-01 , 4.786286704993665e-01 , 5.688888888888889e-01 , 4.786286704993665e-01 , 2.369268850561891e-01];
        x = [-9.061798459386640e-01 , -5.384693101056830e-01 , 0 , 5.384693101056830e-01 , 9.061798459386640e-01];
        
        if dim == 1
            % Algo
            val = 0;
            for i=1:5
                val = val + alpha(i)*f((b(1)-a(1))/2*x(i) + (a(1)+b(1))/2);
            end
            val = val*(b(1)-a(1))/2;
            
        else
            val = 0;
            for i=1:5
                val = val + alpha(i)*Gaussq(a(1:dim-1),b(1:dim-1),f((b(dim)-a(dim))/2*x(i) + (a(dim)+b(dim))/2),k);
            end
            val = val*(b(dim)-a(dim))/2;
        end
        
end

end