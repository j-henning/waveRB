% Calculates the d'Alembert solution of the wave equation on the interval
% [0,1]
% Input:
%       u0 - function for the initial position
%       V0 - integrated function for the initial velocity v0
%       c  - wave speed
%       N  - number of points
%       T  - length of the time interval [0,T]
%       K  - number of time steps
%
% Output:
%       U  - solution matrix (N x K)
function U = dAlembert1D (u0, V0, c, N, T, K)
tau = T/(K-1);

U = zeros(N,K);
for k=1:K
    for j=2:(N-1)
        xLeft = (j-1)/(N-1)-c*tau*(k-1);
        xRight = (j-1)/(N-1)+c*tau*(k-1);
        signumLeft = 1;
        signumRight = 1;
        
        if xLeft < 0
            if mod(ceil(xLeft),2) == 0
                signumLeft = -1;
                xLeft = abs(xLeft) - floor(abs(xLeft));
            else
                xLeft = 1 - (abs(xLeft) - floor(abs(xLeft)));
            end
            
        end
        
        if xRight > 1
            if mod(floor(xRight),2) ~= 0
                signumRight = -1;
                xRight = 1 -(xRight - floor(xRight));
            else
                xRight = xRight - floor(xRight);
            end
        end
        
        U(j,k) = 0.5 * (signumLeft*u0(xLeft) + signumRight *u0(xRight)) ...
            + (V0(xRight)   - V0(xLeft) )/(2*c);
    end
end
end
