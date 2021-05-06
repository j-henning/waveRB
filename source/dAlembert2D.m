% Calculates the d'Alembert solution of the wave equation on the interval
% [0,1]^2

function U = dAlembert2D (u0_x, u0_y, N, T, K)
tau = T/(K-1);

U = zeros(N,N,K);
for k=1:K
    for j=2:(N-1)
        for i=2:(N-1)
            U(i,j,k) = 0; % TODO;
        end
        
    end
end
end
