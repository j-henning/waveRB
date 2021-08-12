% By C. Mollet
%

%% Settings
format long e;

% Since Code is in a subdirectory
addpath('../');
addpath('../Utilities');

% Setting Dimension
dim = 2;

% space-time
% [ [t] ; [x] ; [y] ]
% [Left, Right]
% First Form, Zero AB
offset = [[1,0];[1,1]];

% Order
k = 7;
K = [k,k];

% Level
for lev = 1:7;
    J = [lev,lev];
    
    % Initialize Uniform grid
    T = cell(1,dim);
    for i=1:dim
        T(i)={InitUniNodes(J(i),K(i))};
    end
    
    % Set help variable of (full) dimensions:
    n = zeros(1,dim);
    for i=1:dim
        n(i) = length(T{i}) - K(i);
    end
    
    % Set Collocation Points (Greville)
    %x = zeros(n-offset(1)-offset(2),1);
    x = cell(1,dim);
    for ind = 1:dim
        for i=1+offset(ind,1):n-offset(ind,2)
            x{ind}(i-offset(1)) = sum(T{ind}(i+1:i+K(ind)-1))/(K(ind)-1);
        end
    end
    
    % Assemble Collocation matrices
    T11 = CollMat(T{1},K(1),offset(1,:),x{1},1);
    M21 = CollMat(T{2},K(2),offset(2,:),x{2},0);
    
    M12 = CollMat(T{1},K(1),offset(1,:),x{1},0);
    S22 = (-1)*CollMat(T{2},K(2),offset(2,:),x{2},2);
    
    B = kron(T11,M21) + kron(M12,S22);
    %B = (-1)*CollMat(T{1},K(1),offset(1,:),x{1},2);
    
    
    % RHS
    gt = @(t) (1+4*pi*pi*t);
    gx = @(x) sin(2*pi*x);
    %gt = @(t) t.^0;
    %gx = @(x) x.^0;
    
    ft = gt(x{1}');
    fx = gx(x{2}');
    f = kron(ft,fx);
    
    % Solve
    U = B\f;
    
    
    % Reshaping
    u = reshape(U,n(2)-(offset(2,1)+offset(2,2)),n(1)-(offset(1,1)+offset(1,2)));
    % Preparing boundary
    u_full = zeros(n(2),n(1));
    u_full(1+offset(2,1):n(2)-offset(2,2) , 1+offset(1,1):n(1)-offset(1,2)) = u;
    
    % Calculate error
    val = 0;
    f = @(x)@(t) (abs(sin(2*pi*x)*t - Nev(u_full(:),T,K,[t,x]))).^2;
    t=0:2^(-lev):1-2^(-lev);
    %for t=0:2^(-lev):1-2^(-lev)
    parfor i=1:length(t)
        for x=0:2^(-lev):1-2^(-lev)
            %y=0:2^(-lev):1-2^(-lev);
            val = val + Gaussq([t(i),x],[t(i)+2^(-lev),x+2^(-lev)],f,5);
        end
    end
    
    err(lev) = sqrt(val);
end

%save('Data_Col/L2error_k7_1_7','err');

