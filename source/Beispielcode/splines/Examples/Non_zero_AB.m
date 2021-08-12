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
% First Homogenized form, non-zero AB
offsetsol = [[1,0];[1,1]];
offsettest = [[1,0];[1,1]];

for ind=1:7
    k=2;
    j=ind;
    
    
    % Setting level
    lev = j;
    jsole = [lev,lev]; % discretizaion level solution space
    jteste = [lev+1,lev+1]; % discretization level test space
    
    % Setting order
    ksol = [k,k]; % order solution space
    ktest = [k,k]; % order test space
    
    % Initialize Uniform grid
    % Test spaces
    % Implementation with cell arrays required
    Tsole = cell(1,dim);
    for i=1:dim
        Tsole(i)={InitUniNodes(jsole(i),ksol(i))};
    end
    
    % Test spaces
    % Implementation with cell arrays required
    Tteste = cell(1,dim);
    for i=1:dim
        Tteste(i)={InitUniNodes(jteste(i),ktest(i))};
    end
    
    % Set help variable of (full) dimensions:
    nsole = zeros(1,dim);
    nteste = zeros(1,dim);
    for i=1:dim
        nsole(i) = length(Tsole{i}) - ksol(i);
        nteste(i) = length(Tteste{i}) - ktest(i);
    end
    
    % Setting up matrices
    % Number of terms
    terms = 2;
    
    % % % % % % % % % % %
    
    %% Generating
    
    % Stiff matrix
    T11 = StiffMat(Tsole{1},Tteste{1},ksol(1),ktest(1),offsetsol(1,:),offsettest(1,:),[1,0]);
    M21 = StiffMat(Tsole{2},Tteste{2},ksol(2),ktest(2),offsetsol(2,:),offsettest(2,:),[0,0]);
    
    M12 = StiffMat(Tsole{1},Tteste{1},ksol(1),ktest(1),offsetsol(1,:),offsettest(1,:),[0,0]);
    S22 = StiffMat(Tsole{2},Tteste{2},ksol(2),ktest(2),offsetsol(2,:),offsettest(2,:),[1,1]);
    
    B = kron(T11,M21) + kron(M12,S22);
    
    
    % X = L_2xH^1 \cup H^1xL_2 \subset L_2xH^1 \cup H^1xH^(-1)
    % Equivalent to H^1-Norm on space-time cylinder!
    RX1_1 = RieszMat(Tsole{1}, ksol(1) , [0,0] , 0);
    RX2_1 = RieszMat(Tsole{2}, ksol(2) , [0,0] , 1);
    % Seminorm
    RX1_2 = StiffMat(Tsole{1},Tsole{1},ksol(1),ksol(1),[0,0],[0,0],[1,1]);
    RX2_2 = RieszMat(Tsole{2}, ksol(2) , [0,0] , 0);% ...
    %* inv(RieszMat(Tsole{2}, ksol(2) , [0,0] , 1)) ...
    %* RieszMat(Tsole{2}, ksol(2) , [0,0] , 0);
    RX = kron(RX1_1,RX2_1) + kron(RX1_2,RX2_2);
    %RX = kron(RX1,RX2);
    RX = sparse(RX);
    
    % Y = L_2xH^1
    RY1 = RieszMat(Tteste{1}, ktest(1) , offsettest(1,:) , 0);
    RY2 = RieszMat(Tteste{2}, ktest(2) , offsettest(2,:) , 1);
    RY = kron(RY1,RY2);
    
    
    %%%%%%%%%%%%%%%  Soluiton: (1+t)*x*(1-x) %%%%%%%%%%%%%%%%%
    
    % RHS
    % initial
    h = zeros(nteste(dim)-(offsettest(dim,1)+offsettest(dim,2)),1);
    Tt = Tteste{dim};
    
    for i=1+offsettest(dim,1):nteste(dim)-offsettest(dim,2)
        g = @(x) x.*(x-1) .* Ndiff(Tt,ktest(dim),0,i,x);
        %g = @(x) c(k) .* sin(2*pi*x).* Ndiff(Tt,ktest(dim),0,i,x);
        for step = 0:ktest(dim)-1
            h(i-offsettest(dim,1)) = h(i-offsettest(dim,1)) + Gaussq(Tt(i+step),Tt(i+step+1),g,5);
        end
    end
    f1 = h;
    
    for d = dim-1:-1:1
        h = zeros(nteste(d)-(offsettest(d,1)+offsettest(d,2)),1);
        Tt = Tteste{d};
        
        for i=1+offsettest(d,1):nteste(d)-offsettest(d,2)
            %g = @(t) -2*(1+t).*Ndiff(Tt,ktest(d),0,i,t);
            g = @(t) Ndiff(Tt,ktest(d),0,i,t);
            for step = 0:ktest(d)-1
                h(i-offsettest(d,1)) = h(i-offsettest(d,1)) + Gaussq(Tt(i+step),Tt(i+step+1),g,5);
            end
        end
        
        f1 = kron(h,f1);
    end
    
    % RHS
    % initial
    h = zeros(nteste(dim)-(offsettest(dim,1)+offsettest(dim,2)),1);
    Tt = Tteste{dim};
    
    for i=1+offsettest(dim,1):nteste(dim)-offsettest(dim,2)
        g = @(x) Ndiff(Tt,ktest(dim),0,i,x);
        %g = @(x) c(k) .* sin(2*pi*x).* Ndiff(Tt,ktest(dim),0,i,x);
        for step = 0:ktest(dim)-1
            h(i-offsettest(dim,1)) = h(i-offsettest(dim,1)) + Gaussq(Tt(i+step),Tt(i+step+1),g,5);
        end
    end
    f2 = h;
    
    for d = dim-1:-1:1
        h = zeros(nteste(d)-(offsettest(d,1)+offsettest(d,2)),1);
        Tt = Tteste{d};
        
        for i=1+offsettest(d,1):nteste(d)-offsettest(d,2)
            g = @(t) -2*(1+t).*Ndiff(Tt,ktest(d),0,i,t);
            %g = @(t) Ndiff(Tt,ktest(d),0,i,0);
            for step = 0:ktest(d)-1
                h(i-offsettest(d,1)) = h(i-offsettest(d,1)) + Gaussq(Tt(i+step),Tt(i+step+1),g,5);
            end
        end
        
        f2 = kron(h,f2);
    end
    
    % initial
    h = zeros(nteste(dim)-(offsettest(dim,1)+offsettest(dim,2)),1);
    Tt = Tteste{dim};
    
    for i=1+offsettest(dim,1):nteste(dim)-offsettest(dim,2)
        g = @(x) 2.* Ndiff(Tt,ktest(dim),0,i,x);
        for step = 0:ktest(dim)-1
            h(i-offsettest(dim,1)) = h(i-offsettest(dim,1)) + Gaussq(Tt(i+step),Tt(i+step+1),g,5);
        end
    end
    f3 = h;
    
    for d = dim-1:-1:1
        h = zeros(nteste(d)-(offsettest(d,1)+offsettest(d,2)),1);
        Tt = Tteste{d};
        
        for i=1+offsettest(d,1):nteste(d)-offsettest(d,2)
            g = @(t) Ndiff(Tt,ktest(d),0,i,t);
            %g = @(t) Ndiff(Tt,ktest(d),0,i,0);
            for step = 0:ktest(d)-1
                %h(i-offsettest(d,1)) = Ndiff(Tt,ktest(d),0,i,0);
                h(i-offsettest(d,1)) = h(i-offsettest(d,1)) + Gaussq(Tt(i+step),Tt(i+step+1),g,5);
            end
        end
        
        f3 = kron(h,f3);
    end
    
    f = f1 + f2 + f3;
    
    %% Solving
    
    U = Normal_CG(B , RY , f(:) , 5000, 1e-2 * 2^( -max(ksol) * max(jsole) ) , zeros(size(B,2),1) );
    
    
    
    %% Error estimate
    
    % Reshape
    u = reshape(U,nsole(2)-(offsetsol(2,1)+offsetsol(2,2)),nsole(1)-(offsetsol(1,1)+offsetsol(1,2)));
    % Preparing boundary
    u_full = zeros(nsole(2),nsole(1));
    u_full(1+offsetsol(2,1):nsole(2)-offsetsol(2,2) , 1+offsetsol(1,1):nsole(1)-offsetsol(1,2)) = u;
    % Since Hat functions are nodal!!
    u0 = kron(ones(nsole(1),1) , [Tsole{2}(k:length(Tsole{2})-(ksol(2)-1))].*([Tsole{2}(k:length(Tsole{2})-(ksol(2)-1))] - 1));
    
    
    val = 0;
    f = @(x)@(t) (abs((1+t)*x*(x-1) - Nev(u_full(:)+u0(:),Tsole,ksol,[t,x]))).^2;
    for t=0:2^(-lev):1-2^(-lev)
        for x=0:2^(-lev):1-2^(-lev)
            val = val + Gaussq([t,x],[t+2^(-lev),x+2^(-lev)],f,2);
        end
    end
    
    err(ind) = sqrt(val);
    
    
end;

%save('Data_Non_Zero_AB_Hom/L2error_k2','err');

