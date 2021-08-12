% By C. Mollet
%

%% Settings
format long e;

% Since Code is in a subdirectory
addpath('../');
addpath('../Utilities');
ref = 15;


ran15 = exp( sqrt(2) * ...
    [-4.499990707309391553664 ...
    -3.669950373404452534729 ...
    -2.967166927905603248489 ...
    -2.325732486173857745454 ...
    -1.719992575186488932416 ...
    -1.136115585210920666319 ...
    -0.565069583255575748526 ...
    0 ...
    0.565069583255575748526 ...
    1.136115585210920666319 ...
    1.719992575186488932416 ...
    2.325732486173857745454 ...
    2.967166927905603248489 ...
    3.669950373404452534729 ...
    4.499990707309391553664]);

a = ran15;
c = ones(1,length(ran15));

reflev = 7;
parfor k = 1:length(ran15)
    
    %% Reference
    % Setting Dimension
    dim = 2;
    
    % space-time
    % [ [t] ; [x] ; [y] ]
    % [Left, Right]
    % Firsr formulation
    offsetsol = [[1,0];[1,1]];
    offsettest = [[1,0];[1,1]];
    
    % Setting order
    ksol = [2,2]; % order solution space
    ktest = [2,2]; % order test space
    
    % Setting level
    lev = reflev;
    jsole = [lev,lev]; % discretizaion level solution space
    jteste = [lev+1,lev+1]; % discretization level test space
    
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

    
    %% Generating
    
    % Stiff matrix
    T11 = StiffMat(Tsole{1},Tteste{1},ksol(1),ktest(1),offsetsol(1,:),offsettest(1,:),[1,0]);
    M21 = StiffMat(Tsole{2},Tteste{2},ksol(2),ktest(2),offsetsol(2,:),offsettest(2,:),[0,0]);
    
    M12 = StiffMat(Tsole{1},Tteste{1},ksol(1),ktest(1),offsetsol(1,:),offsettest(1,:),[0,0]);
    S22 = StiffMat(Tsole{2},Tteste{2},ksol(2),ktest(2),offsetsol(2,:),offsettest(2,:),[1,1]);
    %S22 = StiffMat(Tsole{2},Tteste{2},ksol(2),ktest(2),offsetsol(2,:),offsettest(2,:),[1,1]);
    
    B = kron(T11,M21) + a(k) * kron(M12,S22);
    
    % RHS
    % initial
    h = zeros(nteste(dim)-(offsettest(dim,1)+offsettest(dim,2)),1);
    Tt = Tteste{dim};
    
    for i=1+offsettest(dim,1):nteste(dim)-offsettest(dim,2)
        g = @(x) c(k) .* Ndiff(Tt,ktest(dim),0,i,x);
        %g = @(x) c(k) .* sin(2*pi*x).* Ndiff(Tt,ktest(dim),0,i,x);
        for step = 0:ktest(dim)-1
            h(i-offsettest(dim,1)) = h(i-offsettest(dim,1)) + Gaussq(Tt(i+step),Tt(i+step+1),g,5);
        end
    end
    f = h;
    
    for d = dim-1:-1:1
        h = zeros(nteste(d)-(offsettest(d,1)+offsettest(d,2)),1);
        Tt = Tteste{d};
        
        for i=1+offsettest(d,1):nteste(d)-offsettest(d,2)
            g = @(t) Ndiff(Tt,ktest(d),0,i,t);
            %g = @(t) (1+4*pi*pi*t).* Ndiff(Tt,ktest(d),0,i,t);
            for step = 0:ktest(d)-1
                h(i-offsettest(d,1)) = h(i-offsettest(d,1)) + Gaussq(Tt(i+step),Tt(i+step+1),g,5);
            end
        end
        
        f = kron(h,f);
    end
    
    %% Solving
    % Y = L_2xH^1
    RY1 = RieszMat(Tteste{1}, ktest(1) , offsettest(1,:) , 0);
    RY2 = RieszMat(Tteste{2}, ktest(2) , offsettest(2,:) , 1);
    RY = kron(RY1,RY2);
    %RX = sparse(RX);
    
    Uref = Normal_CG(B , RY , f(:) , 5000, 1e-1 * 2^( -max(ksol) * max(jsole) ) , zeros(size(B,2),1) );
 
    % Reshape
    u = reshape(Uref,nsole(2)-(offsetsol(2,1)+offsetsol(2,2)),nsole(1)-(offsetsol(1,1)+offsetsol(1,2)));
    % Preparing boundary
    u_full = zeros(nsole(2),nsole(1));
    u_full(1+offsetsol(2,1):nsole(2)-offsetsol(2,2) , 1+offsetsol(1,1):nsole(1)-offsetsol(1,2)) = u;
    exact = u_full(:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% Discretization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Setting level
    for l = 1:reflev-1;
        jsol = [l,l]; % discretizaion level solution space
        jtest = [l+1,l+1]; % discretization level test space
        
        % Initialize Uniform grid
        % Test spaces
        % Implementation with cell arrays required
        Tsol = cell(1,dim);
        for i=1:dim
            Tsol(i)={InitUniNodes(jsol(i),ksol(i))};
        end
        
        % Test spaces
        % Implementation with cell arrays required
        Ttest = cell(1,dim);
        for i=1:dim
            Ttest(i)={InitUniNodes(jtest(i),ktest(i))};
        end
        
        % Set help variable of (full) dimensions:
        nsol = zeros(1,dim);
        ntest = zeros(1,dim);
        for i=1:dim
            nsol(i) = length(Tsol{i}) - ksol(i);
            ntest(i) = length(Ttest{i}) - ktest(i);
        end

        
        %% Generating
        
        % Stiff matrix
        T11 = StiffMat(Tsol{1},Ttest{1},ksol(1),ktest(1),offsetsol(1,:),offsettest(1,:),[1,0]);
        M21 = StiffMat(Tsol{2},Ttest{2},ksol(2),ktest(2),offsetsol(2,:),offsettest(2,:),[0,0]);
        
        M12 = StiffMat(Tsol{1},Ttest{1},ksol(1),ktest(1),offsetsol(1,:),offsettest(1,:),[0,0]);
        S22 = StiffMat(Tsol{2},Ttest{2},ksol(2),ktest(2),offsetsol(2,:),offsettest(2,:),[1,1]);
        %S22 = StiffMat(Tsole{2},Tteste{2},ksol(2),ktest(2),offsetsol(2,:),offsettest(2,:),[1,1]);
        
        B = kron(T11,M21) + a(k) * kron(M12,S22);
        
        % RHS
        % initial
        h = zeros(ntest(dim)-(offsettest(dim,1)+offsettest(dim,2)),1);
        Tt = Ttest{dim};
        
        for i=1+offsettest(dim,1):ntest(dim)-offsettest(dim,2)
            g = @(x) c(k) .* Ndiff(Tt,ktest(dim),0,i,x);
            %g = @(x) c(k) .* sin(2*pi*x).* Ndiff(Tt,ktest(dim),0,i,x);
            for step = 0:ktest(dim)-1
                h(i-offsettest(dim,1)) = h(i-offsettest(dim,1)) + Gaussq(Tt(i+step),Tt(i+step+1),g,5);
            end
        end
        f = h;
        
        for d = dim-1:-1:1
            h = zeros(ntest(d)-(offsettest(d,1)+offsettest(d,2)),1);
            Tt = Ttest{d};
            
            for i=1+offsettest(d,1):ntest(d)-offsettest(d,2)
                g = @(t) Ndiff(Tt,ktest(d),0,i,t);
                %g = @(t) (1+4*pi*pi*t).* Ndiff(Tt,ktest(d),0,i,t);
                for step = 0:ktest(d)-1
                    h(i-offsettest(d,1)) = h(i-offsettest(d,1)) + Gaussq(Tt(i+step),Tt(i+step+1),g,5);
                end
            end
            
            f = kron(h,f);
        end
        
        %% Solving
        % Y = L_2xH^1
        RY1 = RieszMat(Ttest{1}, ktest(1) , offsettest(1,:) , 0);
        RY2 = RieszMat(Ttest{2}, ktest(2) , offsettest(2,:) , 1);
        RY = kron(RY1,RY2);
        %RX = sparse(RX);
        U = Normal_CG(B , RY , f(:) , 5000, 1e-1 * 2^( -max(ksol) * max(jsol) ) , zeros(size(B,2),1) );

        
        % Reshape
        u = reshape(U,nsol(2)-(offsetsol(2,1)+offsetsol(2,2)),nsol(1)-(offsetsol(1,1)+offsetsol(1,2)));
        % Preparing boundary
        u_full = zeros(nsol(2),nsol(1));
        u_full(1+offsetsol(2,1):nsol(2)-offsetsol(2,2) , 1+offsetsol(1,1):nsol(1)-offsetsol(1,2)) = u;
        uapprox = u_full(:);
        
        % Refinement onto reference grid
        utrafo = zeros(numel(exact),1);
        M1 = RefineMat(l,ksol(1));
        M2 = RefineMat(l,ksol(2));
        M = kron(M1,M2);
        s_old = size(M,1);
        utrafo(1:s_old) = M * uapprox;
        for le=l+1:reflev-1
            M1 = RefineMat(le,ksol(1));
            M2 = RefineMat(le,ksol(2));
            M = kron(M1,M2);
            s_new = size(M,1);
            utrafo(1:s_new) = M * utrafo(1:s_old);
            s_old = s_new;
        end
        
        diffsol = utrafo - exact;
        
        % X = L_2xH^1 \cup H^1xL_2 \subset L_2xH^1 \cup H^1xH^(-1)
        % Equivalent to H^1-Norm on space-time cylinder!
        RX1_1 = RieszMat(Tsole{1}, ksol(1) , [0,0] , 0);
        RX2_1 = RieszMat(Tsole{2}, ksol(2) , [0,0] , 1);
        % Seminorm
        RX1_2 = StiffMat(Tsole{1},Tsole{1},ksol(1),ksol(1),[0,0],[0,0],[1,1]);
        RX2_2 = RieszMat(Tsole{2}, ksol(2) , [0,0] , 0);% ...
        %* inv(RieszMat(Tsole{2}, ksol(2) , [0,0] , 1)) ...
        %* RieszMat(Tsole{2}, ksol(2) , [0,0] , 0);
        
        %RX1 = StiffMat(Tsole{1},Tsole{1},ksol(1),ksol(1),[0,0],[0,0],[1,1]);
        %RX2 = StiffMat(Tsole{2},Tsole{2},ksol(2),ksol(2),[0,0],[0,0],[1,1]);
        %RX1 = RieszMat(Tsole{1},ksol(1),[0,0],1);
        %RX2 = RieszMat(Tsole{2},ksol(2),[0,0],1);
        
        RX = kron(RX1_1,RX2_1) + kron(RX1_2,RX2_2);
        %RX = kron(RX1,RX2);
        RX = sparse(RX);
        
        %e(k,lev) = sqrt(transpose(diffsol) * RX * (diffsol));
        e(k,l) = sqrt(transpose(diffsol) * RX * (diffsol));
        
        
    end

end % end ran loop


for l=1:reflev-1
    L1(l) = GaussHermite(e(:,l),ref);
    L2(l) = GaussHermite(e(:,l).^2,ref)^(1/2);
    L3(l) = GaussHermite(e(:,l).^3,ref)^(1/3);
    L20(l) = GaussHermite(e(:,l).^20,ref)^(1/20);
end

%save('Data_LogNormal/refInt15_RefLev7_Gauss15','e');
%save('Data_LogNormal/refInt15_RefLev7_Gauss15_L1','L1');
%save('Data_LogNormal/refInt15_RefLev7_Gauss15_L2','L2');
%save('Data_LogNormal/refInt15_RefLev7_Gauss15_L3','L3');
%save('Data_LogNormal/refInt15_RefLev7_Gauss15_L20','L20');


