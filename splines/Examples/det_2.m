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
% SECOND form, 0-AB
offsetsol = [[1,0];[1,1]];
offsettest = [[0,1];[1,1]];

s_min =  zeros(6,7);
s_max = zeros(6,7);

parfor j = 1:6
    for k=3:7
        
        lev = j;
        jsole = [lev,lev]; % discretizaion level solution space
        jteste = [lev,lev]; % discretization level test space
        
        % Setting order
        ksol = [1,2]; % order solution space
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
        
        %% Generating
        
        % Stiff matrix
        T11 = (-1)*StiffMat(Tsole{1},Tteste{1},ksol(1),ktest(1),offsetsol(1,:),offsettest(1,:),[0,1]);
        M21 = StiffMat(Tsole{2},Tteste{2},ksol(2),ktest(2),offsetsol(2,:),offsettest(2,:),[0,0]);
        
        M12 = StiffMat(Tsole{1},Tteste{1},ksol(1),ktest(1),offsetsol(1,:),offsettest(1,:),[0,0]);
        S22 = StiffMat(Tsole{2},Tteste{2},ksol(2),ktest(2),offsetsol(2,:),offsettest(2,:),[1,1]);
        
        B = kron(T11,M21) + kron(M12,S22);
        
        % X = L_2xH^1
        RX1 = RieszMat(Tsole{1}, ksol(1) , offsetsol(1,:) , 0);
        RX2 = RieszMat(Tsole{2}, ksol(2) , offsetsol(2,:) , 1);
        % X = L_2xL_2
        %RX1 = RieszMat(Tsole{1}, ksol(1) , offsetsol(1,:) , 0);
        %RX2 = RieszMat(Tsole{2}, ksol(2) , offsetsol(2,:) , 0);
        RX = kron(RX1,RX2);
        
        % Y = L_2xH^1 \cup H^1xH^(-1)
        RY1_1 = RieszMat(Tteste{1}, ktest(1) , offsettest(1,:) , 0);
        RY2_1 = RieszMat(Tteste{2}, ktest(2) , offsettest(2,:) , 1);
        RY1_2 = RieszMat(Tteste{1}, ktest(1) , offsettest(1,:) , 1);
        RY2_2 = RieszMat(Tteste{2}, ktest(2) , offsettest(2,:) , 0) ...
            * inv(RieszMat(Tteste{2}, ktest(2) , offsettest(2,:) , 1)) ...
            * RieszMat(Tteste{2}, ktest(2) , offsettest(2,:) , 0);
        % Y = L_2xH^2 \cup H^1xL_2
        %RY1_1 = RieszMat(Tteste{1}, ktest(1) , offsettest(1,:) , 0);
        %RY2_1 = RieszMat(Tteste{2}, ktest(2) , offsettest(2,:) , 2);
        %RY1_2 = RieszMat(Tteste{1}, ktest(1) , offsettest(1,:) , 1);
        %RY2_2 = RieszMat(Tteste{2}, ktest(2) , offsettest(2,:) , 0);
        RY = kron(RY1_1,RY2_1) + kron(RY1_2,RY2_2);
        
        
        
        
        %% Calculate inf-sup and sup-sup (very inefficient!)
        
        RXi12 = full(inv(RX))^(1/2);
        RYi12 = full(inv(RY))^(1/2);
        Bt = RYi12 * B * RXi12;
        
        s = svd(Bt);
        s_max(j,k) = s(1);
        s_min(j,k) = s(numel(s));
        
    end
end


%save('s_min_jkfix12_gleich_j1-6_k3-7','s_min');
%save('s_max_jkfix12_gleich_j1-6k3-7','s_max');
