% By C. Mollet
%

%% Settings
format long e;

% Since Code is in a subdirectory
addpath('../');
addpath('../Utilities');

% Setting Dimension
dim = 2;

% Setting boundary offset
% 2D space-time
% [ [t] ; [x] ]
% [Left, Right]
offsetsol = [[0,0];[1,1]];
offsettest = [[0,1];[1,1]];


%% Grid and Order
% init
ksol  = zeros(1,dim); % order solution space
ktest = zeros(1,dim); % order test space
jsol  = zeros(1,dim); % discretizaion level solution space
jtest = zeros(1,dim); % discretization level test space

e = zeros(1,8);
quot = zeros(1,7);

for lev = 1:8
    % setting
    ksol = [2,2];
    ktest = [2,2];
    
    jsol = [lev,lev];
    jtest = [lev+1,lev+1];
    
    %% Setting up matrices and solving
    % Here: space-time weak Schroedinger equation in first form
    % Initialize set of uniform nodes corresponding to B-Spline basis
    % with boundary adaptation offset:
    
    % Solution spaces
    Tsol1 = InitUniNodes(jsol(1),ksol(1));
    Tsol2 = InitUniNodes(jsol(2),ksol(2));
    %Tsol3 = InitUniNodes(jsol(3),ksol(3));
    
    % Test spaces
    Ttest1 = InitUniNodes(jtest(1),ktest(1));
    Ttest2 = InitUniNodes(jtest(2),ktest(2));
    %Ttest3 = InitUniNodes(jtest(3),ktest(3));
    
    % (full) dimensions:
    nsol(1) = length(Tsol1) - ksol(1);
    nsol(2) = length(Tsol2) - ksol(2);
    %nsol(3) = length(Tsol3) - ksol(3);
    ntest(1) = length(Ttest1) - ktest(1);
    ntest(2) = length(Ttest2) - ktest(2);
    %ntest(3) = length(Ttest3) - ktest(3);
    
    
    
    % Generate system matrices
    % Mass Matrix in t:
    Mt = StiffMat(Tsol1,Ttest1,ksol(1),ktest(1),offsetsol(1,:),offsettest(1,:),[0,0]);
    % 'Time Derivative Matrix':
    % first form
    % Dt = StiffMat(Tsol1,Ttest1,ksol(1),ktest(1),offsetsol(1,:),offsettest(1,:),[1,0]);
    % second form:
    Dt = (-1) * StiffMat(Tsol1,Ttest1,ksol(1),ktest(1),offsetsol(1,:),offsettest(1,:),[0,1]);
    
    % Mass Matrix in x:
    Mx = StiffMat(Tsol2,Ttest2,ksol(2),ktest(2),offsetsol(2,:),offsettest(2,:),[0,0]);
    % Stiffnes Matrix in x:
    Sx = StiffMat(Tsol2,Ttest2,ksol(2),ktest(2),offsetsol(2,:),offsettest(2,:),[1,1]);
    
    % Tensorization:
    B = kron(Dt,Mx) + 1i*kron(Mt,Sx);
    
    % Assamble right hand side
    % No tensor product
    f = sparse(ntest(2)-(offsettest(2,1)+offsettest(2,2)),ntest(1)-(offsettest(1,1)+offsettest(1,2)));
    for jt=1+offsettest(1,1):ntest(1)-offsettest(1,2)
        for jx=1+offsettest(2,1):ntest(2)-offsettest(2,2)
            % Initial Cond.
            g = @(x) sin(pi*x).*Ndiff(Ttest1,ktest(1),0,jt,0) .* Ndiff(Ttest2,ktest(2),0,jx,x);
            for stepx = 0:ktest(2)-1
                f(jx-offsettest(2,1),jt-offsettest(1,1)) = f(jx-offsettest(2,1),jt-offsettest(1,1)) + Gaussq(Ttest2(jx+stepx), Ttest2(jx+stepx+1),g,5);
            end
        end
    end
    
    
    % Solving:
    % Y = L_2xH^1
    RY1 = RieszMat(Ttest1, ktest(1) , offsettest(1,:) , 0);
    RY2 = RieszMat(Ttest2, ktest(2) , offsettest(2,:) , 1);
    RY = kron(RY1,RY2);
    %RX = sparse(RX);
    
    U = Normal_CG(B , RY , f(:) , 50000, 1e-4 * 2^(-max(ksol) * max(jsol)), zeros(size(B,2),1) );
    
    
    %% Error
    %
    % Reshape:
    u = reshape(U,nsol(2)-(offsetsol(2,1)+offsetsol(2,2)),nsol(1)-(offsetsol(1,1)+offsetsol(1,2)));
    % Preparing boundary:
    u_full = zeros(nsol(2),nsol(1));
    u_full(1+offsetsol(2,1):nsol(2)-offsetsol(2,2) , 1+offsetsol(1,1):nsol(1)-offsetsol(1,2)) = u;
    %
    
    % Plotting error to exact solution:
    exact = @(t,x) sin(pi*x).*exp(-1i*pi*pi*t);
    
    
    je = [max(max(ksol(1),ksol(2))*max(jsol(1),jsol(2))/5 , max(jsol(1),jsol(2))) , max(max(ksol(1),ksol(2))*max(jsol(1),jsol(2))/5 , max(jsol(1),jsol(2))) ];
    
    err = 0;
    g = @(t)@(x) (abs(exact(t,x) - Nev(u_full(:),{Tsol1,Tsol2},[ksol(1),ksol(2)],[t,x]))).^2;
    
    for i = 0:2^(je(1))-1
        t = i/(2^(je(1)));
        for j = 0:2^(je(2))-1
            x = j/(2^(je(2)));
            err = err + Gaussq([t,x],[t+1/(2^(je(1))) , x+1/(2^(je(2)))],g,5);
        end
    end
    e(lev) = err;
    if lev ~= 1
        quot(lev-1) = log(sqrt(e(lev-1) / e(lev)))/log(2);
    end
    
end

%save('Data_Schroedinger/err_L2','e');
%save('Data_Schroedinger/quot','quot');