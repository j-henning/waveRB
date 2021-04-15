% By C. Mollet
%

%% Settings
format long e;

% Since Code is in a subdirectory
addpath('../');
addpath('../Utilities');

% Setting Dimension
dim = 3;

% space-time
% [ [t] ; [x] ; [y] ]
% [Left, Right]
% First Form, Zero AB
offsetsol = [[1,0];[1,1];[1,1]];
offsettest = [[1,0];[1,1];[1,1]];

for ind=1:5
    k=2;
    j=ind;

% Setting level
lev = j;
jsole = [lev,lev,lev]; % discretizaion level solution space
jteste = [lev+1,lev+1,lev+1]; % discretization level test space

% Setting order
ksol = [k,k,k]; % order solution space
ktest = [1,2,2]; % order test space

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
terms = 3;

%% Generating

% Stiff matrix
T11 = StiffMat(Tsole{1},Tteste{1},ksol(1),ktest(1),offsetsol(1,:),offsettest(1,:),[1,0]);
M21 = StiffMat(Tsole{2},Tteste{2},ksol(2),ktest(2),offsetsol(2,:),offsettest(2,:),[0,0]);
M31 = StiffMat(Tsole{3},Tteste{3},ksol(3),ktest(3),offsetsol(3,:),offsettest(3,:),[0,0]);

M12 = StiffMat(Tsole{1},Tteste{1},ksol(1),ktest(1),offsetsol(1,:),offsettest(1,:),[0,0]);
S22 = StiffMat(Tsole{2},Tteste{2},ksol(2),ktest(2),offsetsol(2,:),offsettest(2,:),[1,1]);
M32 = StiffMat(Tsole{3},Tteste{3},ksol(3),ktest(3),offsetsol(3,:),offsettest(3,:),[0,0]);

M13 = StiffMat(Tsole{1},Tteste{1},ksol(1),ktest(1),offsetsol(1,:),offsettest(1,:),[0,0]);
S23 = StiffMat(Tsole{2},Tteste{2},ksol(2),ktest(2),offsetsol(2,:),offsettest(2,:),[0,0]);
M33 = StiffMat(Tsole{3},Tteste{3},ksol(3),ktest(3),offsetsol(3,:),offsettest(3,:),[1,1]);

B = kron(kron(T11,M21),M31) + kron(kron(M12,S22),M32) + kron(kron(M13,S23),M33);

% Y = L_2xH^1
RY1 = RieszMat(Tteste{1}, ktest(1) , offsettest(1,:) , 0);
RY2_1 = RieszMat(Tteste{2}, ktest(2) , offsettest(2,:) , 1);
RY3_1 = RieszMat(Tteste{3}, ktest(3) , offsettest(3,:) , 0);

RY2_2 = RieszMat(Tteste{2}, ktest(2) , offsettest(2,:) , 0);
RY3_2 = RieszMat(Tteste{3}, ktest(3) , offsettest(3,:) , 1);

RY = kron(RY1,kron(RY2_1,RY3_1) + kron(RY2_2,RY3_2));


% RHS
% initial
% y %%
h = zeros(nteste(dim)-(offsettest(dim,1)+offsettest(dim,2)),1);
Tt = Tteste{dim};

for i=1+offsettest(dim,1):nteste(dim)-offsettest(dim,2)
    g = @(y) sin(2*pi*y).*Ndiff(Tt,ktest(dim),0,i,y);
    %g = @(x) c(k) .* sin(2*pi*x).* Ndiff(Tt,ktest(dim),0,i,x);
    for step = 0:ktest(dim)-1
        h(i-offsettest(dim,1)) = h(i-offsettest(dim,1)) + Gaussq(Tt(i+step),Tt(i+step+1),g,5);
    end
end
f = h;

% x %%
%for d = dim-1:-1:1
d=dim - 1; 
    h = zeros(nteste(d)-(offsettest(d,1)+offsettest(d,2)),1);
    Tt = Tteste{d};
    
    for i=1+offsettest(d,1):nteste(d)-offsettest(d,2)
        g = @(x) sin(2*pi*x).*Ndiff(Tt,ktest(d),0,i,x);
        %g = @(t) Ndiff(Tt,ktest(d),0,i,0);
        for step = 0:ktest(d)-1
            h(i-offsettest(d,1)) = h(i-offsettest(d,1)) + Gaussq(Tt(i+step),Tt(i+step+1),g,5);
        end
    end
    
    f = kron(h,f);
    
% t %%
d=1;
    h = zeros(nteste(d)-(offsettest(d,1)+offsettest(d,2)),1);
    Tt = Tteste{d};
    
    for i=1+offsettest(d,1):nteste(d)-offsettest(d,2)
        g = @(t) (8*pi*pi*t+1).*Ndiff(Tt,ktest(d),0,i,t);
        %g = @(t) Ndiff(Tt,ktest(d),0,i,0);
        for step = 0:ktest(d)-1
            h(i-offsettest(d,1)) = h(i-offsettest(d,1)) + Gaussq(Tt(i+step),Tt(i+step+1),g,5);
        end
    end
    
    f = kron(h,f);
%end

%% Solving

 U = Normal_CG(B , RY , f(:) , 5000000, 1e-2 * 2^( -max(ksol) * max(jsole) ) , zeros(size(B,2),1) );
 
 
  %% Error estimate
 
 % Reshape
u = reshape(U,nsole(3)-(offsetsol(3,1)+offsetsol(3,2)),nsole(2)-(offsetsol(2,1)+offsetsol(2,2)),nsole(1)-(offsetsol(1,1)+offsetsol(1,2)));
% Preparing boundary
u_full = zeros(nsole(3),nsole(2),nsole(1));
u_full(1+offsetsol(3,1):nsole(3)-offsetsol(3,2) , 1+offsetsol(2,1):nsole(2)-offsetsol(2,2) , 1+offsetsol(1,1):nsole(1)-offsetsol(1,2)) = u;


 val = 0;
 f = @(y)@(x)@(t) (abs(sin(2*pi*x)*sin(2*pi*y)*t - Nev(u_full(:),Tsole,ksol,[t,x,y]))).^2;
 t=0:2^(-lev):1-2^(-lev);
 %for t=0:2^(-lev):1-2^(-lev)
 parfor i=1:length(t)
     for x=0:2^(-lev):1-2^(-lev)
         %y=0:2^(-lev):1-2^(-lev);
         for y=0:2^(-lev):1-2^(-lev);
         %parfor i=1:length(y)
             %val = val + Gaussq([t,x,y(i)],[t+2^(-lev),x+2^(-lev),y(i)+2^(-lev)],f,4);
             val = val + Gaussq([t(i),x,y],[t(i)+2^(-lev),x+2^(-lev),y+2^(-lev)],f,4);
         end
     end
 end
        
err(ind) = sqrt(val);


end;

%save('Data_3D/L2error_k2_1_5','err');

