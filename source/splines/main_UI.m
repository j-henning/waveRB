% By C. Mollet
%
% Graphical User Interface (GUI) for Petrov-Galerkin approximations
% of parabolic PDEs using B-splines of arbitrary oder, dimension, boundary
% and first form, its homogenization and second form
% Only zero initial conditions implemented yet

format long e;
addpath('Utilities/');

%% Setting Dimension

x = inputdlg('Enter spatial dimension:','Dimension',1,{'1'});
dim = str2double(x)+1;
clear x;

%% Setting boundary offset
% [ [t] ; [x1,x2,...] ]
% [Left, Right]
offsetsol = zeros(dim,2);
offsettest = zeros(dim,2);

type=menu(strcat('Choose full space-time weak formulation:'),...
    'First (H^1)',...
    'First Homogenization (H^1_0)',...
    'Second (L_2)');
switch type
    case 1
        offsetsol(1,:)  = [0,0];
        offsettest(1,:) = [0,0];
    case 2
        offsetsol(1,:)  = [1,0];
        offsettest(1,:) = [1,0];
    case 3
        offsetsol(1,:)  = [0,0];
        offsettest(1,:) = [0,1];
end

for i=2:dim
    
    boundary=menu(strcat('Choose boundary conditions for solution space wrt. spatial dimension No. x_', num2str(i-1),':'),...
        'Neumann/Neumann',...
        'Neumann/Zero',...
        'Zero/Neumann',...
        'Zero/Zero');
    switch boundary
        case 1
            offsetsol(i,:)=[0,0];
        case 2
            offsetsol(i,:)=[0,1];
        case 3
            offsetsol(i,:)=[1,0];
        case 4
            offsetsol(i,:)=[1,1];
    end
    
    boundary=menu(strcat('Choose boundary conditions for test space wrt. spatial dimension No. x_', num2str(i-1) ,':'),...
        'Neumann/Neumann',...
        'Neumann/Zero',...
        'Zero/Neumann',...
        'Zero/Zero');
    switch boundary
        case 1
            offsettest(i,:)=[0,0];
        case 2
            offsettest(i,:)=[0,1];
        case 3
            offsettest(i,:)=[1,0];
        case 4
            offsettest(i,:)=[1,1];
    end
    
end

clear boundary;

%% Setting Order
% init
ksol  = zeros(1,dim); % order solution space
ktest = zeros(1,dim); % order test space

x = inputdlg({'Enter Order of B-Spline for the solution space wrt. time: ' , 'and for the test space:'},'Order',1,{'2','2'});
ksol(1)=str2double(x(1));
ktest(1)=str2double(x(2));

for i = 2:dim
    x = inputdlg({strcat('Enter Order of B-Spline for the solution space wrt. dimension No. x_ ',num2str(i-1),':') , 'and for the test space:'},'Order',1,{'2','2'});
    ksol(i)=str2double(x(1));
    ktest(i)=str2double(x(2));
end
clear x;

%% Setting Level
% init
jsol  = zeros(1,dim); % discretizaion level solution space
jtest = zeros(1,dim); % discretization level test space

x = inputdlg({'Enter refinement level of B-Spline for the solution space wrt. time:' , 'and for the test space:'},'Resolution',1,{'5','6'});
jsol(1)=str2double(x(1));
jtest(1)=str2double(x(2));

for i = 2:dim
    x = inputdlg({strcat('Enter refinement level of B-Spline for the solution space wrt. dimension No. x_',num2str(i-1),':') , 'and for the test space:'},'Resolution',1,{'5','5'});
    jsol(i)=str2double(x(1));
    jtest(i)=str2double(x(2));
end
% In case of first form, I take the same Resolution for the additional
% space H
clear x;

%% Init
% Solution spaces
Tsol = cell(1,dim);
for i=1:dim
    Tsol(i)={InitUniNodes(jsol(i),ksol(i))};
end

% Test spaces
Ttest = cell(1,dim);
for i=1:dim
    Ttest(i)={InitUniNodes(jtest(i),ktest(i))};
end

% (full) dimensions:
nsol = zeros(1,dim);
ntest = zeros(1,dim);
for i=1:dim
    nsol(i) = length(cell2mat(Tsol(i))) - ksol(i);
    ntest(i) = length(cell2mat(Ttest(i))) - ktest(i);
end

%% Setting PDE
if dim > 1
    x = inputdlg('Number of spatial Terms? (A+B+C+...)','PDE',1,{'1'});
    terms = str2double(x);
    clear x;
    
    % Init Matrices as Cells
    Matind = cell(dim,terms+1);
    
    % Initialize temporal Matrix according to type of formulation
    switch type
        
        case 1 %First form
            Matind(1,1) = {[1,0]};
            for n=2:terms+1
                Matind(1,i) = {[0,0]};
            end
            
        case 2 %First form Homogenization
            Matind(1,1) = {[1,0]};
            for n=2:terms+1
                Matind(1,i) = {[0,0]};
            end
            
        case 3 %Second form
            Matind(1,1) = {[0,1]};
            for i=2:terms+1
                Matind(1,i) = {[0,0]};
            end
    end
    
    for i=2:dim
        Matind(i,1) = {[0,0]};
    end
    
    % Setting spatial matrices
    for n=1:terms
        for i=2:dim
            Matind(i,n+1) = {str2double(inputdlg({[strcat('Set up matrix for term No. ',num2str(n)) sprintf('\n') strcat('Dimension No. x_ ',num2str(i-1)) sprintf('\n') 'Order of derivative wrt. solution space:'] , 'Order of derivative wrt. test space:' },'PDE',1,{'1','1'}))};
        end
        Matind(1,n+1) = {[0,0]};
    end
end
%% Riesz Map
if dim > 1
    RieszIndHelp = inputdlg({strcat('Enter order of elliptic Operator A:') , 'and the Sobolev regularity (of the domain):'},'Elliptic Order',1,{'2','1'});
    RieszInd = [str2double(RieszIndHelp{2}),str2double(RieszIndHelp{2})-str2double(RieszIndHelp{1})];
end

%% Setting RHS
func = cell(dim,1);
rhs = inputdlg(strcat('Enter right hand side in tensor format for t',' [Format: F(x)]:'),'RHS',1,{'1'});
r = strcat('@(x) ',rhs{1});
func{1} = str2func(r);
for i=2:dim
    rhs = inputdlg(strcat('Enter right hand side in tensor format for x_,',num2str(i-1),' [Format: F(x)]:'),'RHS',1,{'1'});
    r = strcat('@(x) ',rhs{1});
    func{i} = str2func(r);
end
clear rhs;


% initial
if type == 1
    if dim == 1
        h = sparse(ntest(dim)-(offsettest(dim,1)+offsettest(dim,2)),1);
        Tt = cell2mat(Ttest(dim));
        
        for i=1+offsettest(dim,1):ntest(dim)-offsettest(dim,2)
            g = @(t) func{dim}(t) * Ndiff(Tt,ktest(dim),0,i,t);
            for step = 0:ktest(dim)-1
                h(i-offsettest(dim,1)) = h(i-offsettest(dim,1)) + Gaussq(Tt(i+step),Tt(i+step+1),g,5);
            end
        end
        f = [h;0];
    else
        h = sparse(ntest(dim)-(offsettest(dim,1)+offsettest(dim,2)),1);
        Tt = cell2mat(Ttest(dim));
        
        for i=1+offsettest(dim,1):ntest(dim)-offsettest(dim,2)
            g = @(t) func{dim}(t) * Ndiff(Tt,ktest(dim),0,i,t);
            for step = 0:ktest(dim)-1
                h(i-offsettest(dim,1)) = h(i-offsettest(dim,1)) + Gaussq(Tt(i+step),Tt(i+step+1),g,5);
            end
        end
        f = h;
        
        for d = dim-1:-1:1
            h = sparse(ntest(d)-(offsettest(d,1)+offsettest(d,2)),1);
            Tt = cell2mat(Ttest(d));
            
            for i=1+offsettest(d,1):ntest(d)-offsettest(d,2)
                g = @(t) func{d}(t) * Ndiff(Tt,ktest(d),0,i,t);
                for step = 0:ktest(d)-1
                    h(i-offsettest(d,1)) = h(i-offsettest(d,1)) + Gaussq(Tt(i+step),Tt(i+step+1),g,5);
                end
            end
            
            f = kron(h,f);
        end
        % No offset for additional Cartesian product
        % 0 initial: fex = [f;(u_0,v)_H]
        fex = sparse(size(f,1)+prod(ntest(2:dim)),size(f,2));
        fex(1:size(f,1)) = f;
        f = fex;
    end
else
    h = sparse(ntest(dim)-(offsettest(dim,1)+offsettest(dim,2)),1);
    Tt = cell2mat(Ttest(dim));
    
    for i=1+offsettest(dim,1):ntest(dim)-offsettest(dim,2)
        g = @(t) func{dim}(t) * Ndiff(Tt,ktest(dim),0,i,t);
        for step = 0:ktest(dim)-1
            h(i-offsettest(dim,1)) = h(i-offsettest(dim,1)) + Gaussq(Tt(i+step),Tt(i+step+1),g,5);
        end
    end
    f = h;
    
    for d = dim-1:-1:1
        h = sparse(ntest(d)-(offsettest(d,1)+offsettest(d,2)),1);
        Tt = cell2mat(Ttest(d));
        
        for i=1+offsettest(d,1):ntest(d)-offsettest(d,2)
            g = @(t) func{d}(t) * Ndiff(Tt,ktest(d),0,i,t);
            for step = 0:ktest(d)-1
                h(i-offsettest(d,1)) = h(i-offsettest(d,1)) + Gaussq(Tt(i+step),Tt(i+step+1),g,5);
            end
        end
        
        f = kron(h,f);
    end
end

%% Build Matrix
if dim > 1
    if type == 3
        B = StiffGen(Tsol,Ttest,ksol,ktest,offsetsol,offsettest,Matind,1);
    elseif type == 2
        B = StiffGen(Tsol,Ttest,ksol,ktest,offsetsol,offsettest,Matind);
    else
        Bhelp = StiffGen(Tsol,Ttest,ksol,ktest,offsetsol,offsettest,Matind);
        B = sparse(size(Bhelp,1) + prod(ntest(2:dim)), size(Bhelp,2));
        B(1:size(Bhelp,1) , 1:size(Bhelp,2)) = Bhelp;
        % No offset for additional Cartesian product
        % 0 initial: Bex = [B;M_h,0,...,0] due to B-Spline properties 
        MH = StiffGen(Tsol(2:dim),Ttest(2:dim),ksol(2:dim),ktest(2:dim),offsetsol(2:dim,:),zeros(dim-1,2),Matind(2:dim,2:dim));
        B(size(Bhelp,1)+1:size(Bhelp,1) + prod(ntest(2:dim)), 1:size(MH,2)) ...
            = MH;
    end
else
    if type == 3
        B = (-1) * StiffMat(Tsol{1},Ttest{1},ksol(1),ktest(1),offsetsol,offsettest,[0,1]);
    elseif type == 2
        B = StiffMat(Tsol{1},Ttest{1},ksol(1),ktest(1),offsetsol,offsettest,[1,0]);
    else
        B = sparse(ntest(1)+1-sum(offsettest(1,:)),nsol(1)-sum(offsetsol(1,:)));
        Bhelp = StiffMat(Tsol{1},Ttest{1},ksol(1),ktest(1),offsetsol,offsettest,[1,0]);
        B(1:ntest(1)-sum(offsettest(1,:)),1:nsol(1)-sum(offsetsol(1,:))) = Bhelp;
        B(ntest(1)-sum(offsettest(1,:))+1,1) = 1;
    end
end

% build Riesz Matrix
% Second Form
if type == 3
    %L_2(I)
    Rt0 = RieszMat(Ttest{1},ktest(1),offsettest(1,:),0);
    % H^1(I)
    % Note: This assembles the full norm where the seminorm would be required
    % but they are equivalent
    Rt1 = RieszMat(Ttest{1},ktest(1),offsettest(1,:),1);
    if dim > 1
        % V, A:W -> V', i.e. L_2(I;V) \cap H^1(I;W')
        RInd = abs(RieszInd(1)) * ones(dim-1);
        RV = RieszGen(Ttest(2:dim) , ktest(2:dim) , offsettest(2:dim,:) , RInd);
        if RieszInd(2) > 0
            P = RieszGen(Ttest(2:dim) , ktest(2:dim) , offsettest(2:dim,:) , zeros(dim-1,1));
            RV = P * (RV\P);
        end
        % W'
        RInd = abs(RieszInd(1)) * ones(dim-1);
        RWd = RieszGen(Ttest(2:dim) , ktest(2:dim) , offsettest(2:dim,:) , RInd);
        if RieszInd(1) > 0
            P = RieszGen(Ttest(2:dim) , ktest(2:dim) , offsettest(2:dim,:) , zeros(dim-1,1));
            RWd = P * (RWd\P);
        end
        % Final Riez mapping
        R = kron(Rt0,RV) + kron(Rt1,RWd);
        
        % First Form
    else
        R = Rt1;
    end
elseif type ==2
    %L_2(I)
    Rt0 = RieszMat(Ttest{1},ktest(1),offsettest(1,:),0);
    % V (i.a. dual)
    if dim > 1
        RInd = abs(RieszInd(2)) * ones(dim-1);
        RV = RieszGen(Ttest(2:dim) , ktest(2:dim) , offsettest(2:dim,:) , RInd);
        % A: W -> V', i.e., Test Space V
        if RieszInd(2) > 0
            P = RieszGen(Ttest(2:dim) , ktest(2:dim) , offsettest(2:dim,:) , zeros(dim-1,1));
            RV = P * (RV\P);
        end
        R = kron(Rt0,RV);
    else
        R = Rt0;
    end
else
    %L_2(I)
    Rt0 = RieszMat(Ttest{1},ktest(1),offsettest(1,:),0);
    % V (i.a. dual)
    if dim > 1 
        RInd = abs(RieszInd(2)) * ones(dim-1);
        RV = RieszGen(Ttest(2:dim) , ktest(2:dim) , offsettest(2:dim,:) , RInd);
        % A: W -> V', i.e., Test Space V
        if RieszInd(2) > 0
            P = RieszGen(Ttest(2:dim) , ktest(2:dim) , offsettest(2:dim,:) , zeros(dim-1,1));
            RV = P * (RV\P);
        end
        R = kron(Rt0,RV);
        % Rex = [R,0;0,R_H] due to Cartesian product
        RH = RieszGen(Ttest(2:dim) , ktest(2:dim) , zeros(dim-1,2) , zeros(dim-1,1));
        Rex = sparse(size(R,1) + prod(ntest(2:dim)) , size(R,2) + prod(ntest(2:dim)));
        Rex(1:size(R,1) , 1:size(R,2)) = R;
        Rex(size(R,1)+1:size(R,1) + size(RH,1) , size(R,2)+1:size(R,2)+ size(RH,1)) ...
            = RH;
        R = Rex;
    else
        Rt0help = Rt0;
        Rt0 = sparse(ntest(1)+1-sum(offsettest(1,:)),ntest(1)+1-sum(offsettest(1,:)));
        Rt0(1:ntest(1)-sum(offsettest(1,:)),1:ntest(1)-sum(offsettest(1,:))) = Rt0help;
        Rt0(ntest(1)-sum(offsettest(1,:))+1,ntest(1)-sum(offsettest(1,:))+1) = 1;
        R = Rt0;
    end
    
end

% Works perfect! But theoretically not the correct construction if
% negative Sobolev Order required
%RieszMapTest = RieszGen(Ttest , ktest , offsettest , Matind);

%% Solving

U = Normal_CG(B , R , f(:) , 5000, 1e-0 * 2^( -max(ksol) * max(jsol) ) , zeros(size(B,2),1) );

%% Plot
plotres = ones(dim,1);
x = inputdlg('Enter plot resolution wrt. time:','Plot',1,{'5'});
plotres(1)=str2double(x(1));

for i = 2:dim
    x = inputdlg(strcat('Enter plot resolution wrt. dimension No. x_ ',num2str(i-1),':'),'Plot',1,{'5'});
    plotres(i)=str2double(x(1));
end
clear x;


if dim == 2
    % Standard Plot
    % Reshape in 2D
    u = reshape(U,nsol(2)-(offsetsol(2,1)+offsetsol(2,2)),nsol(1)-(offsetsol(1,1)+offsetsol(1,2)));
    % Preparing boundary for plot
    u_full = zeros(nsol(2),nsol(1));
    u_full(1+offsetsol(2,1):nsol(2)-offsetsol(2,2) , 1+offsetsol(1,1):nsol(1)-offsetsol(1,2)) = u;
    
    figure(1);
    f = plot_Spline(u_full(:),Tsol,ksol,plotres);
    title('Solution');
    xlabel('t');
    ylabel('x');
    
    % ODE
elseif dim == 1
    u_full = zeros(nsol(1),1);
    u_full(1+offsetsol(1,1):nsol(1)-offsetsol(1,2)) = U;
    
    figure(1);
    f = plot_Spline(u_full,Tsol,ksol,plotres);
    title('Solution');
    xlabel('t');
    ylabel('x');
    
    % Movie
elseif dim == 3
    u = reshape(U,nsol(3)-(offsetsol(3,1)+offsetsol(3,2)),nsol(2)-(offsetsol(2,1)+offsetsol(2,2)),nsol(1)-(offsetsol(1,1)+offsetsol(1,2)));
    % Preparing boundary for plot
    u_full = zeros(nsol(3),nsol(2),nsol(1));
    u_full(1+offsetsol(3,1):nsol(3)-offsetsol(3,2) , 1+offsetsol(2,1):nsol(2)-offsetsol(2,2) , 1+offsetsol(1,1):nsol(1)-offsetsol(1,2)) = u;
    
    figure(1);
    f = plot_Spline(u_full(:),Tsol,ksol,plotres);
    title('Solution');
    xlabel('t');
    ylabel('x')
end
