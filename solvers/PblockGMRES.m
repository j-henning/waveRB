function [Z1,Z2,res,mem,orthog,threshold]=PblockGMRES(A,B,C1,C2,m,tol,sigmanL,precond,info,opts)
% function [Z1,Z2,res]=PblockGMRES(A,B,C1,C2,m,tol)
%
% Solve linear matrix equations of the form
% 
%     \sum_{i=1}^p A_i * X * B_i = C_1C2^T
%
% by projection onto the block Krylov subspace 
% K_m=span{vec(C_1*C_2^T),L(vec(C_1*C_2^T)),...,L^{m-1}(vec(C_1*C_2^T))}
% where L=\sum_{i=1}^p B_i^T \otimes A_i   (GMRES)
%
% INPUT
%
% A, B   cells containing the (non) sym coeff.matrices
% C1, C2   block matrices of the rhs
% m   max number of matrix-matrix multiplies
% tol final residual accuracy (backward error)
% normL 2-norm of the operator L
%
% OUTPUT:
%  Z1, Z2   such that Z1*Z2^T is the approx solution
%           (Z is truncated to sing.val >1e-12)
%  res      history of the relative residual norms 

p=length(A);
n=size(C1,1);
n2=size(C2,1);
if p~=length(B)
    error('The number of coefficient matrices in A and B must be the same.\n')
end
if ~exist('opts','var')
    opts.issym=0;
    opts.method='qrsvd';
end

% Go for Lanczos
k_max=2;
% Check if there's any nonsymmetric coefficient matrix
for i=1:p
    normAcoeff=norm(A{i}-A{i}',1);
    normBcoeff=norm(B{i}-B{i}',1);
    
    if normAcoeff>1e-14 || normBcoeff>1e-14
        % in this case switch to full Arnoldi
        k_max = m;
        break
    end
end
k_max=m;


% compute the initial residual norm
beta=sqrt(trace((C1'*C1)*(C2'*C2)));

%preallocation
omega=eye(m+1);
e=beta*eye(m+1,1);
R=zeros(m+1);
H=zeros((m+1),m);
res=ones(1,m);
restilde=ones(1,m);

%%% Preallocation of the matrices for the basis
%%% We set the maximum memory allocation equal to 5 vectors of length n^2 
V_left=spalloc(n,1*n,n);
V_right=spalloc(n2,1*n2,n2);
index=zeros(m,1);

if precond.flex
    V_left_tilde=V_left;
    V_right_tilde=V_right;
    index_tilde=index;
end
index(2)=size(C1,2);
V_left(:,1:index(2))=C1/sqrt(beta);
V_right(:,1:index(2))=C2/sqrt(beta);


% we need to incorporate the norm of L here (few steps of a power method?)
epsilon=sigmanL*tol/m; 
threshold=zeros(m,1);
%mem=zeros(m,1);
%mem(1)=2*size(C1,2);
%bdim(1)=mem(1);
%meminner=zeros(m,1);
threshold(1)=epsilon/beta; 
epsilon_orth=min(threshold(1),tol/m);
    

if strcmp(precond.method,'Sylv_precond')
    BB=(B{2}-B{3}-B{4})/B{1};
    params_precond.maxit=10;
    params_precond.tol=1e-10;
    params_precond.ch=0;   % set to 1 if spectrum is complex
        % estimates for (real) spectral interval of (A_space,M_space)
    params_precond.emin=precond.init_shift_min; 
    params_precond.emax=precond.init_shift_max;
    params_precond.Mtime=B{1};
elseif strcmp(precond.method,'GLEK')
    tolGLEK.outer=1e-2;
    tolGLEK.inner=1e-5;
    tolGLEK.trunc=1e-3;
    maxitGLEK.outer=2;
    maxitGLEK.inner=1;
    method.res=0;
    method.rhs=0;
elseif strcmp(precond.method,'mean-based')
    L_P=chol(A{1},'lower');
elseif strcmp(precond.method,'Ullmann')
    L_P=chol(A{1},'lower');
    G=B{1};
    den=trace(A{1}'*A{1});
    ll=length(B);
    if ll>1
        for i=2:length(B)
            G=G+trace(A{i}'*A{1})/den*B{i};
        end
    end
    L_G=chol(G,'lower');
end

for j=1:m
    
    
    if strcmp(precond.method,'inner-outer')
        % inner-outer preconditiener
        [rhsP1,~,rhsP2,et]=trunc(full(V_left(:,index(j)+1:index(j+1))),eye(index(j+1)-index(j)),full(V_right(:,index(j)+1:index(j+1))),1e-5,opts);
        [V1_tilde,V2_tilde,~,meminner(j+1)]=...
            blockGMRES(A,B,rhsP1,rhsP2,15,threshold(j),sigmanL,0);
        index_tilde(j+1)=index_tilde(j)+size(V1_tilde,2);
        V_left_tilde(:,index_tilde(j)+1:index_tilde(j+1))=V1_tilde;
        V_right_tilde(:,index_tilde(j)+1:index_tilde(j+1))=V2_tilde;
        
    elseif strcmp(precond.method,'Sylv_precond')

        [rhsP1,~,rhsP2,et]=trunc(-full(V_left(:,index(j)+1:index(j+1))),eye(index(j+1)-index(j)),full(V_right(:,index(j)+1:index(j+1))),1e-5,opts);
        
        %rhsP1=full(-V_left(:,index(j)+1:index(j+1)));
        %rhsP2=full(V_right(:,index(j)+1:index(j+1)));
        
      %  [V1_tilde,V2_tilde,res_inner,meminner(j+1)]=...
        %    kpik_sylv(precond.P1,L_P1,U_P1,precond.P2',L_P2,U_P2,rhsP1,rhsP2,10,tol); 
        [V1_tilde,V2_tilde,res_precond]=arn_adapt_SylvMRHS_oneside(-A{1},A{2},...
            -BB,rhsP1,B{1}'\rhsP2,params_precond);
        index_tilde(j+1)=index_tilde(j)+size(V1_tilde,2);
        V_left_tilde(:,index_tilde(j)+1:index_tilde(j+1))=-V1_tilde;
        V_right_tilde(:,index_tilde(j)+1:index_tilde(j+1))=V2_tilde.';
        
    elseif strcmp(precond.method,'GLEK')
        
        [V1_tilde]=efficientGenLyap(A{1},A{3},full(V_left(:,index(j)+1:index(j+1))),tolGLEK,maxitGLEK,method);
        
        index_tilde(j+1)=index_tilde(j)+size(V1_tilde,2);
        V_left_tilde(:,index_tilde(j)+1:index_tilde(j+1))=V1_tilde;
        V_right_tilde(:,index_tilde(j)+1:index_tilde(j+1))=V1_tilde;
    
    elseif strcmp(precond.method,'mean-based')
        % rank preservation
        index_tilde(j+1)=index(j+1);
        %V_left_tilde(:,index_tilde(j)+1:index_tilde(j+1))=L_P'\(L_P\V_left(:,index(j)+1:index(j+1)));
        %V_right_tilde(:,index_tilde(j)+1:index_tilde(j+1))=V_right(:,index(j)+1:index(j+1));
        V_left_tilde=L_P'\(L_P\V_left(:,index(j)+1:index(j+1)));
        V_right_tilde=V_right(:,index(j)+1:index(j+1));
        
    elseif strcmp(precond.method,'Ullmann')
        % rank preservation
        index_tilde(j+1)=index(j+1);
        
        %V_left_tilde(:,index_tilde(j)+1:index_tilde(j+1))=L_P'\(L_P\V_left(:,index(j)+1:index(j+1)));
        %V_right_tilde(:,index_tilde(j)+1:index_tilde(j+1))=L_G'\(L_G\V_right(:,index(j)+1:index(j+1)));
        V_left_tilde=L_P'\(L_P\V_left(:,index(j)+1:index(j+1)));
        V_right_tilde=L_G'\(L_G\V_right(:,index(j)+1:index(j+1)));
    end
    
    % multiply the last basis vector (in low-rank form) by L
    s=index_tilde(j+1)-index_tilde(j);
    left=zeros(n,p*s);
    right=zeros(n2,p*s);

    for i=1:p
        
        if precond.flex
            left(:,(i-1)*s+1:i*s)= A{i}*V_left_tilde(:,index_tilde(j)+1:index_tilde(j+1));
            right(:,(i-1)*s+1:i*s)= B{i}'*V_right_tilde(:,index_tilde(j)+1:index_tilde(j+1));
        else
            left(:,(i-1)*s+1:i*s)= A{i}*V_left_tilde;
            right(:,(i-1)*s+1:i*s)= B{i}'*V_right_tilde;
        end
    end
    
    % left*right' corresponds to the exact multiplication
    % L(vec(VV1{j}*VV2{j}'))
    % we now truncate left*right'
    [V1,~,V2,err]=trunc(left,1,right,threshold(j),opts);
    
    if info>1
        fprintf('It: %3d, MVdim %d, MV trunctol %10.5e, error %2.1e, truncdim %d \n',j,size(left,2),threshold(j),err,size(V1,2))
    end
    
    % ortogonalize V1 and V2 w.r.t the previous basis vectors (re-orth modified gram-schmidt)
    %basesize_inc=size(V1,2)+sum(bdim(max(1,j-k_max):j));
    for l=1:2        
        k_min=max(1,j-k_max);
        Theta=eye(index(j+1)+size(V1,2));
        
        for kk=k_min:j
            gamma=trace((V1'*V_left(:,index(kk)+1:index(kk+1)))*(V_right(:,index(kk)+1:index(kk+1))'*V2));   
            H(kk,j)=H(kk,j)+gamma;
            Theta(size(V1,2)+index(kk)+1:size(V1,2)+index(kk+1),size(V1,2)+index(kk)+1:size(V1,2)+index(kk+1))=-gamma*eye(index(kk+1)-index(kk));
        end
        [V1,~,V2,et]=trunc([V1,V_left(:,1:index(kk+1))],Theta,[V2,V_right(:,1:index(kk+1))],epsilon_orth,opts);
        %pause
        if info>1
            err=err+et; 
            fprintf('It: %3d, dim %d, orth trunctol %10.5e, error %2.1e, truncdim %d \n',j,basesize_old+basesize_inc,epsilon_orth,err,basesize_old+size(V1,2)) 
        end
    end
    
    if info>1
        err=err+et; 
        fprintf('It: %3d, dim %d, orth trunctol %10.5e, error %2.1e, truncdim %d \n',j,basesize_old+basesize_inc,epsilon_orth,err,basesize_old+size(V1,2)) 
    end


    
    % orthonormalize the last basis vector
    if (j<=m)        
        index(j+2)=index(j+1)+size(V1,2);
        norm_V1V2T=sqrt(trace((V1'*V1)*(V2'*V2)));
        V_left(:,index(j+1)+1:index(j+2))=V1/sqrt(norm_V1V2T);
        V_right(:,index(j+1)+1:index(j+2))=V2/sqrt(norm_V1V2T);
        H(j+1,j)=norm_V1V2T;
    end
    
    % apply the Givens rotations to the new column of H
    R(1:j+1,j) = omega(1:j+1,1:j+1)*H(1:j+1,j);
    
    % compute the next rotation and apply it also to the rhs
    GG=givens(R(j,j),R(j+1,j));
    omega(j:j+1,1:j+1)= GG *omega(j:j+1,1:j+1);  
    R(j:j+1,j)= GG *R(j:j+1,j);
    e(j:j+1)= GG *e(j:j+1);

    % solve the linear systems: we need the entire y for computing the
    % upper bound on th true residual norm. In standard implementation of
    % FOM with no truncation we can solve only one linear system once we
    % converge
    y=R(1:j,1:j)\e(1:j);
    
    % Updated inexact residual norm 
    restilde(j+1)=abs(e(j+1));
    
    % Update the upper bound on the true relative residual norm
    res(j+1)=(restilde(j+1)+abs(y)'*threshold(1:j)+j/m*tol*norm(y,1))/beta;
     
    % update the threshold for the truncation
    threshold(j+1)=epsilon/(2*restilde(j+1));
    %bdim(j+1)= index(j+2);
    
    epsilon_orth=min(threshold(j+1),tol/((j+1)*m));
    if info>1
        fprintf('It: %3d, fullbasedim %d, orth trunctol %10.5e, error %2.1e, truncdim %d \n',j,sum(bdim(max(1,j-k_max):j))+basesize_inc,epsilon_orth,err,sum(bdim(max(1,j+1-k_max):j+1)))
    
    end
    if info
        if precond.flex
            fprintf('It: %3d, Rel. Res: %10.5e, Rank Precon. basis %d, Rank Unprecond basis: %d \n',j,res(j+1),index_tilde(j+1),index(j+2))
        else
            fprintf('It: %3d, Rel. Res: %10.5e, Rank basis: %d \n',j,res(j+1),index(j+2))
        end
       
       %  fprintf('It: %3d, Rel. Res: %10.5e, Rank of the whole basis: %3d\n',j,res(j+1),bdim(j+1))

    end

  if (res(j+1)<tol) %|| abs((res(j)-res(j+1)))/res(j+1)<1e-3
     break
  end

end

% % compute the low-rank factors of the solution
if precond.flex
    Theta=eye(index_tilde(j+1));
    for kk=1:j
        Theta(index_tilde(kk)+1:index_tilde(kk+1),index_tilde(kk)+1:index_tilde(kk+1))=y(kk)*eye(index_tilde(kk+1)-index_tilde(kk));
    end
    [Z1,~,Z2,~]=trunc(V_left_tilde(:,1:index_tilde(j+1)),Theta,V_right_tilde(:,1:index_tilde(j+1)),1e-14,opts);
else
    Theta=eye(index(j+1));
    for kk=1:j
        Theta(index(kk)+1:index(kk+1),index(kk)+1:index(kk+1))=y(kk)*eye(index(kk+1)-index(kk));
    end
    if strcmp(precond.method,'Ullmann')
       [Z1,~,Z2,~]=trunc(L_P'\(L_P\V_left(:,1:index(j+1))),Theta,L_G'\(L_G\V_right(:,1:index(j+1))),1e-14,opts);
    elseif strcmp(precond.method,'mean-based')
       [Z1,~,Z2,~]=trunc(L_P'\(L_P\V_left(:,1:index(j+1))),Theta,V_right(:,1:index(j+1)),1e-14,opts);
    end
           
end

res=res(1:j+1);
threshold=threshold(1:j+1);
rank_basis=index(2:j+2)-index(1:j+1);
mem=rank_basis;
memtotal=2*rank_basis;
if precond.flex
    rank_precondbasis=index_tilde(2:j+1)-index_tilde(1:j);
    memtotal=memtotal+2*[0;rank_precondbasis];
    mem=[mem,[0;rank_precondbasis]];
end
    
mem_demand=sum(memtotal);
    

if info
    figure
    subplot(1,2,1)
    semilogy(res,'*-')
    hold on
    semilogy(threshold,'r--o')
    legend('Rel. Res.', '||E_k||')
    xlabel('Iteration k')
    ylabel('Magnitude')
    subplot(1,2,2)
    semilogy(rank_basis,'*-')
    hold on 
    xlabel('Iteration k')
    ylabel('Rank of the kth basis vector')
    semilogy(mem_demand*ones(1,j+1),'--k')
    title(['N. vects stored: ',num2str(mem_demand)])    
    legend('Unpreconditioned basis','Overall memory requirements')
    if precond.flex
        hold on
        semilogy(rank_precondbasis,'g.-')
        legend('Unpreconditioned basis','Overall memory requirements','Preconditioned basis')
    end
    
end
%%% check the orthogonality of the last basis vectors with respect all the
%%% previous ones
orthog=zeros(j+1,1);
for i=1:j+1
    %orthog(i)=trace((VV1{i}'*VV1{j+1})*(VV2{j+1}'*VV2{i})); 
    orthog(i)=trace((V_left(:,index(i)+1:index(i+1))'*V_left(:,index(j+1)+1:index(j+2)))*...
            (V_right(:,index(j+1)+1:index(j+2))'*V_right(:,index(i)+1:index(i+1))));       
end
figure
semilogy(abs(orthog))
title('Orthogonaloty between the last basis vector and the previous one')
return
end


