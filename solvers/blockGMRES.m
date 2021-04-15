function [Z1,Z2,res,mem]=blockGMRES(A,B,C1,C2,m,tol,sigmanL,info)
% function [Z1,Z2,res]=blockGMRES(A,B,C1,C2,m,tol)
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
k_max=m;

% compute the initial residual norm
beta=sqrt(trace((C1'*C1)*(C2'*C2)));

%preallocation
%VV1=cell(m);
%VV2=cell(m);
H=zeros((m+1),m);
omega=eye(m+1);
e=beta*eye(m+1,1);
R=zeros(m+1);
res=ones(1,m);
restilde=ones(1,m);


% fisrt basis vector (in low-rank format)
% v1=vec(VV1{1}*VV2{1}')
%VV1{1}=C1/sqrt(beta);
%VV2{1}=C2/sqrt(beta);

%%% Preallocation of the matrices for the basis
%%% We set the maximum memory allocation equal to 5 vectors of length n^2 
V_left=zeros(n,1*n);
V_right=zeros(n2,1*n2);

index=zeros(m,1);
index(2)=size(C1,2);
V_left(:,1:index(2))=C1/sqrt(beta);
V_right(:,1:index(2))=C2/sqrt(beta);

% we need to incorporate the norm of L here (few steps of a power method?)
epsilon=sigmanL*tol/m; 
threshold=zeros(m,1);
mem=zeros(m,1);
mem(1)=size(C1,2);
threshold(1)=epsilon/beta; 
epsilon_orth=min(threshold(1),tol/m);

for j=1:m

    % multiply the last basis vector (in low-rank form) by L
%     s1=size(VV1{j},2);
%     s2=size(VV2{j},2);
%     
%     left=zeros(n,p*s1);
%     right=zeros(n2,p*s2);
%     
    s=index(j+1)-index(j);
    left=zeros(n,p*s);
    right=zeros(n2,p*s);
    
    for i=1:p
        %left(:,(i-1)*s1+1:i*s1)= A{i}*VV1{j};
        %right(:,(i-1)*s2+1:i*s2)= B{i}'*VV2{j};
        
        left(:,(i-1)*s+1:i*s)= A{i}*V_left(:,index(j)+1:index(j+1));
        right(:,(i-1)*s+1:i*s)= B{i}'*V_right(:,index(j)+1:index(j+1));
    end  
    
    % left*right' corresponds to the exact multiplication
%     % L(vec(VV1{j}*VV2{j}'))
%     % we now truncate left*right'
%     [Q1,R1]=qr(left,0);
%     R1=triu(R1);
%     [Q2,R2]=qr(right,0);
%     R2=triu(R2);
%      
%     [u,s,v]=svd(R1*R2'); 
%     [s,id]=sort(diag(s));
%     u=u(:,id(end:-1:1));
%     v=v(:,id(end:-1:1));    
%     % truncated SVD in Frobenius norm
%     is=sum(cumsum(abs(s))/sum(s)<threshold(j));
%     s=flipud(s);
%     V1=Q1*(u(:,1:end-is)*diag(sqrt(s(1:end-is))));
%     V2=Q2*(v(:,1:end-is)*diag(sqrt(s(1:end-is))));
%     
    [V1,~,V2,err]=trunc(left,1,right,threshold(j),opts);
    
    % ortogonalize V1 and V2 w.r.t the previous basis vectors (re-orth modified gram-schmidt) 
    for l=1:2
        k_min=max(1,j-k_max);
         Theta=eye(index(j+1)+size(V1,2));
         
        for kk=k_min:j
            %gamma=trace((V1'*VV1{kk})*(VV2{kk}'*V2));            
            %H(kk,j)=H(kk,j)+gamma;
            %left=[V1 VV1{kk}];
            %right=[V2 VV2{kk}];
            %[Q1,R1]=qr(left,0);
            %R1=triu(R1);
            %[Q2,R2]=qr(right,0);
            %R2=triu(R2);
            
            gamma=trace((V1'*V_left(:,index(kk)+1:index(kk+1)))*(V_right(:,index(kk)+1:index(kk+1))'*V2));   
            H(kk,j)=H(kk,j)+gamma;
             
            %Theta=eye(size(left,2));
            %Theta(size(V1,2)+1:end,size(V1,2)+1:end)=-gamma*eye(size(VV1{kk},2));
       
            Theta(size(V1,2)+index(kk)+1:size(V1,2)+index(kk+1),size(V1,2)+index(kk)+1:size(V1,2)+index(kk+1))=-gamma*eye(index(kk+1)-index(kk));
            
%             [u,s,v]=svd(R1*Theta*R2'); 
%             [s,id]=sort(diag(s));
%             u=u(:,id(end:-1:1));
%             v=v(:,id(end:-1:1));    
%             % truncated SVD in Frobenius norm
%             is=sum(cumsum(abs(s))/sum(s)<epsilon_orth);
%             s=flipud(s);
%     
%             V1=Q1*(u(:,1:end-is)*diag(sqrt(s(1:end-is))));
%             V2=Q2*(v(:,1:end-is)*diag(sqrt(s(1:end-is))));
%             
         %   err=0;
         %   [V1,~,V2,et]=trunc(left,Theta,right,epsilon_orth,opts);
        end
        err=0;
        
%         if l==2
%         [Q1,R1]=qr([V1,V_left(:,1:index(kk+1))],0);
%         R1=triu(R1);
%         [Q2,R2]=qr([V2,V_right(:,1:index(kk+1))],0);
%         R2=triu(R2);
%         [u,s,v]=svd(R1*Theta*R2');
%         [s,id]=sort(diag(s));
%         u=u(:,id(end:-1:1));
%         v=v(:,id(end:-1:1));    
%         % truncated SVD in Frobenius norm
%         is=sum(cumsum(abs(s))/sum(s)<epsilon_orth);
%         s=flipud(s);
%        
%     
%         V1prova=Q1*(u(:,1:end-is)*diag(sqrt(s(1:end-is))));
%         V2prova=Q2*(v(:,1:end-is)*diag(sqrt(s(1:end-is))));
%         norm_V1V2T=sqrt(trace((V1prova'*V1prova)*(V2prova'*V2prova)));
%         V1prova=V1prova/sqrt(norm_V1V2T);
%         V2prova=V2prova/sqrt(norm_V1V2T);
%         
%         for kk=k_min:j
%             trace((V1prova'*V_left(:,index(kk)+1:index(kk+1)))*(V_right(:,index(kk)+1:index(kk+1))'*V2prova))        
%             pause
%         end    
%         end    
        
        
        [V1,~,V2,et]=trunc([V1,V_left(:,1:index(kk+1))],Theta,[V2,V_right(:,1:index(kk+1))],epsilon_orth,opts);
        
       
    end
   
 
    % orthonormalize the last basis vector
    if (j<=m)
        %norm_V1V2T=norm(s(1:end-is));        
%         norm_V1V2T=sqrt(trace((V1'*V1)*(V2'*V2)));
%         VV1{j+1}=V1/sqrt(norm_V1V2T);
%         VV2{j+1}=V2/sqrt(norm_V1V2T);
%         H(j+1,j)=norm_V1V2T;
        
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
    % GMRES with no truncation we can solve only one linear system once we
    % converge
    y=R(1:j,1:j)\e(1:j);
    
    % compute the 'inexact' residual
    restilde(j+1)=abs(e(j+1));
    
    % Update the upper bound on the true relative residual norm
    res(j+1)=(restilde(j+1)+abs(y)'*threshold(1:j)+j/m*tol*norm(y,1))/beta;
    
    
    
    % update the threshold for the truncation
    %sigma=min(svd(H(1:j,1:j)));
    threshold(j+1)=epsilon/(2*restilde(j+1));
   % mem(j+1)= size(VV1{j+1},2);
    epsilon_orth=min(threshold(j+1),tol/((j+1)*m));
    
    if info
        %fprintf('It: %3d, Rel. Res: %10.5e, Rank current basis vector: %3d\n',j,res(j+1),mem(j+1))
        fprintf('It: %3d, Rel. Res: %10.5e, Rank of the whole basis: %3d\n',j,res(j+1),index(j+1))
    end

        
  if (res(j+1)<tol) %|| (res(j)-res(j+1))/res(j+1)<1e-3
     break
  end

end


Theta=eye(index(j+1));
for kk=1:j
    Theta(index(kk)+1:index(kk+1),index(kk)+1:index(kk+1))=y(kk)*eye(index(kk+1)-index(kk));
end
[Z1,~,Z2,err]=trunc(V_left(:,1:index(j+1)),Theta,V_right(:,1:index(j+1)),1e-12,opts);

% % compute the low-rank factors of the solution
% Z1=[];
% Z2=[];
% for i=1:j
%     left=[Z1 VV1{i}];
%     right=[Z2 VV2{i}];
%     %[Q1,R1]=qr(left,0);
%     %R1=triu(R1);
%     %[Q2,R2]=qr(right,0);
%     %R2=triu(R2);
%     Theta=eye(size(left,2));
%     Theta(size(Z1,2)+1:end,size(Z1,2)+1:end)=y(i)*eye(size(VV1{i},2));
%     
%     %[u,s,v]=svd(R1*Theta*R2'); 
%     %[s,id]=sort(diag(s));
%     %u=u(:,id(end:-1:1));
%     %v=v(:,id(end:-1:1));    
%     % truncated SVD in Frobenius norm
%     %is=sum(cumsum(abs(s))/sum(s)<tol);
%     %s=flipud(s);
%     
%     %Z1=Q1*(u(:,1:end-is)*diag(sqrt(s(1:end-is))));
%     %Z2=Q2*(v(:,1:end-is)*diag(sqrt(s(1:end-is))));
%     
%      [Z1,~,Z2,err]=trunc(left,Theta,right,1e-12,opts);
% end

res=res(1:j+1);
threshold=threshold(1:j+1);
%mem=2*mem(1:j+1);
mem=index(1:j+1)-[0;index(1:j)];

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
    semilogy(mem,'*-')
    hold on 
    vect=p.^(0:j)*size(C1,2);
    index2=find(vect>n,1,'first');
    vect(index2:end)=n;
    semilogy(vect,'o--r')
    xlabel('Iteration k')
    ylabel('Rank of the kth basis vector')
    title(['N. vects stored: ',num2str(2*index(j+1))])
    legend('Inexact','Exact')

    %%% check the orthogonality of the last basis vectors with respect all the
    %%% previous ones
    orthog=zeros(j+1,1);
    for i=1:j
        %orthog(i)=trace((VV1{i}'*VV1{j+1})*(VV2{j+1}'*VV2{i})); 
        orthog(i)=trace((V_left(:,index(i)+1:index(i+1))'*V_left(:,index(j)+1:index(j+1)))*...
            (V_right(:,index(j)+1:index(j+1))'*V_right(:,index(i)+1:index(i+1))));       
        
    end
    figure
    semilogy(abs(orthog))
    title('Orthogonaloty between the last basis vector and the previous one')
end
mem=sum(mem);

return
end
