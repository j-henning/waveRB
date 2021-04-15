function [X1,X2,restot, iter]=Galerkin3(A1,B1,C1,A2,B2,C2,rhs1,rhs2,maxit,tol,info)
% function [X1,X2,restot]=Galerkin3(A1,B1,C1,A2,B2,C2,rhs1,rhs2,maxit,tol)
% Solve
%        A1 X A2 + B1 X B2 + C1 X C2  = rhs1 rhs2'
% by the Galerking method
% Roughly speaking we compute a certain rational Krylov subspace and
% we impose a Galerkin condition on the residual, namely we impose the
% residual to be orthogonal to the space we construct

format short e;
V=orth(rhs1);
W=orth(rhs2);
A1r=zeros((maxit+1)*size(rhs1,2));
rhs1_projected=zeros((maxit+1)*size(rhs1,2),size(rhs1,2));
B1r=A1r;
C1r=A1r;
p1=size(V,2);
A1r(1:p1,1:p1)=V'*(A1*V); 
B1r(1:p1,1:p1)=V'*(B1*V); 
C1r(1:p1,1:p1)=V'*(C1*V);
rhs1_projected(1:p1,:)=V'*rhs1;
p1old=p1;

A2r=zeros((maxit+1)*size(rhs2,2));
rhs2_projected=zeros((maxit+1)*size(rhs2,2),size(rhs2,2));
B2r=A2r;
C2r=A2r;
p2=size(W,2);
A2r(1:p2,1:p2)=W'*(A2*W); 
B2r(1:p2,1:p2)=W'*(B2*W); 
C2r(1:p2,1:p2)=W'*(C2*W);
rhs2_projected(1:p2,:)=W'*rhs2;
p2old=p2;

s11=1;emax1=1e2;
sigma1=[emax1,emax1,s11,s11];
sigmanew1=s11;
s12=1e-2;emax2=1e7;
sigma2=[emax2,emax2,s12,s12];
sigmanew2=s12;
j=0;
restot=[];

if info
    disp('     iter        rel.res.   bckwrd er.res    sigma1      sigma2       dim(V)       dim(W)')
end
tol_trunc=1e-8;
for k=1:maxit
        iter = k; % JH

       % increase left space
       j=j+1; 
       Vnew = [(-C1-sigmanew1*A1)\V(:,k),(-B1-sqrt(sigmanew1)*A1)\V(:,k)];
       Vnew = Vnew - V*(V'*Vnew); Vnew = Vnew - V*(V'*Vnew);
       [uu,ss,vv] = svd(Vnew,0); ns=sum(diag(ss)/ss(1,1)>tol_trunc); Vnew=uu(:,1:ns);
       
       % compute the projection of the coefficient matrices
       p1=size(Vnew,2);
       A1r(1:p1old,p1old+1:p1old+p1)=V'*(A1*Vnew);
       A1r(p1old+1:p1old+p1,1:p1old)=(Vnew'*A1)*V;
       A1r(p1old+1:p1old+p1,p1old+1:p1old+p1)=(Vnew'*A1)*Vnew;
       
       B1r(1:p1old,p1old+1:p1old+p1)=V'*(B1*Vnew);
       B1r(p1old+1:p1old+p1,1:p1old)=(Vnew'*B1)*V;
       B1r(p1old+1:p1old+p1,p1old+1:p1old+p1)=(Vnew'*B1)*Vnew;
       
       C1r(1:p1old,p1old+1:p1old+p1)=V'*(C1*Vnew);
       C1r(p1old+1:p1old+p1,1:p1old)=(Vnew'*C1)*V;
       C1r(p1old+1:p1old+p1,p1old+1:p1old+p1)=(Vnew'*C1)*Vnew;
       
       V = [V, Vnew];
       nV=size(V,2);
       p1old=nV;
       
       %A1r=V'*A1*V; B1r=V'*(B1*V); C1r=V'*(C1*V); 

       % increase right space
       Wnew = [(-A2-sigmanew2*C2)\W(:,k),(-B2-sqrt(sigmanew2)*C2)\W(:,k)];
       Wnew = Wnew - W*(W'*Wnew); Wnew = Wnew - W*(W'*Wnew);
       [uu,ss,vv] = svd(Wnew,0); ns=sum(diag(ss)/ss(1,1)>tol_trunc); Wnew=uu(:,1:ns);
       
       % compute the projection of the coefficient matrices
       p2=size(Wnew,2);
       A2r(1:p2old,p2old+1:p2old+p2)=W'*(A2*Wnew);
       A2r(p2old+1:p2old+p2,1:p2old)=(Wnew'*A2)*W;
       A2r(p2old+1:p2old+p2,p2old+1:p2old+p2)=(Wnew'*A2)*Wnew;
       
       B2r(1:p2old,p2old+1:p2old+p2)=W'*(B2*Wnew);
       B2r(p2old+1:p2old+p2,1:p2old)=(Wnew'*B2)*W;
       B2r(p2old+1:p2old+p2,p2old+1:p2old+p2)=(Wnew'*B2)*Wnew;
       
       C2r(1:p2old,p2old+1:p2old+p2)=W'*(C2*Wnew);
       C2r(p2old+1:p2old+p2,1:p2old)=(Wnew'*C2)*W;
       C2r(p2old+1:p2old+p2,p2old+1:p2old+p2)=(Wnew'*C2)*Wnew;
       
       W = [W, Wnew];
       nW=size(W,2);
       p2old=nW;
       
       rhsr=rhs1_projected(1:p1old,:)*rhs2_projected(1:p2old,:)'; 

       % solve the projected problem
       T=kron(A2r(1:p2old,1:p2old)',A1r(1:p1old,1:p1old))+kron(B2r(1:p2old,1:p2old)',B1r(1:p1old,1:p1old))...
           +kron(C2r(1:p2old,1:p2old)',C1r(1:p1old,1:p1old)); 
       y = T\rhsr(:);
       Y = reshape(y,nV,nW);

       % compute the residual norm
       fact1=[A1*(V*Y),B1*(V*Y),C1*(V*Y),-rhs1];
       [~,R1]=qr(fact1,0);
       fact2=[A2'*W,B2'*W,C2'*W,rhs2];
       [~,R2]=qr(fact2,0);
       
       res_bw=norm(R1*R2','fro')/(norm(Y,'fro')*norm(A1,'fro')*norm(A2,'fro')+norm(rhs1)*norm(rhs2));
       if info
           res_rel=norm(R1*R2','fro')/(norm(rhs1)*norm(rhs2));
           disp([k,res_rel,res_bw,sigmanew1,sigmanew2,nV,nW])
       end
       
       if res_bw < tol
           break
       end

       % Compute the new shift for the left space
       eH=-eig(C1r(1:p1old,1:p1old));
       eHpoints = sort([emax1; s11; -real(eH)]);
       gs=kron(sigma1(2:end),ones(1,nW))';
       gs=sigma1(2:length(eH)+1)';
       sigmanew1 = real(newpolei(eHpoints,eH,gs)); %JH
       sigma1=[sigma1,sigmanew1,sigmanew1];

       % Compute the new shift for the right space
       eH=-eig(A2r(1:p2old,1:p2old));
       eHpoints = sort([emax2; s12; -real(eH)]);
       gs=kron(sigma2(2:end),ones(1,nW))';
       gs=sigma2(2:length(eH)+1)';
       sigmanew2 = real(newpolei(eHpoints,eH,gs)); %JH
       sigma2=[sigma2,sigmanew2,sigmanew2];

end
X1=V*Y; X2=W;


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r=ratfun(x,eH,s)

for j=1:length(x)
r(j)=abs(prod( (x(j)-s)./(x(j)-eH) ));
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r=ratfuni(x,eH,s)

for j=1:length(x)
r(j)=1./abs(prod( (x(j)-s)./(x(j)-eH) ));
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function snew=newpolei(eHpoints,eH,s)

for j=1:length(eHpoints)-1
%   snew(j) = fminbnd( @(x) ratfuni(x,eH,s), eHpoints(j),eHpoints(j+1),optimset('TolX',1e-3));
   %sval=logspace(log10(eHpoints(j)),log10(eHpoints(j+1)),200);
     sval=linspace(eHpoints(j),eHpoints(j+1),200);
    [sf,jx] = max (abs(ratfun(sval,eH,s)));
    snew(j)=sval(jx);
end
[sn,jx]=max(abs(ratfun(snew,eH,s)));
snew=snew(jx);
return

