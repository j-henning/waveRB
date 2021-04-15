function X=pcg_fun2(funA,b,x,params,maxit,tol,n1,n2,flagexact)
%function X=pcg_fun1(funA,b,x,params,maxit,tol)
% 
%  PCG with functional preconditioning. 
%  funA  =  functional of coeff. matrix
%  b  =  rhs of size (n1*n2) x 1
%  x  =  starting approx soln (e.g., x=0), (n1*n2) x 1
%  maxit = max number of its. (e.g., maxit=length(b))
%  tol = stopping tolerance for the residual (e.g., tol=1e-5)
%  n1, n2 = size matrix problem
%
% Preconditioner:  Z = f(F)
%   Solve    AA Z + Z BB = F    AA = Ls\params.Q_space/Ls', 
%                               BB = Lt\params.Q_time/Lt'
%                                F = Ls\F/Lt'
% (Ls=chol(params.M_space,'lower'); Lt=chol(params.M_time,'lower');
%
%   flagexact=1, use explicit Bartels-Stewart
%   flagexact=0  use iterative solver (rational Krylov)


%path(path,'/home/valeria/matlab/Lyap/') % JH

n=length(b);
b=reshape(b,n1,n2);
X=reshape(x,n1,n2);
beta=0; x=zeros(n1,n2);
R=b; AP=X;
gamma=sum(sum(R.*R));
res0=sqrt(gamma);
res=full(res0);
normterm=norm(params.M_time,1)*norm(params.Q_space,1)+norm(params.M_space,1)*norm(params.Q_time,1);
backward_error=res0;

k=0;
Ls=chol(params.M_space,'lower');
Lt=chol(params.M_time,'lower');
G=Ls\params.Q_space/Ls'; AA=(G+G')/2;
G=Lt\params.Q_time/Lt'; BB=(G+G')/2;
  if (flagexact)
    [X1,E1]=eig(full(AA));
    [X2,E2]=eig(full(BB));
    LL=1./(diag(E1)+diag(E2)');
    V1=X1'/Ls;  V2=Lt'\X2;
    U1=Ls'\ X1; U2=X2'/Lt;
  else
    eAmin=eigs(AA,1,'sm','Tolerance',1e-2); eAmax=eigs(AA,1,'lm','Tolerance',1e-2);
    eBmin=eigs(BB,1,'sm','Tolerance',1e-2); eBmax=eigs(BB,1,'lm','Tolerance',1e-2);
  [eAmin,eAmax,eBmin,eBmax] 
  end

%X=Ls'\( X1*( (X1'*(Ls\F/Lt')*X2).*LL)*X2' )/Lt;


while (res/backward_error > tol & k<maxit)

  %z=P(r);
  if (flagexact)
     Z=U1*( (V1*R*V2).*LL)*U2;
    %Z=Ls'\( X1*( (X1'*(Ls\R/Lt')*X2).*LL)*X2' )/Lt;
  else
    [uu,ss,vv]=svds(Ls\R/Lt',1); R1=uu*sqrt(ss); R2=vv*sqrt(ss);
    % solve -AZ-ZB+R=0
    [Z1,Z2,~,~]=arn_adapt_SylvMRHS(-AA,-BB,R1,R2,n1,2e-5,eAmin,eAmax,eBmin,eBmax,0,1);
    Z=(Ls'\Z1)*(Z2'/Lt);
  end
  k=k+1;
  gamma=sum(sum(R.*Z));
  if (k==1), P=Z;else, beta=gamma/gamma0;P=Z+beta*P;end
  AP=funA(P);
  delta=sum(sum(P.*AP));
  alfa = gamma/delta;
  X = X + alfa*P;
  R = b-funA(X);  %R=R - alfa*AP;
  gamma0=gamma;
  res(k)=full(sqrt(sum(sum(R.*R))));
  backward_error=res0+norm(X,1)*normterm;
  disp([k,res(k)/backward_error])
end	
 fprintf('pcg_fun2: Dim %d %d .   CG its %d  rel.res %d\n',n1,n2,k,full(res(end)/backward_error))
%res

