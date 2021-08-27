function x=pcg_fun1(a,b,x,params,maxit,tol)
%function x=pcg_fun1(a,b,x,params,maxit,tol)
% 
%  PCG with functional preconditioning. 
%  a  =  coeff. matrix
%  b  =  rhs
%  x  =  starting approx soln (e.g., x=0)
%  maxit = max number of its. (e.g., maxit=size(A,1))
%  tol = stopping tolerance for the residual (e.g., tol=1e-8)
%



n=length(b);
beta=0; x=zeros(n,1);
r=b; ap=x;
gamma=r.'*r;
res0=norm(r);
res=res0;
AP=[];
R=[];
k=0;
Ls=chol(params.M_space,'lower');
Lt=chol(params.M_time,'lower');
G=Ls\params.Q_space/Ls'; G=(G+G')/2;
[X1,E1]=eig(full(G));
G=Lt\params.Q_time/Lt'; G=(G+G')/2;
[X2,E2]=eig(full(G));
LL=1./(diag(E1)+diag(E2)');

%X=Ls'\( X1*( (X1'*(Ls\F/Lt')*X2).*LL)*X2' )/Lt;

nh=sqrt(n);



while (res/res0 > tol & k<maxit)

  %z=P(r);
  R=reshape(r,nh,nh);
  z=reshape(Ls'\( X1*( (X1'*(Ls\R/Lt')*X2).*LL)*X2' )/Lt,n,1);
  k=k+1;
  gamma=r.'*z;
  if (k==1), p=z;else, beta=gamma/gamma0;p=z+beta*p;end
  ap=a*p;
  delta=p.'*ap;
  alfa = gamma/delta;
  x = x + alfa*p;
  r = r - alfa*ap;
  gamma0=gamma;
  res(k)=norm(r);

end	
 fprintf('CG its %d  rel.res %d\n',k,res(end)/res0)
 res

