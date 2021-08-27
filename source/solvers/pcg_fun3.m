function X=pcg_fun1(funA,b,x,params,maxit,tol,n1,n2,flagexact)
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
% Prec operator:   P z = r,  P=(Dt x Ms + Mt x As)' inv(Mt x Ms) (Dtx Ms + Mt x As)
%
%   flagexact=1, use explicit preconditioning operator


%path(path,'/home/valeria/matlab/Lyap/') JH

n=length(b);
b=reshape(b,n1,n2);
X=reshape(x,n1,n2);
beta=0; %x=zeros(n1,n2);
  R = b-funA(X);  %R=R - alfa*AP;
%R=b; %AP=X;
gamma=sum(sum(R.*R));
res0=sqrt(gamma);
res=full(res0);
k=0;
normterm=norm(params.M_time,1)*norm(params.Q_space,1)+norm(params.M_space,1)*norm(params.Q_time,1);
backward_error=res0;


%LMs=chol(params.M_space,'lower'); 
if isfield(params,'d')
    switch params.d
        case 1
            [xxs,ees]=eig(full(params.M_spacelocal)); 
            EQs=diag(xxs'*params.Q_spacelocal*xxs);
            qq1=xxs;
            rr1=diag(sqrt(EQs)./sqrt(diag(ees)));
            S1=qq1*diag(1./sqrt(ees));

        case 2
            [xxs,ees]=eig(full(params.M_spacelocal)); 
            EQs_local=diag(xxs'*params.Q_spacelocal*xxs);
            ENs_local=diag(xxs'*params.N_spacelocal*xxs);
            EQs = kron(EQs_local, diag(ees)) + 2*kron(ENs_local, ENs_local)   + kron(diag(ees), EQs_local);
            ees=kron(diag(ees),diag(ees));
            qq1=kron(xxs,xxs);
            rr1=diag(sqrt(EQs)./sqrt(ees));
            S1=qq1*diag(1./sqrt(ees));

        case 3
            [xxs,ees]=eig(full(params.M_spacelocal)); 
            EQs_local=diag(xxs'*params.Q_spacelocal*xxs);
            ENs_local=diag(xxs'*params.N_spacelocal*xxs);
            
            EQs = kron(EQs_local,kron(diag(ees), diag(ees))) + kron(ENs_local,kron( diag(ees), ENs_local)) + kron(ENs_local,kron(ENs_local, diag(ees))) ...
               + kron(ENs_local,kron( diag(ees), ENs_local)) + kron(diag(ees), kron(diag(ees), EQs_local)) + kron(diag(ees),kron(ENs_local,ENs_local)) ...
               + kron(ENs_local,kron(ENs_local, diag(ees))) + kron(diag(ees),kron(ENs_local, ENs_local)) + kron(diag(ees), kron(EQs_local, diag(ees)));
            ees=kron(kron(diag(ees),diag(ees)),diag(ees));
            qq1=kron(kron(xxs,xxs),xxs);
            %EQs2=diag(qq1'*params.Q_space*qq1);
            rr1=diag(sqrt(EQs)./sqrt(ees));
    %       S1=qq1*diag(1./sqrt(ees));
sees=1./sqrt(ees);sees=sees(:)';
    S1=qq1.*repmat(sees,size(qq1,1),1);

    end
else
    [xxs,ees]=eig(full(params.M_space)); 
    % Davide: params.Q_space and params.M_space can be simultaneously
    % diagonalized, we can compute only one eigendecomposition
    %EQs=diag(xxs'*params.Q_space*xxs);
    EQs=sum(xxs.*(params.Q_space*xxs))';
    qq1=xxs;
    rr1=spdiags(sqrt(EQs)./sqrt(diag(ees)),0,n1,n1);
    %S1=qq1*diag(1./sqrt(diag(ees)));
    S1=qq1.*repmat(1./sqrt(diag(ees))',n1,1);

end
    
    %LMs=xxs*sqrt(ees)*xxs';
%[XQs,EQs]=eig(full(params.Q_space+params.Q_space')/2);
%[qq1_schur,rr1_schur]=schur(full(XQs*sqrt(EQs)/sqrt(ees)*xxs'));
%qq1=qq1_schur;
%rr1=rr1_schur;
%rr1(1:10,1:10)
%rr1_schur(1:10,1:10)
%norm(sort(diag(rr1))-sort(diag(rr1_schur)))
%pause
%norm(sort(diag(EQs))-sort(diag(xxs'*params.Q_space*xxs)))
%diag(EQs)
%ll=xxs'*params.Q_space*xxs;
%ll(1:10,1:10)
%norm(XQs'*xxs)
%pause

%As=LMs*XQs*sqrt(EQs)*XQs';
%eig(full(As)),pause
LQt=chol(params.Q_time,'lower'); 
%[XQt,EQt]=eig(full(params.Q_time+params.Q_time')/2);
%LQt=XQt*sqrt(EQt)*XQt';


%[XQt,EQt]=eig(full(params.Q_time+params.Q_time')/2);
%LQt=XQt*sqrt(EQt)*XQt';
LMt=chol(params.M_time,'lower');  
%[xxt,eet]=eig(full(params.M_time)); 
%LMt_inv=xxt*diag(1./sqrt(diag(eet)))*xxt';

%Dt=LMt*LQt';
%[XQt,EQt]=eig(full(params.Q_time+params.Q_time')/2);

%eig(full(params.D_time'/params.M_time*params.D_time), full(Dt'/params.M_time*Dt))

%[qq1,rr1]=schur(full(params.M_space\As)); [qq2,rr2]=schur(full(Dt'/params.M_time));
%[qq3,rr3]=schur(full(params.M_space\As')); [qq4,rr4]=schur(full(Dt/params.M_time));

%[qq1,rr1]=schur(full(LMs\As/LMs')); [qq2,rr2]=schur(full(LMt\Dt'/LMt'));




%norm(full(XQs*sqrt(EQs)*XQs'/LMs)-full(XQs*sqrt(EQs)*XQs'/LMs)',1)
%pause

[qq2,rr2]=schur(full(LMt\LQt));
%[qq2,rr2]=schur(full(LMt_inv*LQt));
%rr2(1:10,1:10)
%spy(rr2>1e-10)
%norm(full(LMt\LQt)-full(LMt\LQt)',1)
%pause

%[qq1,rr1]=eig(full(LMs\As/LMs')); [qq2,rr2]=eig(full(LMt\Dt'/LMt'));
%[qq3,rr3]=schur(full(LMs\As'/LMs')); [qq4,rr4]=schur(full(LMt\Dt/LMt'));
%qq3=qq1; 
rr3=rr1'; %qq4=qq2; 
rr4=rr2';
%S1=LMs'\qq1;
S2=qq2'/LMt;
%S2=qq2'*LMt_inv;

%S3=qq3'/LMs; 
S3=S1';
%S4=LMt'\qq4;
S4=S2';
II=eye(size(rr4));
   
%S1=LMs'\qq1; S2=inv(qq2)/LMt; S3=inv(qq3)/LMs; S4=LMt'\qq4;

while (res/backward_error > tol & k<maxit)

  %z=P(r);
  if (flagexact)

   %Z1= LMs'\(qq3*lyap(rr3,rr4,-qq3'*(LMs\R/LMt')*qq4)*qq4')/LMt;
   %Z1 = params.M_space*Z1*params.M_time;
   %Z = LMs'\(qq1*lyap(rr1,rr2,-qq1'*(LMs\Z1/LMt')*qq2)*qq2')/LMt;
   %Z1= (qq3*lyap(rr3,rr4,-qq3'*(LMs\R/LMt')*qq4)*qq4');
   %Z = LMs'\(qq1*lyap(rr1,rr2,-qq1'*Z1*qq2)*qq2')/LMt;
   %Z1= qq3*lyap(rr3,rr4,-S3*R*S4)*qq4';
   %Z = S1*lyap(rr1,rr2,-qq1'*Z1*qq2)*S2;
%   Z1= qq3*lyap(rr3,rr4,-S3*R*S4)*qq4';
    Z = S1*lyap(rr1,rr2,-lyap(rr3,rr4,-S3*R*S4))*S2;
    
%{
    % Davide: solve 2l shifted linear systems (l=# spatial nodes) instead of two
    % lyap
    rhs1=S3*R*S4;
    Z=zeros(size(rhs1));
    for ii=1:size(rr1,1)
        Z(ii,:)=(rhs1(ii,:)/(rr4+rr1(ii,ii)*II))/(rr2+rr1(ii,ii)*II);
    end
    Z=S1*Z*S2;    
%}

  else
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
   disp([k,res(k),res(k)/backward_error]) %JH

end	
 fprintf('pcg_fun3: Dim %d %d .   CG its %d  backward error %d\n',n1,n2,k,full(res(end)/backward_error))
%fprintf('pcg_fun3: Dim %d %d .   CG its %d  rel.res %d\n',n1,n2,k,full(res(end)/res0))
%res/res0

