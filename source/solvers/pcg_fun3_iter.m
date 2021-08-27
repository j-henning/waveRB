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
backward_error=1;

if (flagexact)
% Direct preconditioner

%LMs=chol(params.M_space,'lower'); 
%LQt=chol(params.Q_time,'lower');
%Lt=sqrtm(full(params.M_time));
[qq1,rr1]=schur(full(params.M_space\params.A_space));
%S1=qq1;
    

% if isfield(params,'d')
%     switch params.d
%         case 1
%             [xxs,ees]=eig(full(params.M_spacelocal)); 
%             EQs=diag(xxs'*params.Q_spacelocal*xxs);
%             qq1=xxs;
%             rr1=diag(sqrt(EQs)./sqrt(diag(ees)));
%             S1=qq1*diag(1./sqrt(ees));
% 
%         case 2
%             [xxs,ees]=eig(full(params.M_spacelocal)); 
%             EQs_local=diag(xxs'*params.Q_spacelocal*xxs);
%             ENs_local=diag(xxs'*params.N_spacelocal*xxs);
%             EQs = kron(EQs_local, diag(ees)) + 2*kron(ENs_local, ENs_local)   + kron(diag(ees), EQs_local);
%             ees=kron(diag(ees),diag(ees));
%             qq1=kron(xxs,xxs);
%             rr1=diag(sqrt(EQs)./sqrt(ees));
%             S1=qq1*diag(1./sqrt(ees));
% 
%         case 3
%             [xxs,ees]=eig(full(params.M_spacelocal)); 
%             EQs_local=diag(xxs'*params.Q_spacelocal*xxs);
%             ENs_local=diag(xxs'*params.N_spacelocal*xxs);
%             
%             EQs = kron(EQs_local,kron(diag(ees), diag(ees))) + kron(ENs_local,kron( diag(ees), ENs_local)) + kron(ENs_local,kron(ENs_local, diag(ees))) ...
%                + kron(ENs_local,kron( diag(ees), ENs_local)) + kron(diag(ees), kron(diag(ees), EQs_local)) + kron(diag(ees),kron(ENs_local,ENs_local)) ...
%                + kron(ENs_local,kron(ENs_local, diag(ees))) + kron(diag(ees),kron(ENs_local, ENs_local)) + kron(diag(ees), kron(EQs_local, diag(ees)));
%             ees=kron(kron(diag(ees),diag(ees)),diag(ees));
%             qq1=kron(kron(xxs,xxs),xxs);
%             %EQs2=diag(qq1'*params.Q_space*qq1);
%             rr1=diag(sqrt(EQs)./sqrt(ees));
%     %       S1=qq1*diag(1./sqrt(ees));
%             sees=1./sqrt(ees);sees=sees(:)';
%             S1=qq1.*repmat(sees,size(qq1,1),1);
%         end
%      else
%          [xxs,ees]=eig(full(params.M_space)); 
%          % Davide: params.Q_space and params.M_space can be simultaneously
%          % diagonalized, we can compute only one eigendecomposition
%          %EQs=diag(xxs'*params.Q_space*xxs);
%          EQs=sum(xxs.*(params.Q_space*xxs))';
%          qq1=xxs;
%          rr1=spdiags(sqrt(EQs)./sqrt(diag(ees)),0,n1,n1);
%          %S1=qq1*diag(1./sqrt(diag(ees)));
%          S1=qq1.*repmat(1./sqrt(diag(ees))',n1,1);
% 
%      end

     %LQt=chol(params.Q_time,'lower'); 
     %LMt=chol(params.M_time,'lower');  
     %[qq2,rr2]=schur(full(LMt\LQt));
     %[qq2,rr2]=schur(full(params.M_time\params.Q_time));
     LQt=chol(params.Q_time,'lower');
     Lt=sqrtm(full(params.M_time));
     hatD=Lt*LQt';
     BB1=params.M_time\hatD;
     [S2,rr4]=schur(full(BB1));
     BB2=params.M_time\hatD';
     [S4,rr2]=schur(full(BB2));
      
     
     
     %rr3=rr1'; %qq4=qq2; 
     %rr4=rr2';
     %S2=qq2';
     %S3=S1';
     %S4=S2';
     %II=eye(size(rr4));

else
% Iterative preconditioner
 LQt=chol(params.Q_time,'lower');
 Lt=sqrtm(full(params.M_time));
 hatD=Lt*LQt';

 hatA=params.A_space;
 BB1=params.M_time\hatD;
 BB2=params.M_time\hatD';
%AA2=hatA/params.M_space; AA1=AA2;
 
 eAmin=abs(eigs(hatA',params.M_space,1,'sm','Tolerance',1e-1)); eAmax=abs(eigs(hatA',params.M_space,1,'lm','Tolerance',1e-1));
 disp([eAmin,eAmax])
 LM=chol(params.M_space,'lower');

end
   
while (backward_error > tol & k<maxit)

  if (flagexact)
% Exact preconditioner (schur decomposition)
    %Z = params.M_space\(S1*lyap(rr1,rr2,-lyap(rr1,rr4,-S1'*R*S2))*S4)/params.M_time';
   % Z = params.M_space\(lyap(params.M_space\params.A_space,BB2,...
   %    -lyap(params.A_space'/params.M_space',BB1,-R)))/params.M_time';
    %Z = params.M_space\(qq1*lyap(rr1,BB2,-lyap(rr1',BB1,-qq1'*R)))/params.M_time';
    Z = params.M_space\(qq1*lyap(rr1,rr2,-lyap(rr1',rr4,-qq1'*R*S2)*S2'*S4)*S4')/params.M_time';
    
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

% Inexact preconditioner (one-sided rational krylov iteration)

    nr=4; droptol=1e-5; tol_inner=1e-10;
    [uu,ss,vv]=svds(R,min([nr,n1])); is = sum(diag(ss)/ss(1,1)>droptol);
    R1=uu(:,1:is)*sqrt(ss(1:is,1:is)); R2=vv(:,1:is)*sqrt(ss(1:is,1:is));
    % solve -AZ-ZB+R=0
    [Z1,Z2,~,~]=arn_adapt_SylvMRHS_oneside(-hatA,params.M_space,LM,-BB1,R1,R2,2*n1,tol_inner,eAmin,eAmax,0,1);
% the svd below works well. the double svd on Z1, Z2 did not work as well...
    [uu,ss,vv]=svds(Z1*Z2',min([nr,size(Z1,2)])); is=sum(diag(ss)/ss(1,1)>droptol); 
    Z1=uu(:,1:is)*sqrt(ss(1:is,1:is)); Z2=vv(:,1:is)*sqrt(ss(1:is,1:is));
    % solve -AZ-ZB+R=0
    [Z1,Z2,~,~]=arn_adapt_SylvMRHS_oneside(-hatA,params.M_space,LM,-BB2,Z1,Z2,2*n1,tol_inner,eAmin,eAmax,0,1);
     Z=(LM'\(LM\Z1))*(Z2'/params.M_time);
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
  backward_error=res(k)/(res0+norm(X,1)*normterm); 
   disp([k,res(k),backward_error]) %JH

end	
 fprintf('pcg_fun3: Dim %d %d .   CG its %d  backward error %d\n',n1,n2,k,backward_error)

