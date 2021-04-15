function [X, iter]=pcg_fun4(funA,b,x,params,maxit,tol,n1,n2,flagexact)
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

if ~isfield(params,'info')
    params.info=0;
end

b=reshape(b,n1,n2);
X=reshape(x,n1,n2);
R = b-funA(X);  %R=R - alfa*AP;
%R=b; %AP=X;
gamma=sum(sum(R.*R));
res0=sqrt(gamma);
res=full(res0);
k=0;
normterm=norm(params.M_time,1)*norm(params.Q_space,1)+norm(params.M_space,1)*norm(params.Q_time,1);
backward_error=1;

if strcmp(params.precond,'lyap')
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
    end
    
elseif strcmp(params.precond,'optimal')

    if (flagexact)
 
        % Direct preconditioner
        % Compute the Schur decomposition once and for all
        [qq1,rr1]=schur(full(params.M_space\params.A_space));
        LQt=chol(params.Q_time,'lower');
        Lt=sqrtm(full(params.M_time));
        hatD=Lt*LQt';
        BB1=params.M_time\hatD;
        [S2,rr4]=schur(full(BB1));
        BB2=params.M_time\hatD';
        [S4,rr2]=schur(full(BB2));
        
    else
        
        % Iterative preconditioner
        LQt=chol(params.Q_time,'lower');
        Lt=sqrtm(full(params.M_time));
        hatD=Lt*LQt';

        hatA=params.A_space;
        BB1=params.M_time\hatD;
        BB2=params.M_time\hatD';
 
        % prepare all we need for the iterative solution of the
        % preconditioning step
        eAmin=abs(eigs(hatA',params.M_space,1,'sm','Tolerance',1e-1)); eAmax=abs(eigs(hatA',params.M_space,1,'lm','Tolerance',1e-1));
        LM=chol(params.M_space,'lower');
    end
end
   
% start the actual loop
while (backward_error > tol && k<maxit)
    
    % preconditioning step
    if strcmp(params.precond,'lyap')
        if (flagexact)
            Z=U1*( (V1*R*V2).*LL)*U2;        
        else
            % low-rank truncation of the rhs 
            nr=4; droptol=1e-5; tol_inner=1e-10;
            [uu,ss,vv]=svds(R,min([nr,n1])); is = sum(diag(ss)/ss(1,1)>droptol);
            R1=uu(:,1:is)*sqrt(ss(1:is,1:is)); R2=vv(:,1:is)*sqrt(ss(1:is,1:is)); 
            [Z1,Z2,~,~]=arn_adapt_SylvMRHS_oneside(-AA,speye(n1),speye(n1),-BB,R1,R2,n1,tol_inner,eAmin,eAmax,0,1);
            
            Z=(Ls'\Z1)*(Z2'/Lt);
        end
    
    elseif strcmp(params.precond,'optimal')
        
        if (flagexact)
      
            % Exact preconditioner (schur decomposition)
            Z = params.M_space\(qq1*lyap(rr1,rr2,-lyap(rr1',rr4,-qq1'*R*S2)*S2'*S4)*S4')/params.M_time';
        else
            
            % Inexact preconditioner (one-sided rational krylov iteration)
            
            % low-rank truncation of the rhs of the firs Sylvetser eq.
            nr=4; droptol=1e-5; tol_inner=1e-10;
            [uu,ss,vv]=svds(R,min([nr,n1])); is = sum(diag(ss)/ss(1,1)>droptol);
            R1=uu(:,1:is)*sqrt(ss(1:is,1:is)); R2=vv(:,1:is)*sqrt(ss(1:is,1:is));
            % solve the first equation
            [Z1,Z2,~,~]=arn_adapt_SylvMRHS_oneside(-hatA,params.M_space,LM,-BB1,R1,R2,2*n1,tol_inner,eAmin,eAmax,0,1);
            % Further reduce the rank of Z1*Z2' to make the solution of the
            % second equation cheaper
            [uu,ss,vv]=svds(Z1*Z2',min([nr,size(Z1,2)])); is=sum(diag(ss)/ss(1,1)>droptol); 
            Z1=uu(:,1:is)*sqrt(ss(1:is,1:is)); Z2=vv(:,1:is)*sqrt(ss(1:is,1:is));
            % solve the second equation
            [Z1,Z2,~,~]=arn_adapt_SylvMRHS_oneside(-hatA,params.M_space,LM,-BB2,Z1,Z2,2*n1,tol_inner,eAmin,eAmax,0,1);
            Z=(LM'\(LM\Z1))*(Z2'/params.M_time);
        end
    end

  k=k+1;
  gamma=sum(sum(R.*Z));
  if (k==1)
      P=Z;
  else
      beta=gamma/gamma0;
      P=Z+beta*P;
  end
  AP=funA(P);
  delta=sum(sum(P.*AP));
  alfa = gamma/delta;
  X = X + alfa*P;
  R = b-funA(X);  %R=R - alfa*AP;
  gamma0=gamma;
  res(k)=full(sqrt(sum(sum(R.*R))));
  backward_error=res(k)/(res0+norm(X,1)*normterm); 
  if params.info
      disp([k,res(k),backward_error]) %JH
  end

end	
iter = k;
disp(['pcg (',params.precond,')'])
fprintf('\nDim space: %d time: %d .  CG its %d  backward error %d\n',n1,n2,k,backward_error)

