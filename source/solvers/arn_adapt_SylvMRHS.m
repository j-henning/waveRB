function [Z1,Z2,nrmres,sdim]=arn_adapt_SylvMRHS(A,M,B,C1,C2,m,tol,s1,smax,s1B,smaxB,ch,period);
%function [Z1,Z2,nrmres,sdim]=arn_adapt_SylvMRHS(A,M,B,C1,C2,m,tol,s1,smax,s1B,smaxB,ch,period);
%      
% Approximately Solve  
%                A X  +  X B = C1 C2' 
%
% by the Rational Krylov subspace method 
% (Galerkin condition onto the Rational Krylov subspace)
%
% Input:  TO BE REVISED
%
% A, B  coeff. matrices.   BOTH  Stable matrices  (spec(A),spec(B) in C^-)
% C1,C2     rhs factors
% m       max space dimension allowed
% tol     stopping tolerance (backward error)
% s1,smax estimates for real spectral interval of A (associated with F(A))
% s1B,smaxB estimates for real spectral interval of B (associated with F(B))
% ch      ch=1  complex poles  ch=0 real poles
% period  how often the projected problem is solved
%
% Output:
%
% Z1, Z2    factors of approximate solution  X = Z1 Z2'
% nrmres    final residual norm
% sdim      max space dim
%
% Hints:
% 1) Before the call permute entries of A, B and C1,C2 so as to 
%    limit fill-in in the system solves
% 2) Provide "comfortable" (loose bounds) estimates s1, emax, s1B, emaxB


ttime=cputime;
tet = 10.; hinorm=0.;
% Done outside
% perm=symamd(A); A=A(perm,perm);E=E(perm,perm); B=B(perm,:);
%eA=eig(full(A));
[n,n]=size(A);
[n,p]=size(C1);
[nB,nB]=size(B);
[nB,pB]=size(C2);
I=speye(p);O=0*I;
In=speye(n);
InB=speye(nB);
Lres=C1;
LresB=C2;
[V,rr]=qr(Lres,0); nrmb=norm(rr,'fro');
beta=V'*Lres; 
[W,rrB]=qr(LresB,0); nrmb=nrmb*norm(rrB,'fro');
betaB=W'*LresB; 
beta2=beta*betaB';
%nrmb0=norm(B,'fro')^2;

VV=V;
WW=W;
H=zeros(p*(2*m),p*(2*m));
T=zeros(p*(2*m),p*(2*m));
nrmrestot=[];
nrma=norm(A,'fro');
%nrmb=norm(B,'fro');

%if (norm(A-E-(A-E)',1)<1e-14), symm=1; else symm=0;end
i=0;
iB=0;

% sequence in A
 newAv=A*(M\V);
 K=full(V'*newAv);   % K=V'*(L\(A*(L'\V)));
 s=s1;
 eH=eig(K);
 snew=select_new_pole(eH,[s,s],s1,smax,ch,p,i);
%eHpoints = sort([s1,emax]);
%snew=newpolei(eHpoints,eH,s1*ones(p,1));
 if real(snew)<0, snew=-real(snew)+sqrt(-1)*imag(snew);end
 s=[s,snew];

% sequence in B
 newBw=B'*W;
 KT=full(W'*newBw);   % K=V'*(L\(A*(L'\V)));
 sB=s1B;
 eH=eig(KT);
 snewB=select_new_pole(eH,[sB,sB],s1B,smaxB,ch,p,iB);
%eHpoints = sort([s1B,emaxB]);
%snewB=newpolei(eHpoints,eH,s1B*ones(p,1));
 if real(snewB)<0, snewB=-real(snewB)+sqrt(-1)*imag(snewB);end
 sB=[sB,snewB];


% additional steps
cmplxflag=0;
itsinner=0;

cmplxflagB=0;
itsinnerB=0;

while i < m

  i=i+1;
  iB=iB+1;

  paired=0;
  itp=1;
if (size(VV,2)<n)
%Sequence in A
  while (paired==0),
    i1=i+1; it = 0; t=0.;

    wrk = M*((A-snew*M)\V); 
% Gram-Schmidt step
    jms=(i-1)*p+1;j1s=(i+1)*p;js=i*p;js1=js+1;
    for it=1:2,
      for kk=1:i
        k1=(kk-1)*p+1; k2=kk*p; 
        gamma=VV(1:n,k1:k2)'*wrk;
        H(k1:k2,jms:js) = H(k1:k2,jms:js)+ gamma;
        wrk = wrk - VV(:,k1:k2)*gamma;
      end
    end
    [V,H(js1:j1s,jms:js)]=qr(wrk,0);
    if (cmplxflag), snew=conj(snew); s=[s,snew];cmplxflag=0; %newAv=EL\(A*(EU\V(:,i1))); 
       newAv=A*(M\V);
       D = kron(spdiag(s(2:end)),I);
       g = VV'*newAv; g1 = g; g2 = V'*(A*(M\VV)); g3 = V'*newAv; %(A*(M\V);
       K = [K g1; g2, g3];
       VV=[VV,V];
       i=i+1; itp=itp+1;
    else, 
       paired=1; 
    end
  end

    ih1=i1; ih=i;
    newAv=A*(M\V);
    D = kron(spdiag(s(2:end)),I);
    g = VV'*newAv;
end

%K=VV'*A*VV;
%  if (symm), K=(K+K')/2; end

%Sequence in B
if (size(WW,2)<nB)
  pairedB=0;
  while (pairedB==0),

    iB1=iB+1; it = 0; t=0.;
    wrk = (B'-snewB*InB)\W;
% Gram-Schmidt step
    jms=(iB-1)*p+1;j1s=(iB+1)*p;jsB=iB*p;js1=jsB+1;
    for it=1:2,
      for kk=1:iB
        k1=(kk-1)*p+1; k2=kk*p;
        gamma=WW(1:nB,k1:k2)'*wrk;
        T(k1:k2,jms:jsB) = T(k1:k2,jms:jsB)+ gamma;
        wrk = wrk - WW(1:nB,k1:k2)*gamma;
      end
    end
    [W,T(js1:j1s,jms:jsB)]=qr(wrk,0);
    if (cmplxflagB), snewB=conj(snewB); sB=[sB,snewB];cmplxflagB=0; %newAv=EL\(A*(EU\V(:,i1))); 
    nw1=size(WW,2); nw2=size(W,2);
    if nw1+nw2>nB, nw2=nB-nw1;end
    newBw=B'*W(:,1:nw2);
    D = kron(spdiag(sB(2:end)),I);
    gB = WW'*newBw; g1 = gB; g2 = W'*(B'*WW); g3 = W'*(B'*W);
    KT = [KT g1; g2, g3];
    WW=[WW,W];
    iB=iB+1; itp=itp+1;
    else, 
     pairedB=1; 
    end
  end
    ih1=i1; ih=i;
    newBw=B'*W;
    D = kron(spdiag(sB(2:end)),I);
    gB = WW'*newBw;

end

if (rem(i,period)==0 | i >= m | max(js,jsB) >= min([n,nB]) )
     rhs2=speye(js,p)*beta2*speye(jsB,p)';
     Y = lyap(full(K(1:js,1:js)),full(KT(1:jsB,1:jsB)'),rhs2);
%    nrmx = norm(Y,'fro');

% computed residual   (exact, in exact arithmetic)
%    u1=newAv-VV*g;
%    d=-VV*(Y*(H(1:ih*p,1:ih*p)'\[sparse(p*(ih-1),p);I])*H(p*ih+1:p*ih1,p*ih-p+1:p*ih)');
%    U=[-V*s(end),  d u1 ];
%    rr=qr(full(U),0); rr=triu(rr(1:size(rr,2),:));
%    nrmres=norm(rr*sparse([O I O; I O I; O I O ])*rr','fro')/nrmb;

% This can be made cheaper!
     [qq,rr]=qr([VV(:,1:js),A*(M\VV(:,1:js))],0);
     [qq,rrB]=qr([WW(:,1:jsB),B'*WW(:,1:jsB)],0);
     nrmres=norm(rr*sparse([rhs2 Y; Y sparse(js,jsB) ])*rrB','fro')/nrmb;

%backward error
%   %nrmres=norm(rr*sparse([O I O; I O I; O I O ])*rr','fro')/(nrmb+singE*nrma*nrmx);
%    nrmrestot=[nrmrestot,nrmres];

% True residual
 %X=VV(:,1:js)*Y*WW(:,1:jsB)';
 %nrmres=norm(A*X+X*B+C1*C2','fro')/nrmb;

%     disp([size(VV),size(WW),i,nrmres]) 
     if (nrmres<tol | max(js,jsB) >= n), sdim=max(size(VV,2),size(WW,2));break,end

end

% New poles and zeros  in A
     eH=sort(eig(K));
     snew=select_new_pole(eH,s,s1,smax,ch,p,i);
     if real(snew)<0, snew=-real(snew)+sqrt(-1)*imag(snew);end
% If pole is complex, include its conjugate
     if (imag(snew) ~=0), cmplxflag=1;end
     s=[s,snew];
     g1 = g; g2 = V'*(A*(M\VV)); g3 = V'*(A*(M\V));
     K = [K g1; g2, g3];
     VV=[VV,V];

% New poles and zeros  in B'
     eH=sort(eig(KT));
     snewB=select_new_pole(eH,sB,s1B,smaxB,ch,p,iB);
     if real(snewB)<0, snewB=-real(snewB)+sqrt(-1)*imag(snewB);end
% If pole is complex, include its conjugate
     if (imag(snewB) ~=0), cmplxflagB=1;end
if size(KT,2)<nB,
     sB=[sB,snewB];
     g1 = gB; g2 = W'*(B'*WW); g3 = W'*(B'*W);
     KT = [KT g1; g2, g3];
     WW=[WW,W];
end

end;
fprintf('inner sylv solver %d  %d \n',i,nrmres)

% factored solution 
[uY,sY,vY]=svd(Y); 
is=sum(diag(sY)>1e-14);
Y1 = uY(:,1:is)*(sqrt(sY(1:is,1:is))); 
Y2 = vY(:,1:is)*(sqrt(sY(1:is,1:is))); 
sdim=max(size(VV,2),size(WW,2));
Z1 = real(VV(:,1:size(Y,1))*Y1); 
Z2 = real(WW(:,1:size(Y,2))*Y2);
%final_rank=is;
%RKStotal_time=cputime-ttime;
%
%semilogy(p:p:p*length(nrmrestot),nrmrestot,'b-')
%xlabel('space dimension')
%
%fprintf('Space dim %d  Solution rank %d   \n',sdim,max(size(Y)));

%X=Z*Z';
%R=A*X+X*A'+B*B';
%nrmx=norm(X,'fro');
%norm(R,'fro')/(nrmb+singE*nrma*nrmx)
%fprintf('Poles computed during the iteration\n')

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function r=ratfun(x,eH,s)
% 
%for j=1:length(x)
%r(j)=abs(prod( (x(j)-s)./(x(j)-eH) ));
%end
%
%return
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function r=ratfuni(x,eH,s)
%
%for j=1:length(x)
%r(j)=1./abs(prod( (x(j)-s)./(x(j)-eH) ));
%end
%
%return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function snew=newpolei(eHpoints,eH,s)
%
%for j=1:length(eHpoints)-1
%%   snew(j) = fminbnd( @(x) ratfuni(x,eH,s), eHpoints(j),eHpoints(j+1),optimset('TolX',1e-3));
%    sval=linspace(eHpoints(j),eHpoints(j+1),200);
%    [sf,jx] = max (abs(ratfun(sval,eH,s)));
%    snew(j)=sval(jx);
%end
%[sn,jx]=max(abs(ratfun(snew,eH,s)));
%snew=snew(jx);
%return

