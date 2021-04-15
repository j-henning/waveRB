function [Z1,Z2,nrmres,sdim]=arn_adapt_SylvMRHS(A1,A2,B,C1,C2,params);
%function [Z1,Z2,nrmres,sdim]=arn_adapt_SylvMRHS(A1,A2,B,C1,C2,params);
%      
% Approximately Solve  
%                A1 X  + A2 X B = C1 C2' 
%
% by the Rational Krylov subspace method 
% (Galerkin condition onto the Rational Krylov subspace)
%
% WARNING:
% 1) this version of the method performs cholesky factorization of 
% permuted version of A2
%
%
m=params.maxit;
tol=params.tol;
s1=params.emin;
smax=params.emax;
ch=params.ch;

ttime=cputime;
tet = 10.; hinorm=0.;
[n,n]=size(A1);
[n,p]=size(C1);
[nB,nB]=size(B);
[nB,pB]=size(C2);
nrmb=sqrt(trace((C2'*C2)*(C1'*C1)));
I=speye(p);O=0*I;
In=speye(n);
pp=symamd(A2);
A2=A2(pp,pp);
L2=chol(A2,'lower');
A1=A1(pp,pp); 
C1=C1(pp,:);

Lres=L2\C1;  %Lres=A2\C1;
[V,rr]=qr(Lres,0); %nrmb=norm(rr,'fro');
beta=V'*Lres; 
beta2=beta*C2';

VV=V;
H=sparse(p*(2*m),p*(2*m));
nrmrestot=[];
nrma=norm(A1,'fro');

i=0;

% sequence in A
 newAv=L2\(A1*(L2'\V));
 K=full(V'*newAv);   % K=V'*(L\(A*(L'\V)));
 s=s1;
 eH=eig(K);
 eHpoints = sort([s1,smax]);
 snew=newpolei(eHpoints,eH,s1*ones(p,1));

 if real(snew)<0, snew=-real(snew)+sqrt(-1)*imag(snew);end
 s=[s,snew];


% additional steps
cmplxflag=0;
itsinner=0;

while i < m

  i=i+1;

  paired=0;
  itp=1;
  while (paired==0),
    i1=i+1; it = 0; t=0.;

%Sequence in A
    wrk = L2'*( (A1-snew*A2)\(L2*V)); 

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
    newAv=L2\(A1*(L2'\V));
    D = kron(spdiag(s(2:end)),I);
    g = VV'*newAv; g1 = g; g2 = V'*(L2\(A1*(L2'\VV))); g3 = V'*(L2\(A1*(L2'\V)));
    K = [K g1; g2, g3];
    VV=[VV,V];
    i=i+1; itp=itp+1;
    else, 
    paired=1; 
    end
  end

    ih1=i1; ih=i;
    newAv=L2\(A1*(L2'\V));
    D = kron(diag(s(2:end)),I);
    g = VV'*newAv;


     rhs2=speye(js,p)*beta2;  
     Y = lyap(full(K(1:js,1:js)),full(B),rhs2);
%    nrmx = norm(Y,'fro');

% computed residual   (exact, in exact arithmetic)
%    u1=newAv-VV*g;
%    d=-VV*(Y*(H(1:ih*p,1:ih*p)'\[sparse(p*(ih-1),p);I])*H(p*ih+1:p*ih1,p*ih-p+1:p*ih)');
%    U=[-V*s(end),  d u1 ];
%    rr=qr(full(U),0); rr=triu(rr(1:size(rr,2),:));
%    nrmres=norm(rr*sparse([O I O; I O I; O I O ])*rr','fro')

%    [qq,rr]=qr([VV(:,1:js),L2\(A1*VV(:,1:js))],0);
%    nrmres=norm(rr*sparse([rhs2 Y; Y sparse(js,jsB) ])*rrB','fro')/nrmb;

%backward error
%   %nrmres=norm(rr*sparse([O I O; I O I; O I O ])*rr','fro')/(nrmb+singE*nrma*nrmx);
%    nrmrestot=[nrmrestot,nrmres];

     %X=L2'\(VV(:,1:js)*Y);
     %nrmrestrue=norm(A1*X+A2*X*B+C1*C2','fro')/nrmb;
     
     
    % compute the residual 
    gg=V*(H(p*ih+1:p*ih1,1:p*ih)/H(1:ih*p,1:ih*p));
    factor1=snew*gg;
    factor2=L2\(A1*(L2'\gg));
    factor3=VV*(VV'*factor2);
    [~,C]=qr(L2*(factor1-factor2+factor3),0);
    nrmres=norm(triu(C)*Y*params.Mtime,'fro')/nrmb;
    %abs(nrmres-nrmrestrue)/nrmrestrue
     %fprintf('Iteration %d,  relative residual %d\n',i,nrmres)

     %disp([i,nrmres]) 
     if (nrmres<tol || i >= m), sdim=max(size(VV,2),size(VV,2));break,end


% New poles and zeros  in A
     eH=sort(eig(K));
  eHorig=eH;

  if (ch)                     % Complex poles. Compute set for next complex pole of r_m

    if (any(imag(eH)) ~=0 & max(abs(imag(eH)))>1e-5 & length(eH)>2) % Roots lambdas come from convex hull too
     eH=[eH;-smax];
      ij=convhull(real(eH),imag(eH)); eH=eH(ij);
      ieH=length(eH); missing=ih*p-ieH;
      while missing>0,                         % include enough points from the border
        neweH=(eH(1:ieH-1)+eH(2:ieH))/2;missing=ih*p-length(eH);
        eH=[eH;neweH];
      end
    % eH=eH(1:ih);
      eHpoints=-eH;
      eH=eHorig;
    else                                  % if all real eigs, no convex hull possible
      eHpoints = sort([s1; smax.';-real(eH)]);
    end


  else   % Real poles s from real set. Compute complex roots of r_m via Ritz convex hull

     if (any(imag(eH)) ~=0 & length(eH)>2)    % Roots lambdas come from convex hull too
       eH=[eH;-s1;-smax.'];
       ij=convhull(real(eH),imag(eH)); eH=eH(ij);
       ieH=length(eH); missing=ih*p-ieH;
       while missing>0, % include enough points from the border
         neweH=(eH(1:ieH-1)+eH(2:ieH))/2;
         eH=[eH;neweH];
         missing=ih*p-length(eH);
       end
       eH=eH(1:ih*p);
     end
      eHpoints = sort([s1; smax.';-real(eH)]);
      eH=eHorig;
  end

    
  gs=kron(s(2:end),ones(1,p))';
  snew = newpolei(eHpoints,eH,gs);
     
     if real(snew)<0, snew=-real(snew)+sqrt(-1)*imag(snew);end

% If pole is complex, include its conjugate
     if (imag(snew) ~=0), cmplxflag=1;end
     s=[s,snew];
     g1 = g; g2 = V'*(L2\(A1*(L2'\VV))); g3 = V'*(L2\(A1*(L2'\V)));
     K = [K g1; g2, g3];
     
    %  size(VV)
    
    %u1=L2\(A1*(L2'\VV))-VV*g;
    %d=-VV*(Y*(H(1:ih*p,1:ih*p)'\[sparse(p*(ih-1),p);I])*H(p*ih+1:p*ih1,p*ih-p+1:p*ih)');
    %U=L2*[-V*s(end),  d u1 ];    
    %[~,rr]=qr(full(U),0); 
          
    %rr=triu(rr);
     % abs residual
  % norm(rr*sparse([O I O; I O I; O I O ])*rr','fro');
   % [~,C]=qr(L2*(factor1-factor2+factor3),0);
    %abs(nrmres-rr/norm(C1,'fro')/norm(C2,'fro'))/nrmres
    %pause

    VV=[VV,V];
     

end;


sdim=size(VV,2);

%fprintf('\n')
% factored solution 
[uY,sY,vY]=svd(Y);
is=sum(diag(sY)>1e-12);
Z1(pp,:)=L2'\(VV*uY(:,1:is));
Z2 = sY(1:is,1:is)*vY(:,1:is)';
%fprintf('Space dim %d  Solution rank %d   \n',sdim,is);
%fprintf('\n')

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
%
for j=1:length(eHpoints)-1
%   snew(j) = fminbnd( @(x) ratfuni(x,eH,s), eHpoints(j),eHpoints(j+1),optimset('TolX',1e-3));
    sval=linspace(eHpoints(j),eHpoints(j+1),200);
    [sf,jx] = max (abs(ratfun(sval,eH,s)));
    snew(j)=sval(jx);
end
[sn,jx]=max(abs(ratfun(snew,eH,s)));
snew=snew(jx);
return

