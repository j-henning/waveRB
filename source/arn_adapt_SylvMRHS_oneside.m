function [Z1,Z2,nrmres,sdim]=arn_adapt_SylvMRHS_oneside(A1,A2,L2,B,C1,C2,m,tol,s1,smax,ch,period)
%function [Z1,Z2,nrmres,sdim]=arn_adapt_SylvMRHS_oneside(A1,A2,L2,B,C1,C2,m,tol,s1,smax,ch,period);
%      
% Approximately Solve  
%                A1 ((L2 L2')\X)  +  X B = C1 C2' 
%
% by the Rational Krylov subspace method 
% (Galerkin condition onto the Rational Krylov subspace)
%
% Input:  TO BE REVISED
%
% A, E  coeff. matrices.  E is spd
% EL, EU  (lower and upper) Cholesky factors of E
% B     rhs factor
% m       max space dimension allowed
% tol     stopping tolerance (backward error)
% s1,smax estimates for real spectral interval
%         associated with field of values of (A,E)
% ch      ch=1  complex poles  ch=0 real poles
%
% Output:
%
% Z    factor of approximate solution  X = Z Z'
% VV   orthonormal basis for the approx space
% K    VV'*(EL\A/EU) VV
% s    sequence of generated poles
%
% Hints:
% 1) Before the call permute entries of A, E and B so as to 
%    limit fill-in in the system solves
% 2) Provide "comfortable" (loose bounds) estimates s1, smax


ttime=cputime;
tet = 10.; hinorm=0.;
% Done outside
% perm=symamd(A); A=A(perm,perm);E=E(perm,perm); B=B(perm,:);
%eA=eig(full(A));
[n,n]=size(A1);
[n,p]=size(C1);
[nB,nB]=size(B);
[nB,pB]=size(C2);
I=speye(p);O=0*I;
In=speye(n);
%L2=chol(A2,'lower');

Lres=C1;  %Lres=A2\C1;
[V,rr]=qr(Lres,0); nrmb=norm(rr,'fro');
beta=V'*Lres; 
beta2=beta*C2';

VV=V;
H=sparse(p*(2*m),p*(2*m));
nrmrestot=[];
nrma=norm(A1,'fro');

i=0;

% sequence in A
 LV=L2'\(L2\V);
 newAv=A1*LV;   %(L2'\(L2\V));
%newAv=(A1*(A2\V));
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
    wrk = A2*((A1-snew*A2)\V); 

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
    sH=svd(full(H(js1:j1s,jms:js)));
    if sH(1)<1e-14,break,end
    if (cmplxflag), snew=conj(snew); s=[s,snew];cmplxflag=0; %newAv=EL\(A*(EU\V(:,i1))); 
      LVnew=L2'\(L2\V);
      newAv=A1*LVnew;  %(L2'\(L2\V));
   %newAv=A1*(A2\V);
      D = kron(spdiag(s(2:end)),I);
      g = VV'*newAv; g1 = g; g2 = V'*(A1*LV); g3 = V'*newAv; %(A1*(A2\V));
      %g = VV'*newAv; g1 = g; g2 = V'*(A1*(L2'\(L2\VV))); g3 = V'*newAv; %(A1*(A2\V));
      K = [K g1; g2, g3];
      VV=[VV,V];
      LV=[LV,LVnew];
      i=i+1; itp=itp+1;
      else, 
      paired=1; 
      end
  end

    ih1=i1; ih=i;
    %newAv=A1*(A2\V);
    LVnew=L2'\(L2\V);
    newAv=A1*LVnew; %(L2'\(L2\V));
    D = kron(diag(s(2:end)),I);
    g = VV'*newAv;


if (rem(i,period)==0 | i >= m | max(js,js) >= n)
     rhs2=speye(js,p)*beta2;  
     Y = lyap(full(K(1:js,1:js)),full(B),rhs2);
%    nrmx = norm(Y,'fro')

% computed residual   (exact, in exact arithmetic)
%    u1=newAv-VV*g;
%    d=-VV*(Y*(H(1:ih*p,1:ih*p)'\[sparse(p*(ih-1),p);I])*H(p*ih+1:p*ih1,p*ih-p+1:p*ih)');
%    U=[-V*s(end),  d u1 ];
%    rr=qr(full(U),0); rr=triu(rr(1:size(rr,2),:));
%    nrmres=norm(rr*sparse([O I O; I O I; O I O ])*rr','fro')

% This could be made cheaper
     [qq,rr]=qr([VV(:,1:js),A1*LV(:,1:js)],0);
     %[qq,rr]=qr([VV(:,1:js),A1*(L2'\(L2\VV(:,1:js)))],0);
     nrmres=norm(rr*sparse([rhs2+Y*B; Y ]),'fro')/norm(C1,'fro')/norm(C2,'fro');

%backward error
    %nrmres=norm(rr*sparse([O I O; I O I; O I O ])*rr','fro')/(nrmb+singE*nrma*nrmx);
%    nrmrestot=[nrmrestot,nrmres];

    if sH(1)<1e-14,
    %X=(VV(:,1:js)*Y);
    %nrmresx=norm(A1*(L2'\(L2\X))+X*B+C1*C2','fro')/norm(C1,'fro')/norm(C2,'fro')
     i=m; break,
    end
%    fprintf('Iteration %d,  relative residual %d\n',i,nrmres)

%     disp([i,nrmres]) 
     if (nrmres<tol | max(js,js) >= n), sdim=max(size(VV,2),size(VV,2));break,end

end

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
     LVnew=L2'\(L2\V);
     g1 = g; g2 = V'*(A1*LV); g3 = V'*(A1*LVnew);
     %g1 = g; g2 = V'*(A1*(A2\VV)); g3 = V'*(A1*(A2\V));
     K = [K g1; g2, g3];
     VV=[VV,V];
     LV=[LV,LVnew];


end;


% factored solution 
%Y(1:min([5,size(Y,1)]),1:5)
%normY=norm(Y)
[uY,sY,vY]=svd(Y);
is=sum(diag(sY)>1e-14);
Y1 = uY(:,1:is)*((sY(1:is,1:is)));
Y2 = vY(:,1:is);
sdim=size(VV,2);
Z1 = VV(:,1:size(Y1,1))*Y1;
Z2 = Y2;
%X=Z1*Z2';
%     nrmresx=norm(A1*(L2'\(L2\X))+X*B+C1*C2','fro')/norm(C1,'fro')/norm(C2,'fro')

%fprintf('Inner its %d,  rel. res. %d  Space dim %d  Soln rank %d   \n',i, nrmres,size(VV,2),size(Y1,2))

%X=Z*Z';
%R=A*X+X*A'+B*B';
%nrmx=norm(X,'fro');
%norm(R,'fro')/(nrmb+singE*nrma*nrmx)
%fprintf('Poles computed during the iteration\n')

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

