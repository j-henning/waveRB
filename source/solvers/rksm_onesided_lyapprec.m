% ???
% NOTE:
% This code was provided by Valeria Simoncini 
% http://www.dm.unibo.it/~simoncin/welcome.html
function [Z,nrmrestot]=rksm_onesided_lyapprec(A,E,EL,B,C1,C2,m,s1,emax,ch,tol)
%function [Z,nrmrestot]=rksm(A,E,EL,B,m,tol,s1,emax,ch,tolY)
%      
% Approximately Solve  
%                A X E + E X A' + BB' = 0
%
% by the Rational Krylov subspace method 
% (Galerkin condition onto the Rational Krylov subspace)
% This code performs system solves with (A-s E)
%
% Author of the matlab code:  Valeria Simoncini
% version 1.0
%  The code uses the function lyap.m of the Control Matlab Toolbox
%
%
% Input:
%
% A, E  coeff. matrices. A<0,  E is spd
% EL   Cholesky lower factor of E
% B     rhs factor
% m       max space dimension allowed
% tol     stopping tolerance, with stopping criterion
%          ||LE\A X LE  + LE' X A'/LE'-LE\BB'/LE'||
%          ----------------------------------------  < tol
%      ||LE\BB'/LE'|| + ||E^{-1}|| ||A|| ||LE'X LE ||
%         computed in a cheap manner
% s1,smax estimates for real spectral interval
%         associated with field of values of (A,E)
% ch      ch=1  complex poles  ch=0 real poles
% tolY    truncation tolerance for final solution, e.g., tolY=1e-12
%
% Output:
%
% Z    factor of approximate solution  X = Z Z'
%      If ch=1,  Z may be complex
% nrmrestot  history of residuals
%
% Hints:
% 1) Before the call permute entries of A, E and B so as to 
%    limit fill-in in the system solves
% 2) Provide "comfortable" (loose bounds) estimates s1, emax
%
%
%  Please contact V. Simoncini for any problem you may encouter when
%  running the code
%
%  When using this code, please cite the following reference:
%
%  V. Druskin and V. Simoncini, 
%  Adaptive rational Krylov subspaces for large-scale dynamical systems
%  Systems & Control Letters, 60 (2011), pp. 546-560.
%
%
%THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
%IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS 
%FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR 
%COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER 
%IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
%CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
%  


tic;
[n,~]=size(A);
B=full(B);
p=size(C1,2);
I=speye(p);O=0*I;
uno=ones(1,p);
Lres=EL\C1;
[V,irr]=qr(Lres,0);rr=inv(irr);
nrmb=norm(inv(rr),'fro')^2; beta=V'*Lres; beta2=beta*C2';
s=zeros(m+2,1);

%fprintf('     no its     backward error\n')
VV=zeros(n,p*(m+2));
VV(1:n,1:p)=V;
H=zeros(p*(m+2),p*(m+1));
nrmrestot=[];
nrma=norm(A,'fro');
%if (norm(E-speye(n),'fro')>1e-14)
%  condestE=condest(E);
%  singE=condestE/norm(E,'fro');
%else
%  singE=1;
%end

if (norm(A-E-(A-E)',1)<1e-14), symm=1; else symm=0;end

 newAv=EL\(A*(EL'\V));
    Q=newAv;
 K=V'*newAv;
 s(1)=s1;
 eH=eig(K);
 eHpoints = sort([s1,emax]);
 snew=newpolei(eHpoints,eH,s1*uno');
 if real(snew)<0, snew=-real(snew)+sqrt(-1)*imag(snew);end
 s(2)=snew;

% additional steps
cmplxflag=0;

i=0;

while i < m

  i=i+1;

  paired=0;
  while (paired==0),

    i1=i+1; 
    w=EL*V;
    wrk = (A-snew*E)\w; 
    wrk= EL'*wrk;  

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
    if (cmplxflag), 
    snew=conj(snew); s(i+2)=snew;cmplxflag=0;
    newAv=EL\(A*(EL'\V));
    Q=[Q, newAv];
    g = VV(1:n,1:js)'*newAv;
    g1 = g; 
    Vwrk=newAv; %EL\(A*(EL'\V));
    g2 = Vwrk'*VV(1:n,1:js);
   %g2 = V'*(EL\(A*(EL'\VV(1:n,1:js))));
    g3 = V'*Vwrk;
   %g3 = V'*(EL\(A*(EL'\V)));
    K = [K g1; g2, g3];
    VV(1:n,js+1:j1s)=V;  
    i=i+1;
    else, paired=1; end
  end

    ih1=i1; ih=i;
    newAv=EL\(A*(EL'\V));
    Q=[Q, newAv];
    g = VV(1:n,1:js)'*newAv;

    if (symm), K=(K+K')/2; end
   %rhs2=speye(ih*p,p)*beta2*speye(ih*p,p)';
    rhs2=speye(ih*p,p)*beta2;
    Y = lyap(K,B,rhs2);
jY=size(Y,1);
%X=VV(:,1:size(Y,1))*Y;
%WW=EL\(A*(EL'\X))+X*B+Lres*C2';
%nrmres=norm(WW,'fro')/nrmb;
%      [qq,rr]=qr([EL\(A*(EL'\VV(:,1:jY))),VV(:,1:jY), Lres],0);
       [qq,rr]=qr([Q(:,1:jY),VV(:,1:jY), Lres],0);
        %[qq,rr]=qr([VV(:,1:js),A1*(L2'\(L2\VV(:,1:js)))],0);
        nrmres=norm(rr*sparse([Y; Y*B; C2' ]),'fro')/norm(Lres,'fro')/norm(C2,'fro');
%disp([i,nrmres])

if (nrmres<tol), break,end
%      nrmx = norm(Y,'fro');
% 
% % computed residual   (exact, in exact arithmetic)
%      u1=newAv-VV(1:n,1:js)*g;
%      d=-VV(1:n,1:js)*(Y*(H(1:ih*p,1:ih*p)'\[sparse(p*(ih-1),p);I])*H(p*ih+1:p*ih1,p*ih-p+1:p*ih)');
%      U=[-V*s(end),  d u1 ];
%      rr=qr(full(U),0); rr=triu(rr(1:size(rr,2),:));
%      nrmres=norm(rr*sparse([O I O; I O I; O I O ])*rr','fro')/(nrmb+singE*nrma*nrmx);
%      nrmrestot=[nrmrestot,nrmres];
  
%      disp([i,nrmres]) 
%     if (nrmres<tol), break,end

% New poles and zeros
    eH=sort(eig(K)); eHorig=eH;

 if (ch)                     % Complex poles. Compute set for next complex pole of r_m

    if (any(imag(eH)) ~=0 && length(eH)>2) % Roots lambdas come from convex hull too
     eH=full([eH;-emax]);
      ij=convhull(real(eH),imag(eH)); eH=eH(ij);
      ieH=length(eH); missing=ih*p-ieH;
      while missing>0,                         % include enough points from the border
        neweH=(eH(1:ieH-1)+eH(2:ieH))/2;missing=ih*p-length(eH);
        eH=[eH;neweH];
      end
      eHpoints=-eH;
      eH=eHorig;
    else                                  % if all real eigs, no convex hull possible
      eHpoints = sort([s1; emax.';-real(eH)]);
    end


 else   % Real poles s from real set. Compute complex roots of r_m via Ritz convex hull
     if (any(imag(eH)) ~=0 && length(eH)>2)    % Roots lambdas come from convex hull too
       eH=full([eH;-s1;-emax.']);
       ij=convhull(real(eH),imag(eH)); eH=eH(ij);
      %ij=convhull(real(eH),imag(eH),{'Qt','Pp'}); eH=eH(ij);
       ieH=length(eH); missing=ih*p-ieH;
       while missing>0, % include enough points from the border
         neweH=(eH(1:ieH-1)+eH(2:ieH))/2;
         eH=[eH;neweH];
         missing=ih*p-length(eH);
       end
       eH=eH(1:ih*p);
     end
      eHpoints = sort([s1; emax.';-real(eH)]);
      eH=eHorig;
 end


 gs=kron(s(2:i+1).',uno)'; 
 snew = newpolei(eHpoints,eH,gs);
 if real(snew)<0, snew=-real(snew)+sqrt(-1)*imag(snew);end  %safeguard strategy

% If pole is complex, include its conjugate
 if (imag(snew) ~=0), cmplxflag=1;end
 s(i+2)=snew; 

   %newAv=EL\(A*(EL'\V));
   %g = VV(1:n,1:js)'*newAv;
 g1 = g; 
    Vwrk=newAv;  %EL\(A*(EL'\V));
 %g2 = V'*(EL\(A*(EL'\VV(1:n,1:js))));
    g2 = g1'; % Vwrk'*VV(1:n,1:js);
 %g3 = V'*(EL\(A*(EL'\V)));
    g3 = V'*Vwrk;
 K = [K g1; g2, g3];
 VV(1:n,js+1:j1s)=V;  

end;

% Done
% Reduce rank of solution, if needed
% [uY,sY]=eig(Y); [sY,id]=sort(diag(sY));
% sY=flipud(sY); uY=uY(:,id(end:-1:1));
% is=sum(abs(sY)>tolY);
% Y0 = uY(:,1:is)*diag(sqrt(sY(1:is))); 
Z = EL'\(VV(:,1:size(Y,1))*Y); 
RKStotal_time=toc;


%fprintf('Space dim %d  Solution rank %d  time %d \n',j1s,is,RKStotal_time);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary routines

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r=ratfun(x,eH,s)
r = abs(prod((x-s) ./(x-eH), 1));
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function snew=newpolei(eHpoints,eH,s)

for j=1:length(eHpoints)-1
    sval=linspace(eHpoints(j),eHpoints(j+1),200);
    [sf,jx] = max (abs(ratfun(sval,eH,s)));
    snew(j)=sval(jx);
end
[sn,jx]=max(abs(ratfun(snew,eH,s)));
snew=snew(jx);
return
