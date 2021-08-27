function snew=select_new_pole(eH,s,s1,emax,ch,p,ih);


% New poles and zeros
   % eH=sort(eig(K)); 
    eHorig=eH;
 if (ch)                     % Complex poles. Compute set for next complex pole of r_m
    if (any(imag(eH)) ~=0 & length(eH)>2) % Roots lambdas come from convex hull too
     eH=[eH;-emax];
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
      eHpoints = sort([s1; emax.';-real(eH)]);
    end


 else   % Real poles s from real set. Compute complex roots of r_m via Ritz convex hull
     if (any(imag(eH)) ~=0 & length(eH)>2)    % Roots lambdas come from convex hull too
       eH=[eH;-s1;-emax.'];
       ij=convhull(real(eH),imag(eH)); eH=eH(ij);
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


 gs=kron(s(2:end),ones(1,p))';

 snew = newpolei(eHpoints,eH,gs);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r=ratfun(x,eH,s)

r=zeros(length(x),1);
for j=1:length(x)
r(j)=abs(prod( (x(j)-s)./(x(j)-eH) ));
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function r=ratfuni(x,eH,s)
%
%for j=1:length(x)
%r(j)=1./abs(prod( (x(j)-s)./(x(j)-eH) ));
%end
%
%return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function snew=newpolei(eHpoints,eH,s)

for j=1:length(eHpoints)-1
%   snew(j) = fminbnd( @(x) ratfuni(x,eH,s), eHpoints(j),eHpoints(j+1),optimset('TolX',1e-3));
    sval=linspace(eHpoints(j),eHpoints(j+1),10);
    [sf,jx] = max (abs(ratfun(sval,eH,s)));
    snew(j)=sval(jx);
end
[sn,jx]=max(abs(ratfun(snew,eH,s)));
snew=snew(jx);
return

