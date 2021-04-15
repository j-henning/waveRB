function [G,D,H,err] = trunc(G,D,H,trunctol,opts)
% truncation of a factored matrix

% G0,D0,H0 - factors of matrix G0*D0*H0'
% trunctol - truncation tolerance: all svals & svecs with 
%           sigma/sigma_max<trunctol are neglected
% opts      - various options

% G,D,H - factors after truncation
if ~exist('opts','var'),
    opts.issym=0;
    opts.method='qrsvd';
    opts.maxcol=size(G,2);
end
if ~exist('H0', 'var') || ~isfield(opts,'issym'),
    issym = false;
else
    issym=opts.issym;
end
if ~exist('D','var'), D=1; end


switch(opts.method)
    case 'qrsvd'
        [G,Sigma] = qr(full(G),0);
        if ~issym
            [H,Sigma2] = qr(full(H),0);  
            [Q,s,Q2] = svd(full(Sigma*(D*Sigma2')));
           
%             s = diag(s);
        else
            [Q,s]=eig(Sigma*D*Sigma'); 
            [~, ord] = sort(diag(abs(s)), 'descend');
            s = s(ord, ord); Q = Q(:, ord);
        end
        maxsig = max(abs(diag(s)));
        rk = sum(diag(abs(s)) >= trunctol*maxsig);
        rk = min([rk,size(Q,2),size(Q2,2)]); %safe guard
       
        if rk<size(s,2); err = sum(abs(diag(s(rk+1,rk+1))));else, err=0; end

        if ~issym
            G = G*Q(:,1:rk)*diag(sqrt(diag(s(1:rk,1:rk))));
            H = H*Q2(:,1:rk)*diag(sqrt(diag(s(1:rk,1:rk))));
            D = eye(rk);    
        else
            G = G*Q(:,1:rk); D=s(1:rk,1:rk);
        end
    case 'itersvd' %todo
    case 'rsvd' %randQR-svd type compression
        if ~isfield(opts,'ell'), opts.ell=min([size(G,2),100]); end %rank estimate,limits
        if ~isfield(opts,'k'), opts.k=5; end %oversampling para
        if ~isfield(opts,'subit'), opts.subit=2; end %subspace iters
        [G,Sigma]=randQR(G,opts.ell,opts.k,opts.subit);
        if ~issym
            [H,Sigma2]=randQR(H,opts.ell,opts.k,opts.subit);
            [Q,s,Q2] = svd(full(Sigma*(D*Sigma2')),0);
        else
            [Q,s]=eig(Sigma*D*Sigma'); 
            [~, ord] = sort(diag(abs(s)), 'descend');
            s = s(ord, ord); Q = Q(:, ord);
        end
         maxsig = max(abs(diag(s)));
        rk = sum(diag(abs(s)) >= trunctol*maxsig);
        rk = min([rk,size(Q,2),size(Q2,2)]); %safe guard
        if rk<size(s,2)
%             pause
%             [size(s),rk]
            err = sum(abs(diag(s(rk+1,rk+1))));
        else, 
            err=0; 
        end

        if ~issym
            G = G*Q(:,1:rk)*diag(sqrt(diag(s(1:rk,1:rk))));
            H = H*Q2(:,1:rk)*diag(sqrt(diag(s(1:rk,1:rk))));
            D = eye(rk);    
        else
            G = G*Q(:,1:rk); D=s(1:rk,1:rk);
        end
%     case 'rsvd2'
%         if ~isfield(opts,'ell'), opts.ell=100; end %rank estimate,limits
%         if ~isfield(opts,'k'), opts.k=5; end %oversampling para
%         if ~isfield(opts,'subit'), opts.subit=2; end %subspace iters
%         tmp=randn(size(H,1),min(size(H,1),opts.ell+opts.k));
%         
%         [Q,~] = qr(G*D*(H'*tmp),0);
%         for j=1:opts.subit,
%             [Q,~] = qr(H*D'*(G'*Q),0);
%             [Q,~] = qr(G*D*(H'*Q),0);
%         end
%         Gn=Q;
%         [t1,r1]=qr(Q'*G,0); [t2,r2]=qr(H,0);
%         [Q,s,Q2]=svd(r1*D*r2',0);
% %         M=Q*UG;
% %         end
% %         if ~issym
% %             [H,Sigma2]=randQR(H,opts.ell,opts.k,opts.subit);
% %             [Q,s,Q2] = svd(full(Sigma*(D*Sigma2')),0);
% %         else
% %             [Q,s]=eig(Sigma*D*Sigma'); 
% %             [~, ord] = sort(diag(abs(s)), 'descend');
% %             s = s(ord, ord); Q = Q(:, ord);
% %         end
%          maxsig = max(abs(diag(s)));
%         rk = sum(diag(abs(s)) > trunctol*maxsig);
%         rk = min([rk,size(Q,2),size(Q2,2)]); %safe guard
%         if rk<size(s,2)
% %             pause
% %             [size(s),rk]
%             err = sum(abs(diag(s(rk+1,rk+1))));
%         else, 
%             err=0; 
%         end
% 
%         if ~issym
%             G = Gn*t1*Q(:,1:rk)*diag(sqrt(diag(s(1:rk,1:rk))));
%             H = H*t2*Q2(:,1:rk)*diag(sqrt(diag(s(1:rk,1:rk))));
%             D = eye(rk);    
%         else
%             H = G*t1*Q(:,1:rk); D=s(1:rk,1:rk);
%         end
%     otherwise
end
end

function [M,s]=randQR(M,ell,k,subit)
        [n1,m1]=size(M);
        tmp=randn(m1,min(m1,ell+k));        
        [Q,~] = qr(M*tmp,0);
        for j=1:subit,
            [Q,~] = qr(M'*Q,0);
            [Q,~] = qr(M*Q,0);
        end
        [UG,s]=qr(Q'*M,0);
        M=Q*UG;
end
