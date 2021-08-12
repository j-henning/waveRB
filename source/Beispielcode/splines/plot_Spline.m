% By C. Mollet
%
% Function to plot a tensor product Spline on a uniform
% grid (but possibly different refinement levels)
%
% Input:
% ------
% 'c'     : n-D cektorized expansion coefficients
% 'T'     : Set of knots of B-Splines given as a cell element {T1,T2,...,Tn}
% 'k'     : Order of B-Splines in each direction
% 'lev'   : level of (uniform) refinement, e.g., lev = [4,5] results in a grid
%           with equidistant spacing 2^(-4) in x and 2^(-5) in y
%
% Output:
% -------
% 'Val'   : Vector/Matrix/Movieframes of Data points to plot

function Val = plot_Spline(c,T,k,lev)

% 1D case
if length(k)==1
    % Convert cell to vector
    T1d = cell2mat(T(1));
    % spacing
    h = 2^(-lev);
    j = 1;
    f = zeros(length(T1d(1):h:T1d(length(T1d))-h),1);
    for x = T1d(1):h:T1d(length(T1d))-h;
        f(j) = Nev(c,T,k,x);
        j = j+1;
    end
    plot(T1d(1):h:T1d(length(T1d))-h , f);
    Val = f;
end

% 2D case
if length(k)==2
    % Convert cell to vector
    T1 = cell2mat(T(1));
    T2 = cell2mat(T(2));
    % spacing
    h1 = 2^(-lev(1));
    h2 = 2^(-lev(2));
    f = zeros(length(T2(1):h2:T2(length(T2))-h2),length(T1(1):h1:T1(length(T1))-h1));
    for i = 0:2^(lev(1))-1
        t = i/(2^(lev(1)));
        for j = 0:2^(lev(2))-1
            x = j/(2^(lev(2)));
            f(j+1,i+1) = Nev(c,T,k,[t,x]);
        end
    end
    [X,Y] =  meshgrid(0:2^(-lev(1)):1-2^(-lev(1)),0:2^(-lev(2)):1-2^(-lev(2)));
    surf(X,Y,f);
    Val = f;
end

% Not optimized yet
% 3D case
if length(k)==3
    
    % Convert cell to vector
    T2 = cell2mat(T(2));
    T3 = cell2mat(T(3));
    % spacing
    h2 = 2^(-lev(2));
    h3 = 2^(-lev(3));
    % Special movie settings
    figure('Renderer','zbuffer');
    axis tight;
    set(gca,'NextPlot','replaceChildren');
    % Preallocate the struct array for the struct returned by getframe
    F(2^(lev(1))) = struct('cdata',[],'colormap',[]);
    
    for s=0:2^(lev(1))-1
        t = s/(2^(lev(1)));
        f = zeros(length(T3(1):h3:T3(length(T3))-h3),length(T2(1):h2:T2(length(T2))-h2));
        for i = 0:2^(lev(2))-1
            x = i/(2^(lev(2)));
            for j = 0:2^(lev(3))-1
                y = j/(2^(lev(3)));
                f(j+1,i+1) = Nev(c,T,k,[t,x,y]);
            end
        end
        [X,Y] =  meshgrid(0:2^(-lev(2)):1-2^(-lev(2)),0:2^(-lev(3)):1-2^(-lev(3)));
        surf(X,Y,f);
        xlabel('x');
        ylabel('y');
        %title(['Solution at time t=',num2str(t)]);
        F(s+1) = getframe;
    end
    movie(F,5,1);
    Val = F;
end

end