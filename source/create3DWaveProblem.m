% Gets a problem configuration and creates a problem with all matrices, the
% right hand side and various variables for retrieving the solution
function [p] = create3DWaveProblem(pC)


p = pC; % Copy all values from the problem configuration into the problem
p.T_time = InitUniNodes(p.refinementLevel_time, p.bSplineOrder_time);

p.Q_time = StiffMat(p.T_time, ... % Tsol
    p.T_time, ... % Ttest
    p.bSplineOrder_time, ... %ksol
    p.bSplineOrder_time, ... %ktest
    p.offset_time_ansatz,... %offsetsol
    p.offset_time_test, ... %offsettest
    [2, 2] ... %diff
    );

p.M_time = StiffMat(p.T_time, ... % Tsol
    p.T_time, ... % Ttest
    p.bSplineOrder_time, ... %ksol
    p.bSplineOrder_time, ... %ktest
    p.offset_time_ansatz,... %offsetsol
    p.offset_time_test, ... %offsettest
    [0, 0] ... %diff
    );

% D_time also known as N_time
p.D_time = StiffMat(p.T_time, ... % Tsol
    p.T_time, ... % Ttest
    p.bSplineOrder_time, ... %ksol
    p.bSplineOrder_time, ... %ktest
    p.offset_time_ansatz,... %offsetsol
    p.offset_time_test, ... %offsettest
    [2, 0] ... %diff
    );

%% Space matrices
% Use the same refinement level in every space dimension. This could be
% altered but does not bring any advantage for the current use case

p.T_space = InitUniNodes(p.refinementLevel_space, p.bSplineOrder_space);
% offset_sol_space = [1, 1];
% offset_test_space = [1, 1];

p.Q_space_local = StiffMat(p.T_space, ... % Tsol
    p.T_space, ... % Ttest
    p.bSplineOrder_space, ... %ksol
    p.bSplineOrder_space, ... %ktest
    p.offset_space_ansatz,... %offsetsol
    p.offset_space_test, ... %offsettest
    [2, 2] ... %diff
    );


p.M_space_local = StiffMat(p.T_space, ... % Tsol
    p.T_space, ... % Ttest
    p.bSplineOrder_space, ... %ksol
    p.bSplineOrder_space, ... %ktest
    p.offset_space_ansatz,... %offsetsol
    p.offset_space_test, ... %offsettest
    [0, 0] ... %diff
    );


% A_space_local also known as N_space_local
% IMPORTANT: FLIP THE SIGN
% TODO: Check sign
p.A_space_local = StiffMat(p.T_space, ... % Tsol
    p.T_space, ... % Ttest
    p.bSplineOrder_space, ... %ksol
    p.bSplineOrder_space, ... %ktest
    p.offset_space_ansatz,... %offsetsol
    p.offset_space_test, ... %offsettest
    [2, 0] ... %diff
    );

p.Q_space = p.mu^2 * (kron3(p.Q_space_local, p.M_space_local, p.M_space_local) ...
    + kron3(p.A_space_local, p.M_space_local, p.A_space_local') ...
    + kron3(p.A_space_local, p.A_space_local', p.M_space_local) ...
    + kron3(p.A_space_local', p.M_space_local, p.A_space_local) ...
    + kron3(p.M_space_local, p.M_space_local, p.Q_space_local) ...
    + kron3(p.M_space_local, p.A_space_local', p.A_space_local) ...
    + kron3(p.A_space_local', p.A_space_local, p.M_space_local) ...
    + kron3(p.M_space_local, p.A_space_local, p.A_space_local') ...
    + kron3(p.M_space_local, p.Q_space_local, p.M_space_local));

p.M_space = kron3(p.M_space_local, p.M_space_local, p.M_space_local);

p.A_space = -p.mu *  ( kron3(p.M_space_local, p.M_space_local, p.A_space_local) ...
    + kron3(p.M_space_local, p.A_space_local, p.M_space_local) ...
    + kron3(p.A_space_local, p.M_space_local, p.M_space_local));

%% rhs

p.ntest_time = length(p.T_time) - p.bSplineOrder_time;
p.ntest_space = length(p.T_space) - p.bSplineOrder_space;
p.rhs = zeros((p.ntest_space-sum(p.offset_space_test))^3,p.ntest_time - sum(p.offset_time_test));

% if length(p.f_time) ~= length(p.f_space_x) || length(p.f_time) ~= length(p.f_space_y)
%     error('f_time and f_space should have the same dimensions')
% end

p.ntest_time = length(p.T_time) - p.bSplineOrder_time;
p.ntest_space = length(p.T_space) - p.bSplineOrder_space;


if ~isempty(p.f_time)
    
    if size(p.f_space,2) == 3
        % The right hand side is provided in a nested form, so we can use
        % the much faster integration
        
        for i=1:length(p.f_time)
            
            % time component
            
            rhs_time = zeros(p.ntest_time-(p.offset_time_test(1)+p.offset_time_test(2)),1);
            for j=1+p.offset_time_test(1):p.ntest_time-p.offset_time_test(2)
                g = @(t) p.f_time{i}(t) * Ndiff(p.T_time, p.bSplineOrder_time,0,j,(t==0)*(t+eps) + (t==1)*(t-eps) + (t > 0 && t < 1)*t);
                for step = 0:p.bSplineOrder_time-1
                    rhs_time(j-p.offset_time_test(1)) = rhs_time(j-p.offset_time_test(1)) ...
                        + Gaussq(p.T_time(j+step),p.T_time(j+step+1),g,5);
                end
            end
            
            % first space component
            rhs_space_x = zeros(p.ntest_space-(p.offset_space_test(1)+p.offset_space_test(2)),1);
            
            for j=1+p.offset_space_test(1):p.ntest_space-p.offset_space_test(2)
                g = @(x) p.f_space{i,1}(x) * Ndiff(p.T_space, p.bSplineOrder_space,0,j,(x==0)*(x+eps) + (x==1)*(x-eps) + (x > 0 && x < 1)*x);
                for step = 0:p.bSplineOrder_space-1
                    rhs_space_x(j-p.offset_space_test(1)) = rhs_space_x(j-p.offset_space_test(1)) ...
                        + Gaussq(p.T_space(j+step),p.T_space(j+step+1),g,5);
                end
            end
            
            % second space component
            p.ntest_space = length(p.T_space) - p.bSplineOrder_space;
            rhs_space_y = zeros(p.ntest_space-(p.offset_space_test(1)+p.offset_space_test(2)),1);
            
            for j=1+p.offset_space_test(1):p.ntest_space-p.offset_space_test(2)
                g = @(y) p.f_space{i,2}(y) * Ndiff(p.T_space, p.bSplineOrder_space,0,j,(y==0)*(y+eps) + (y==1)*(y-eps) + (y > 0 && y < 1)*y);
                for step = 0:p.bSplineOrder_space-1
                    rhs_space_y(j-p.offset_space_test(1)) = rhs_space_y(j-p.offset_space_test(1)) ...
                        + Gaussq(p.T_space(j+step),p.T_space(j+step+1),g,5);
                end
            end
            
            % third space component
            p.ntest_space = length(p.T_space) - p.bSplineOrder_space;
            rhs_space_z = zeros(p.ntest_space-(p.offset_space_test(1)+p.offset_space_test(2)),1);
            
            for j=1+p.offset_space_test(1):p.ntest_space-p.offset_space_test(2)
                g = @(y) p.f_space{i,3}(y) * Ndiff(p.T_space, p.bSplineOrder_space,0,j,(y==0)*(y+eps) + (y==1)*(y-eps) + (y > 0 && y < 1)*y);
                for step = 0:p.bSplineOrder_space-1
                    rhs_space_z(j-p.offset_space_test(1)) = rhs_space_z(j-p.offset_space_test(1)) ...
                        + Gaussq(p.T_space(j+step),p.T_space(j+step+1),g,5);
                end
            end
            
            if i == 1
                p.rhs = kron(rhs_time', kron(rhs_space_z, kron(rhs_space_y, rhs_space_x)));
            else
                p.rhs = p.rhs + kron(rhs_time', kron(rhs_space_z, kron(rhs_space_y, rhs_space_x)));
            end
        end
    
        
    
    else % The right hand side is in nested form, so we have to do more work
    for i=1:size(p.f_time,1)
        
        % time component
        
        rhs_time = zeros(p.ntest_time-(p.offset_time_test(1)+p.offset_time_test(2)),1);
        for j=1+p.offset_time_test(1):p.ntest_time-p.offset_time_test(2)
            g = @(t) p.f_time{i}(t) * Ndiff(p.T_time, p.bSplineOrder_time,0,j,(t==0)*(t+eps) + (t==1)*(t-eps) + (t > 0 && t < 1)*t);
            for step = 0:p.bSplineOrder_time-1
                rhs_time(j-p.offset_time_test(1)) = rhs_time(j-p.offset_time_test(1)) ...
                    + Gaussq(p.T_time(j+step),p.T_time(j+step+1),g,5);
            end
        end
        
        % Space components
        F_space = zeros(p.ntest_space-sum(p.offset_space_test), ...
            p.ntest_space-sum(p.offset_space_test), ...
            p.ntest_space-sum(p.offset_space_test));
        
        
        for j=1+p.offset_space_test(1):p.ntest_space-p.offset_space_test(2)
            for k=1+p.offset_space_test(1):p.ntest_space-p.offset_space_test(2)
                for l=1+p.offset_space_test(1):p.ntest_space-p.offset_space_test(2)
                    
                    
                    g = @(z)@(y)@(x) p.f_space{i}(x,y,z) .* ...
                        Ndiff(p.T_space, p.bSplineOrder_space,0,j,x) ...
                        .* Ndiff(p.T_space, p.bSplineOrder_space,0,k,y) ...
                        .* Ndiff(p.T_space, p.bSplineOrder_space,0,l,z);
                    
                    % Integrate
                    for stepX = 0:p.bSplineOrder_space-1
                        for stepY = 0:p.bSplineOrder_space-1
                            for stepZ = 0:p.bSplineOrder_space-1
                                
                                a = [p.T_space(j+stepX), p.T_space(k+stepY), p.T_space(l+stepZ)];
                                b = [p.T_space(j+stepX+1), p.T_space(k+stepY+1), p.T_space(l+stepZ+1)];
                                
                                
                                % Check if u0 is zero in the center of the
                                % inteval, as well as on the boundary
                                % if yes assume that the integral is zero.
                                % This is done for performance reasons
                                center = 0.5 * (a+b);
                                if abs(p.f_space{i}(center(1),center(2),center(3))) < eps && ...
                                        abs(p.f_space{i}(a(1),a(2),a(3))) < eps && ...
                                        abs(p.f_space{i}(a(1),a(2),b(3))) < eps && ...
                                        abs(p.f_space{i}(a(1),b(2),a(3))) < eps && ...
                                        abs(p.f_space{i}(a(1),b(2),b(3))) < eps && ...
                                        abs(p.f_space{i}(b(1),b(2),b(3))) < eps && ...
                                        abs(p.f_space{i}(b(1),b(2),a(3))) < eps && ...
                                        abs(p.f_space{i}(b(1),a(2),b(3))) < eps && ...
                                        abs(p.f_space{i}(b(1),a(2),a(3))) < eps
                                    continue;
                                end
                                
                                % Integrate from a to b
                                F_space(j-p.offset_space_test(1), ...
                                    k-p.offset_space_test(1), ...
                                    l-p.offset_space_test(1)) = ...
                                    F_space(j-p.offset_space_test(1), ...
                                    k-p.offset_space_test(1), ...
                                    l-p.offset_space_test(1))  + ...
                                    Gaussq(a,b,g,3);
                            end
                        end
                    end
                    
                end
            end
        end

        p.rhs = p.rhs + kron(rhs_time', F_space(:));
        
    end
    end
end











if ~isempty(p.u_0)
    
    
    
    U0 = zeros(p.ntest_space-sum(p.offset_space_test), ...
        p.ntest_space-sum(p.offset_space_test), ...
        p.ntest_space-sum(p.offset_space_test));
    
    
    parfor j=1+p.offset_space_test(1):p.ntest_space-p.offset_space_test(2)
        U0j = zeros(p.ntest_space-sum(p.offset_space_test), ...
        p.ntest_space-sum(p.offset_space_test));
    
        for k=1+p.offset_space_test(1):p.ntest_space-p.offset_space_test(2)
            for l=1+p.offset_space_test(1):p.ntest_space-p.offset_space_test(2)
                
                
                g = @(z)@(y)@(x) p.u_0(x,y,z) .* ...
                    Ndiff(p.T_space, p.bSplineOrder_space,0,j,x) ...
                    .* Ndiff(p.T_space, p.bSplineOrder_space,0,k,y) ...
                    .* Ndiff(p.T_space, p.bSplineOrder_space,0,l,z);
                
                % Integrate
                for stepX = 0:p.bSplineOrder_space-1
                    for stepY = 0:p.bSplineOrder_space-1
                        for stepZ = 0:p.bSplineOrder_space-1
                            
                            a = [p.T_space(j+stepX), p.T_space(k+stepY), p.T_space(l+stepZ)];
                            b = [p.T_space(j+stepX+1), p.T_space(k+stepY+1), p.T_space(l+stepZ+1)];
                            
                            
                            % Check if u0 is zero in the center of the
                            % inteval, as well as on the boundary
                            % if yes assume that the integral is zero.
                            % This is done for performance reasons
                            center = 0.5 * (a+b);
                            if abs(p.u_0(center(1),center(2),center(3))) < eps && ...
                                    abs(p.u_0(a(1),a(2),a(3))) < eps && ...
                                    abs(p.u_0(a(1),a(2),b(3))) < eps && ...
                                    abs(p.u_0(a(1),b(2),a(3))) < eps && ...
                                    abs(p.u_0(a(1),b(2),b(3))) < eps && ...
                                    abs(p.u_0(b(1),b(2),b(3))) < eps && ...
                                    abs(p.u_0(b(1),b(2),a(3))) < eps && ...
                                    abs(p.u_0(b(1),a(2),b(3))) < eps && ...
                                    abs(p.u_0(b(1),a(2),a(3))) < eps
                                continue;
                            end
                            
                            % Integrate from a to b
%                             U0(j-p.offset_space_test(1), ...
%                                 k-p.offset_space_test(1), ...
%                                 l-p.offset_space_test(1)) = ...
%                                 U0(j-p.offset_space_test(1), ...
%                                 k-p.offset_space_test(1), ...
%                                 l-p.offset_space_test(1))  + ...
%                                 Gaussq(a,b,g,3);

U0j(k-p.offset_space_test(1), ...
    l-p.offset_space_test(1)) = ...
    U0j(k-p.offset_space_test(1), ...
    l-p.offset_space_test(1))  + ...
    Gaussq(a,b,g,3);
                        end
                    end
                end
                U0(j,:,:) = U0j;
            end
        end
    end
    
    
    for i=1+p.offset_time_test(1):length(p.T_time) - p.bSplineOrder_time - p.offset_time_test(2)
        v_dot0 = Ndiff(p.T_time, p.bSplineOrder_time,1,i,0);
        
        if v_dot0 == 0
            continue;
        end
        
        p.rhs(:,i) = p.rhs(:,i) + reshape(-v_dot0 * U0, [numel(U0) 1]);
    end
    
end

% And now the same for u_1 (for performance reason it is better to use
% different loops, with the drawback that we have to copy the whole
% code

if ~isempty(p.u_1)
    
    
    
    U1 = zeros(p.ntest_space-sum(p.offset_space_test), ...
        p.ntest_space-sum(p.offset_space_test), ...
        p.ntest_space-sum(p.offset_space_test));
    
    for j=1+p.offset_space_test(1):p.ntest_space-p.offset_space_test(2)
        for k=1+p.offset_space_test(1):p.ntest_space-p.offset_space_test(2)
            for l=1+p.offset_space_test(1):p.ntest_space-p.offset_space_test(2)
                
                
                g = @(z)@(y)@(x) p.u_1(x,y,z) .* ...
                    Ndiff(p.T_space, p.bSplineOrder_space,0,j,x) ...
                    .* Ndiff(p.T_space, p.bSplineOrder_space,0,k,y) ...
                    .* Ndiff(p.T_space, p.bSplineOrder_space,0,l,z);
                
                % Integrate
                for stepX = 0:p.bSplineOrder_space-1
                    for stepY = 0:p.bSplineOrder_space-1
                        for stepZ = 0:p.bSplineOrder_space-1
                            
                            a = [p.T_space(j+stepX), p.T_space(k+stepY), p.T_space(l+stepZ)];
                            b = [p.T_space(j+stepX+1), p.T_space(k+stepY+1), p.T_space(l+stepZ+1)];
                            
                            
                            % Check if u0 is zero in the center of the
                            % inteval, as well as on the boundary
                            % if yes assume that the integral is zero.
                            % This is done for performance reasons
                            center = 0.5 * (a+b);
                            if abs(p.u_1(center(1),center(2),center(3))) < eps && ...
                                    abs(p.u_1(a(1),a(2),a(3))) < eps && ...
                                    abs(p.u_1(a(1),a(2),b(3))) < eps && ...
                                    abs(p.u_1(a(1),b(2),a(3))) < eps && ...
                                    abs(p.u_1(a(1),b(2),b(3))) < eps && ...
                                    abs(p.u_1(b(1),b(2),b(3))) < eps && ...
                                    abs(p.u_1(b(1),b(2),a(3))) < eps && ...
                                    abs(p.u_1(b(1),a(2),b(3))) < eps && ...
                                    abs(p.u_1(b(1),a(2),a(3))) < eps
                                continue;
                            end
                            
                            % Integrate from a to b
                            U1(j-p.offset_space_test(1), ...
                                k-p.offset_space_test(1), ...
                                l-p.offset_space_test(1)) = ...
                                U1(j-p.offset_space_test(1), ...
                                k-p.offset_space_test(1), ...
                                l-p.offset_space_test(1))  + ...
                                Gaussq(a,b,g,3);
                        end
                    end
                end
                
            end
        end
    end
    
    
    for i=1+p.offset_time_test(1):length(p.T_time) - p.bSplineOrder_time - p.offset_time_test(2)
        v_0 = Ndiff(p.T_time, p.bSplineOrder_time,0,i,0);
        
        if v_0 == 0
            continue;
        end
        
        p.rhs(:,i) = p.rhs(:,i) + reshape(v_0 * U0, [numel(U0) 1]);
    end
    
end

% Calculate some additional values
p.nsol_time = length(p.T_time) - p.bSplineOrder_time;
p.nsol_space = length(p.T_space) - p.bSplineOrder_space;
p.dim_time =  p.nsol_time-(p.offset_time_ansatz(1)+p.offset_time_ansatz(2));
end

function sol = kron3(A,B,C)
sol = kron(A, kron(B,C));
end