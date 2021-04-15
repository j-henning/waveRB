function [Q_time, M_time, N_time, Q_space, M_space, N_space, rhs, ...
    M_space_local,Q_space_local,N_space_local] ...
    = computeMatrices2D(refinementLevel_time, bSplineOrder_time, ...
    refinementLevel_space, bSplineOrder_space, f_time, f_space_x, ...
    f_space_y, u_0_x, u_0_y, u_1_x, u_1_y)

T_time = InitUniNodes(refinementLevel_time, bSplineOrder_time);
% offset_sol_time = [0, bSplineOrder_time-1];
% offset_time_test = [0, bSplineOrder_time-1];

offset_time_ansatz = [0, 2];
offset_time_test = [0, 2];
offset_space_ansatz = [1, 1];
offset_space_test = [1, 1];

Q_time = StiffMat(T_time, ... % Tsol
    T_time, ... % Ttest
    bSplineOrder_time, ... %ksol
    bSplineOrder_time, ... %ktest
    offset_time_ansatz,... %offsetsol
    offset_time_test, ... %offsettest
    [2, 2] ... %diff
    );

M_time = StiffMat(T_time, ... % Tsol
    T_time, ... % Ttest
    bSplineOrder_time, ... %ksol
    bSplineOrder_time, ... %ktest
    offset_time_ansatz,... %offsetsol
    offset_time_test, ... %offsettest
    [0, 0] ... %diff
    );
N_time = StiffMat(T_time, ... % Tsol
    T_time, ... % Ttest
    bSplineOrder_time, ... %ksol
    bSplineOrder_time, ... %ktest
    offset_time_ansatz,... %offsetsol
    offset_time_test, ... %offsettest
    [2, 0] ... %diff
    );



%% Space matrices
% Use the same refinement level in every space dimension. This could be
% altered but does not bring any advantage for the current use case

T_space = InitUniNodes(refinementLevel_space, bSplineOrder_space);
% offset_sol_space = [1, 1];
% offset_test_space = [1, 1];

Q_space_local = StiffMat(T_space, ... % Tsol
    T_space, ... % Ttest
    bSplineOrder_space, ... %ksol
    bSplineOrder_space, ... %ktest
    offset_space_ansatz,... %offsetsol
    offset_space_test, ... %offsettest
    [2, 2] ... %diff
    );


M_space_local = StiffMat(T_space, ... % Tsol
    T_space, ... % Ttest
    bSplineOrder_space, ... %ksol
    bSplineOrder_space, ... %ktest
    offset_space_ansatz,... %offsetsol
    offset_space_test, ... %offsettest
    [0, 0] ... %diff
    );



N_space_local = StiffMat(T_space, ... % Tsol
    T_space, ... % Ttest
    bSplineOrder_space, ... %ksol
    bSplineOrder_space, ... %ktest
    offset_space_ansatz,... %offsetsol
    offset_space_test, ... %offsettest
    [2, 0] ... %diff
    );



Q_space = kron(Q_space_local, M_space_local) + kron(N_space_local, N_space_local') ...
    + kron(N_space_local', N_space_local) ...
    + kron(M_space_local, Q_space_local);

M_space = kron(M_space_local, M_space_local);
N_space = -kron(N_space_local, M_space_local) - kron(M_space_local, N_space_local);

%% rhs
if length(f_time) ~= length(f_space_x) || length(f_time) ~= length(f_space_y)
    error('f_time and f_space should have the same dimensions')
end

rhs = [];

if ~isempty(f_time)
    
    ntest_time = length(T_time) - bSplineOrder_time;
    ntest_space = length(T_space) - bSplineOrder_space;
    
    for i=1:length(f_time)
        
        % time component
        
        rhs_time = zeros(ntest_time-(offset_time_test(1)+offset_time_test(2)),1);
        for j=1+offset_time_test(1):ntest_time-offset_time_test(2)
            g = @(t) f_time{i}(t) * Ndiff(T_time, bSplineOrder_time,0,j,(t==0)*(t+eps) + (t==1)*(t-eps) + (t > 0 && t < 1)*t);
            for step = 0:bSplineOrder_time-1
                rhs_time(j-offset_time_test(1)) = rhs_time(j-offset_time_test(1)) ...
                    + Gaussq(T_time(j+step),T_time(j+step+1),g,5);
            end
        end
        
        % first space component
        rhs_space_x = zeros(ntest_space-(offset_space_test(1)+offset_space_test(2)),1);
        
        for j=1+offset_space_test(1):ntest_space-offset_space_test(2)
            g = @(x) f_space_x{i}(x) * Ndiff(T_space, bSplineOrder_space,0,j,(x==0)*(x+eps) + (x==1)*(x-eps) + (x > 0 && x < 1)*x);
            for step = 0:bSplineOrder_space-1
                rhs_space_x(j-offset_space_test(1)) = rhs_space_x(j-offset_space_test(1)) ...
                    + Gaussq(T_space(j+step),T_space(j+step+1),g,5);
            end
        end
        
        % second space component
        ntest_space = length(T_space) - bSplineOrder_space;
        rhs_space_y = zeros(ntest_space-(offset_space_test(1)+offset_space_test(2)),1);
        
        for j=1+offset_space_test(1):ntest_space-offset_space_test(2)
            g = @(y) f_space_y{i}(y) * Ndiff(T_space, bSplineOrder_space,0,j,(y==0)*(y+eps) + (y==1)*(y-eps) + (y > 0 && y < 1)*y);
            for step = 0:bSplineOrder_space-1
                rhs_space_y(j-offset_space_test(1)) = rhs_space_y(j-offset_space_test(1)) ...
                    + Gaussq(T_space(j+step),T_space(j+step+1),g,5);
            end
        end
        
        if i == 1
            rhs = kron(rhs_time', kron(rhs_space_x, rhs_space_y));
        else
            rhs = rhs + kron(rhs_time', kron(rhs_space_x, rhs_space_y));
        end
    end
    
    %% Add the inital condtions
    for i=1+offset_time_test(1):length(T_time) - bSplineOrder_time - offset_time_test(2)
        
        if Ndiff(T_time, bSplineOrder_time,0,i,0) == 0 && ...
                 Ndiff(T_time, bSplineOrder_time,1,i,0) == 0
             break
        end
        
        % first space component
        rhs_x_1 = sparse(ntest_space-(offset_space_test(1)+offset_space_test(2)),1);
        rhs_x_2 = sparse(ntest_space-(offset_space_test(1)+offset_space_test(2)),1);
        
        for j=1+offset_space_test(1):ntest_space-offset_space_test(2)
            g1 = @(x) u_1_x(x) * Ndiff(T_time, bSplineOrder_time,0,i,0) ...
                .* Ndiff(T_space, bSplineOrder_space,0,j,x);
            g2 = @(x) -u_0_x(x) * Ndiff(T_time, bSplineOrder_time,1,i,0) ...
                .* Ndiff(T_space, bSplineOrder_space,0,j,x);
            for step = 0:bSplineOrder_space-1
                rhs_x_1(j-offset_space_test(1)) = rhs_x_1(j-offset_space_test(1)) ...
                    + Gaussq(T_space(j+step),T_space(j+step+1),g1,5);
                rhs_x_2(j-offset_space_test(1)) = rhs_x_2(j-offset_space_test(1)) ...
                    + Gaussq(T_space(j+step),T_space(j+step+1),g2,5);
            end
        end
        
        % second space component
        rhs_y_1 = sparse(ntest_space-(offset_space_test(1)+offset_space_test(2)),1);
        rhs_y_2 = sparse(ntest_space-(offset_space_test(1)+offset_space_test(2)),1);
        
        for j=1+offset_space_test(1):ntest_space-offset_space_test(2)
            g1 = @(y) u_1_y(y) .* Ndiff(T_space, bSplineOrder_space,0,j,y);
            g2 = @(y) u_0_y(y) .* Ndiff(T_space, bSplineOrder_space,0,j,y);
            for step = 0:bSplineOrder_space-1
                rhs_y_1(j-offset_space_test(1)) = rhs_y_1(j-offset_space_test(1)) ...
                    + Gaussq(T_space(j+step),T_space(j+step+1),g1,5);
                rhs_y_2(j-offset_space_test(1)) = rhs_y_2(j-offset_space_test(1)) ...
                    + Gaussq(T_space(j+step),T_space(j+step+1),g2,5);
            end
        end
 
        rhs(:,i) = rhs(:,i) + kron(rhs_x_1, rhs_y_1) + kron(rhs_x_2, rhs_y_2);
        
    end
    
    
end
end