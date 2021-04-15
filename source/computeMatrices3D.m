function [Q_time, M_time, N_time, Q_space, M_space, N_space, rhs,M,Q,N] ...
    = computeMatrices3D(refinementLevel_time, bSplineOrder_time, ...
    refinementLevel_space, bSplineOrder_space, ...
    f_time, f_space_x, f_space_y, f_space_z, ...
    u_0_x, u_0_y, u_0_z, u_1_x, u_1_y, u_1_z)

% if refinementLevel_time <= 13 && bSplineOrder_time == 4
%     % Load the already assembled time matrices
%     load(['Matrices/time' num2str(refinementLevel_time)])
% else
    T_time = InitUniNodes(refinementLevel_time, bSplineOrder_time);
%     offset_time_ansatz = [0, bSplineOrder_time-1];
%     offset_time_test = [0, bSplineOrder_time-1];
    offset_time_ansatz = [0, 2];
    offset_time_test = [0, 2];
    
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
% end


%% Space matrices
% if refinementLevel_space <= 5 && bSplineOrder_space == 4
%     % Load the already assembled time matrices
%     load(['Matrices/space' num2str(refinementLevel_space)])
% else
    T_space = InitUniNodes(refinementLevel_space, bSplineOrder_space);
    offset_space_ansatz = [1, 1];
    offset_space_test = [1, 1];
    
    Q = StiffMat(T_space, ... % Tsol
        T_space, ... % Ttest
        bSplineOrder_space, ... %ksol
        bSplineOrder_space, ... %ktest
        offset_space_ansatz,... %offsetsol
        offset_space_test, ... %offsettest
        [2, 2] ... %diff
        );
    
    M = StiffMat(T_space, ... % Tsol
        T_space, ... % Ttest
        bSplineOrder_space, ... %ksol
        bSplineOrder_space, ... %ktest
        offset_space_ansatz,... %offsetsol
        offset_space_test, ... %offsettest
        [0, 0] ... %diff
        );
    
    N = -StiffMat(T_space, ... % Tsol
        T_space, ... % Ttest
        bSplineOrder_space, ... %ksol
        bSplineOrder_space, ... %ktest
        offset_space_ansatz,... %offsetsol
        offset_space_test, ... %offsettest
        [2, 0] ... %diff
        );
    
    
    Q_space = kron3(Q, M, M) + kron3(N, M, N') + kron3(N, N', M) ...
        + kron3(N', M, N) + kron3(M, M, Q) + kron3(M, N', N) ...
        + kron3(N', N, M) + kron3(M, N, N') + kron3(M, Q, M);
    
    M_space = kron3(M, M, M);
    N_space = kron3(M, M, N) ...
        + kron3(M, N, M) ...
        + kron3(N, M, M);
% end
%% rhs
if length(f_time) ~= length(f_space_x) ...
        || length(f_time) ~= length(f_space_y) ...
        || length(f_time) ~= length(f_space_z)
    error('f_time and f_space should have the same dimensions')
end

rhs = [];

if ~isempty(f_time)
    
    ntest_time = length(T_time) - bSplineOrder_time;
    
    
    for i=1:length(f_time)
        rhs_time = zeros(ntest_time-(offset_time_test(1)+offset_time_test(2)),1);
        % time component
        
        
        for j=1+offset_time_test(1):ntest_time-offset_time_test(2)
            g = @(t) f_time{i}(t) * Ndiff(T_time, bSplineOrder_time,0,j,t);
            for step = 0:bSplineOrder_time-1
                rhs_time(j-offset_time_test(1)) = rhs_time(j-offset_time_test(1)) ...
                    + Gaussq(T_time(j+step),T_time(j+step+1),g,5);
            end
        end
        
        % first space component
        ntest_space = length(T_space) - bSplineOrder_space;
        rhs_space_x = zeros(ntest_space-(offset_space_test(1)+offset_space_test(2)),1);
        
        for j=1+offset_space_test(1):ntest_space-offset_space_test(2)
            g = @(x) f_space_x{i}(x) * Ndiff(T_space, bSplineOrder_space,0,j,x);
            for step = 0:bSplineOrder_space-1
                rhs_space_x(j-offset_space_test(1)) = rhs_space_x(j-offset_space_test(1)) ...
                    + Gaussq(T_space(j+step),T_space(j+step+1),g,5);
            end
        end
        
        % second space component
        rhs_space_y = zeros(ntest_space-(offset_space_test(1)+offset_space_test(2)),1);
        
        for j=1+offset_space_test(1):ntest_space-offset_space_test(2)
            g = @(y) f_space_y{i}(y) * Ndiff(T_space, bSplineOrder_space,0,j,y);
            for step = 0:bSplineOrder_space-1
                rhs_space_y(j-offset_space_test(1)) = rhs_space_y(j-offset_space_test(1)) ...
                    + Gaussq(T_space(j+step),T_space(j+step+1),g,5);
            end
        end
        
        % third space component
        rhs_space_z = zeros(ntest_space-(offset_space_test(1)+offset_space_test(2)),1);
        
        for j=1+offset_space_test(1):ntest_space-offset_space_test(2)
            g = @(z) f_space_z{i}(z) * Ndiff(T_space, bSplineOrder_space,0,j,z);
            for step = 0:bSplineOrder_space-1
                rhs_space_z(j-offset_space_test(1)) = rhs_space_z(j-offset_space_test(1)) ...
                    + Gaussq(T_space(j+step),T_space(j+step+1),g,5);
            end
        end
        
        if i == 1
             rhs = kron(rhs_time', kron(rhs_space_z, kron(rhs_space_y, rhs_space_x)));
        else
           rhs = rhs + kron(rhs_time', kron(rhs_space_z, kron(rhs_space_y, rhs_space_x)));
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
        
         % third space component
        rhs_z_1 = sparse(ntest_space-(offset_space_test(1)+offset_space_test(2)),1);
        rhs_z_2 = sparse(ntest_space-(offset_space_test(1)+offset_space_test(2)),1);
        
        for j=1+offset_space_test(1):ntest_space-offset_space_test(2)
            g1 = @(z) u_1_z(z) .* Ndiff(T_space, bSplineOrder_space,0,j,z);
            g2 = @(z) u_0_z(z) .* Ndiff(T_space, bSplineOrder_space,0,j,z);
            for step = 0:bSplineOrder_space-1
                rhs_z_1(j-offset_space_test(1)) = rhs_z_1(j-offset_space_test(1)) ...
                    + Gaussq(T_space(j+step),T_space(j+step+1),g1,5);
                rhs_z_2(j-offset_space_test(1)) = rhs_z_2(j-offset_space_test(1)) ...
                    + Gaussq(T_space(j+step),T_space(j+step+1),g2,5);
            end
        end
 
        rhs(:,i) = rhs(:,i) + kron(rhs_x_1, kron(rhs_y_1, rhs_z_1)) ...
            + kron(rhs_x_2, kron(rhs_y_2, rhs_z_2));
        
    end
end
end

function sol = kron3(A,B,C)
sol = kron(A, kron(B,C));
end