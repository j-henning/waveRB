function [Q_time, M_time, N_time, Q_space, M_space, N_space, rhs] ...
    = computeMatrices1D(refinementLevel_time, bSplineOrder_time, ...
    refinementLevel_space, bSplineOrder_space, f_time, f_space, ...
    u_0_x, u_1_x)

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

Q_space = StiffMat(T_space, ... % Tsol
    T_space, ... % Ttest
    bSplineOrder_space, ... %ksol
    bSplineOrder_space, ... %ktest
    offset_space_ansatz,... %offsetsol
    offset_space_test, ... %offsettest
    [2, 2] ... %diff
    );


M_space = StiffMat(T_space, ... % Tsol
    T_space, ... % Ttest
    bSplineOrder_space, ... %ksol
    bSplineOrder_space, ... %ktest
    offset_space_ansatz,... %offsetsol
    offset_space_test, ... %offsettest
    [0, 0] ... %diff
    );



N_space = -StiffMat(T_space, ... % Tsol
    T_space, ... % Ttest
    bSplineOrder_space, ... %ksol
    bSplineOrder_space, ... %ktest
    offset_space_ansatz,... %offsetsol
    offset_space_test, ... %offsettest
    [2, 0] ... %diff
    );




%% rhs
if length(f_time) ~= length(f_space)
    error('f_time and f_space should have the same dimensions')
end

rhs = [];

if ~isempty(f_time)
    
    ntest_time = length(T_time) - bSplineOrder_time;
    ntest_space = length(T_space) - bSplineOrder_space;
    
    for i=1:length(f_time)
        
        % time component
        
        rhs_time = sparse(ntest_time-(offset_time_test(1)+offset_time_test(2)),1);
        for j=1+offset_time_test(1):ntest_time-offset_time_test(2)
            g = @(t) f_time{i}(t) * Ndiff(T_time, bSplineOrder_time,0,j,t);
            for step = 0:bSplineOrder_time-1
                rhs_time(j-offset_time_test(1)) = rhs_time(j-offset_time_test(1)) ...
                    + Gaussq(T_time(j+step),T_time(j+step+1),g,5);
            end
        end
        
        % space component
        rhs_space_x = sparse(ntest_space-(offset_space_test(1)+offset_space_test(2)),1);
        
        for j=1+offset_space_test(1):ntest_space-offset_space_test(2)
            g = @(x) f_space{i}(x) * Ndiff(T_space, bSplineOrder_space,0,j,x);
            for step = 0:bSplineOrder_space-1
                rhs_space_x(j-offset_space_test(1)) = rhs_space_x(j-offset_space_test(1)) ...
                    + Gaussq(T_space(j+step),T_space(j+step+1),g,5);
            end
        end
        
                
        if i == 1
            rhs = kron(rhs_time', rhs_space_x);
        else
            rhs = rhs + kron(rhs_time', rhs_space_x);
        end
    end
    
    %% Add the inital condtions
    for i=1+offset_time_test(1):length(T_time) - bSplineOrder_time - offset_time_test(2)
        
        if Ndiff(T_time, bSplineOrder_time,0,i,0) == 0 && ...
                 Ndiff(T_time, bSplineOrder_time,1,i,0) == 0
             break
        end
        
        % space component
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
        
       
 
        rhs(:,i) = rhs(:,i) + rhs_x_1 + rhs_x_2;
        
    end
    
    
end
end