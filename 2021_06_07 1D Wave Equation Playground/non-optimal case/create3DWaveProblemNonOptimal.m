% Gets a problem configuration and creates a problem with all matrices, the
% right hand side and various variables for retrieving the solution 
function [p] = create3DWaveProblemNonOptimal(pC)


p = pC; % Copy all values from the problem configuration into the problem
p.T_time_ansatz = InitUniNodes(p.refinementLevel_time, p.bSplineOrder_time_ansatz);
p.T_time_test = InitUniNodes(p.refinementLevel_time, p.bSplineOrder_time_test);


p.M_time = StiffMat(p.T_time_ansatz, ... % Tsol
    p.T_time_test, ... % Ttest
    p.bSplineOrder_time_ansatz, ... %ksol
    p.bSplineOrder_time_test, ... %ktest
    p.offset_time_ansatz,... %offsetsol
    p.offset_time_test, ... %offsettest
    [0, 0] ... %diff
 );

% D_time also known as N_time
p.D_time = StiffMat(p.T_time_ansatz, ... % Tsol
    p.T_time_test, ... % Ttest
    p.bSplineOrder_time_ansatz, ... %ksol
    p.bSplineOrder_time_test, ... %ktest
    p.offset_time_ansatz,... %offsetsol
    p.offset_time_test, ... %offsettest
    [0, 2] ... %diff
    );

%% Space matrices
% Use the same refinement level in every space dimension. This could be
% altered but does not bring any advantage for the current use case

p.T_space_ansatz = InitUniNodes(p.refinementLevel_space, p.bSplineOrder_space_ansatz);
p.T_space_test = InitUniNodes(p.refinementLevel_space, p.bSplineOrder_space_test);
% offset_sol_space = [1, 1];
% offset_test_space = [1, 1];





p.M_space_local = StiffMat(p.T_space_ansatz, ... % Tsol
    p.T_space_test, ... % Ttest
    p.bSplineOrder_space_ansatz, ... %ksol
    p.bSplineOrder_space_test, ... %ktest
    p.offset_space_ansatz,... %offsetsol
    p.offset_space_test, ... %offsettest
    [0, 0] ... %diff
    );


% A_space also known as N_space
p.A_space_local = -StiffMat(p.T_space_ansatz, ... % Tsol
    p.T_space_test, ... % Ttest
    p.bSplineOrder_space_ansatz, ... %ksol
    p.bSplineOrder_space_test, ... %ktest
    p.offset_space_ansatz,... %offsetsol
    p.offset_space_test, ... %offsettest
    [0, 2] ... %diff
    );

p.M_space = kron3(p.M_space_local, p.M_space_local, p.M_space_local);
p.A_space = p.c^2 * kron3(p.A_space_local, p.M_space_local, p.M_space_local) ...
    + kron3(p.M_space_local, p.A_space_local, p.M_space_local) ...
    + kron3(p.M_space_local, p.M_space_local, p.A_space_local);







%% rhs

if ~isempty(p.f_time)
    
    ntest_time = length(p.T_time_test) - p.bSplineOrder_time_test;
    
    
    for i=1:length(p.f_time)
        rhs_time = zeros(ntest_time-(p.offset_time_test(1)+p.offset_time_test(2)),1);
        % time component
        
        
        for j=1+p.offset_time_test(1):ntest_time-p.offset_time_test(2)
            g = @(t) p.f_time{i}(t) * Ndiff(p.T_time_test, p.bSplineOrder_time_test,0,j,t);
            for step = 0:p.bSplineOrder_time_test-1
                rhs_time(j-p.offset_time_test(1)) = rhs_time(j-p.offset_time_test(1)) ...
                    + Gaussq(p.T_time_test(j+step),p.T_time_test(j+step+1),g,5);
            end
        end
        
        % first space component
        ntest_space = length(p.T_space_test) - p.bSplineOrder_space_test;
        rhs_space_x = zeros(ntest_space-(p.offset_space_test(1)+p.offset_space_test(2)),1);
        
        for j=1+p.offset_space_test(1):ntest_space-p.offset_space_test(2)
            g = @(x) p.f_space_x{i}(x) * Ndiff(p.T_space_test, p.bSplineOrder_space_test,0,j,x);
            for step = 0:p.bSplineOrder_space_test-1
                rhs_space_x(j-p.offset_space_test(1)) = rhs_space_x(j-p.offset_space_test(1)) ...
                    + Gaussq(p.T_space_test(j+step),p.T_space_test(j+step+1),g,5);
            end
        end
        
        % second space component
        rhs_space_y = zeros(ntest_space-(p.offset_space_test(1)+p.offset_space_test(2)),1);
        
        for j=1+p.offset_space_test(1):ntest_space-p.offset_space_test(2)
            g = @(y) p.f_space_y{i}(y) * Ndiff(p.T_space_test, p.bSplineOrder_space_test,0,j,y);
            for step = 0:p.bSplineOrder_space_test-1
                rhs_space_y(j-p.offset_space_test(1)) = rhs_space_y(j-p.offset_space_test(1)) ...
                    + Gaussq(p.T_space_test(j+step),p.T_space_test(j+step+1),g,5);
            end
        end
        
        % third space component
        rhs_space_z = zeros(ntest_space-(p.offset_space_test(1)+p.offset_space_test(2)),1);
        
        for j=1+p.offset_space_test(1):ntest_space-p.offset_space_test(2)
            g = @(z) p.f_space_z{i}(z) * Ndiff(p.T_space_test, p.bSplineOrder_space_test,0,j,z);
            for step = 0:p.bSplineOrder_space_test-1
                rhs_space_z(j-p.offset_space_test(1)) = rhs_space_z(j-p.offset_space_test(1)) ...
                    + Gaussq(p.T_space_test(j+step),p.T_space_test(j+step+1),g,5);
            end
        end
        
        if i == 1
            p.rhs = kron(rhs_time', kron(rhs_space_z, kron(rhs_space_y, rhs_space_x)));
        else
            p.rhs = p.rhs + kron(rhs_time', kron(rhs_space_z, kron(rhs_space_y, rhs_space_x)));
        end
    end
    
    
    %% Add the inital condtions
    for i=1+p.offset_time_test(1):length(p.T_time_test) - p.bSplineOrder_time_test - p.offset_time_test(2)
        
        if Ndiff(p.T_time_test, p.bSplineOrder_time_test,0,i,0) == 0 && ...
                Ndiff(p.T_time_test, p.bSplineOrder_time_test,1,i,0) == 0
            break
        end
        
        % first space component
        rhs_x_1 = sparse(ntest_space-(p.offset_space_test(1)+p.offset_space_test(2)),1);
        rhs_x_2 = sparse(ntest_space-(p.offset_space_test(1)+p.offset_space_test(2)),1);
        
        for j=1+p.offset_space_test(1):ntest_space-p.offset_space_test(2)
            g1 = @(x) p.u_1_x(x) * Ndiff(p.T_time_test, p.bSplineOrder_time_test,0,i,0) ...
                .* Ndiff(p.T_space_test, p.bSplineOrder_space_test,0,j,x);
            g2 = @(x) -p.u_0_x(x) * Ndiff(p.T_time_test, p.bSplineOrder_time_test,1,i,0) ...
                .* Ndiff(p.T_space_test, p.bSplineOrder_space_test,0,j,x);
            for step = 0:p.bSplineOrder_space_test-1
                rhs_x_1(j-p.offset_space_test(1)) = rhs_x_1(j-p.offset_space_test(1)) ...
                    + Gaussq(p.T_space_test(j+step),p.T_space_test(j+step+1),g1,5);
                rhs_x_2(j-p.offset_space_test(1)) = rhs_x_2(j-p.offset_space_test(1)) ...
                    + Gaussq(p.T_space_test(j+step),p.T_space_test(j+step+1),g2,5);
            end
        end
        
        % second space component
        rhs_y_1 = sparse(ntest_space-(p.offset_space_test(1)+p.offset_space_test(2)),1);
        rhs_y_2 = sparse(ntest_space-(p.offset_space_test(1)+p.offset_space_test(2)),1);
        
        for j=1+p.offset_space_test(1):ntest_space-p.offset_space_test(2)
            g1 = @(y) p.u_1_y(y) .* Ndiff(p.T_space_test, p.bSplineOrder_space_test,0,j,y);
            g2 = @(y) p.u_0_y(y) .* Ndiff(p.T_space_test, p.bSplineOrder_space_test,0,j,y);
            for step = 0:p.bSplineOrder_space_test-1
                rhs_y_1(j-p.offset_space_test(1)) = rhs_y_1(j-p.offset_space_test(1)) ...
                    + Gaussq(p.T_space_test(j+step),p.T_space_test(j+step+1),g1,5);
                rhs_y_2(j-p.offset_space_test(1)) = rhs_y_2(j-p.offset_space_test(1)) ...
                    + Gaussq(p.T_space_test(j+step),p.T_space_test(j+step+1),g2,5);
            end
        end
        
        % third space component
        rhs_z_1 = sparse(ntest_space-(p.offset_space_test(1)+p.offset_space_test(2)),1);
        rhs_z_2 = sparse(ntest_space-(p.offset_space_test(1)+p.offset_space_test(2)),1);
        
        for j=1+p.offset_space_test(1):ntest_space-p.offset_space_test(2)
            g1 = @(z) p.u_1_z(z) .* Ndiff(p.T_space_test, p.bSplineOrder_space_test,0,j,z);
            g2 = @(z) p.u_0_z(z) .* Ndiff(p.T_space_test, p.bSplineOrder_space_test,0,j,z);
            for step = 0:p.bSplineOrder_space_test-1
                rhs_z_1(j-p.offset_space_test(1)) = rhs_z_1(j-p.offset_space_test(1)) ...
                    + Gaussq(p.T_space_test(j+step),p.T_space_test(j+step+1),g1,5);
                rhs_z_2(j-p.offset_space_test(1)) = rhs_z_2(j-p.offset_space_test(1)) ...
                    + Gaussq(p.T_space_test(j+step),p.T_space_test(j+step+1),g2,5);
            end
        end
        
        p.rhs(:,i) = p.rhs(:,i) + kron(rhs_x_1, kron(rhs_y_1, rhs_z_1)) ...
            + kron(rhs_x_2, kron(rhs_y_2, rhs_z_2));
        
    end
end

% Calculate some additional values
     p.nsol_time = length(p.T_time_ansatz) - p.bSplineOrder_time_ansatz;
     p.nsol_space = length(p.T_space_ansatz) - p.bSplineOrder_space_ansatz;
     p.dim_time =  p.nsol_time-(p.offset_time_ansatz(1)+p.offset_time_ansatz(2));
end

function D = kron3(A,B,C)
D = kron(A,kron(B,C));
end