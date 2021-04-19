% Gets a problem configuration and creates a problem with all matrices, the
% right hand side and various variables for retrieving the solution
function [p] = create3DWaveProblem(pC)

p = pC; % Copy all values from the problem configuration into the problem
p.T_time = InitUniNodes(p.refinementLevel_time, p.bSplineOrder_time);

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

p.T_space = InitUniNodes(p.refinementLevel_space, p.bSplineOrder_space);
p.offset_space_ansatz = [1, 1];
p.offset_space_test = [1, 1];

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
% TODO: Sign is flipped, which is no problem for B, but maybe for the
% preconditioner?
p.A_space_local = -StiffMat(p.T_space, ... % Tsol
    p.T_space, ... % Ttest
    p.bSplineOrder_space, ... %ksol
    p.bSplineOrder_space, ... %ktest
    p.offset_space_ansatz,... %offsetsol
    p.offset_space_test, ... %offsettest
    [2, 0] ... %diff
    );


p.Q_space = kron3(p.Q_space_local, p.M_space_local, p.M_space_local) ...
    + kron3(p.A_space_local, p.M_space_local, p.A_space_local') ...
    + kron3(p.A_space_local, p.A_space_local', p.M_space_local) ...
    + kron3(p.A_space_local', p.M_space_local, p.A_space_local) ...
    + kron3(p.M_space_local, p.M_space_local, p.Q_space_local) ...
    + kron3(p.M_space_local, p.A_space_local', p.A_space_local) ...
    + kron3(p.A_space_local', p.A_space_local, p.M_space_local) ...
    + kron3(p.M_space_local, p.A_space_local, p.A_space_local') ...
    + kron3(p.M_space_local, p.Q_space_local, p.M_space_local);

p.M_space = kron3(p.M_space_local, p.M_space_local, p.M_space_local);
p.A_space = kron3(p.M_space_local, p.M_space_local, p.A_space_local) ...
    + kron3(p.M_space_local, p.A_space_local, p.M_space_local) ...
    + kron3(p.A_space_local, p.M_space_local, p.M_space_local);
%% rhs
if length(p.f_time) ~= length(p.f_space_x) ...
        || length(p.f_time) ~= length(p.f_space_y) ...
        || length(p.f_time) ~= length(p.f_space_z)
    error('p.f_time and p.f_space should have the same dimensions')
end

rhs = [];

if ~isempty(p.f_time)
    
    ntest_time = length(p.T_time) - p.bSplineOrder_time;
    
    
    for i=1:length(p.f_time)
        rhs_time = zeros(ntest_time-(p.offset_time_test(1)+p.offset_time_test(2)),1);
        % time component
        
        
        for j=1+p.offset_time_test(1):ntest_time-p.offset_time_test(2)
            g = @(t) p.f_time{i}(t) * Ndiff(p.T_time, p.bSplineOrder_time,0,j,t);
            for step = 0:p.bSplineOrder_time-1
                rhs_time(j-p.offset_time_test(1)) = rhs_time(j-p.offset_time_test(1)) ...
                    + Gaussq(p.T_time(j+step),p.T_time(j+step+1),g,5);
            end
        end
        
        % first space component
        ntest_space = length(p.T_space) - p.bSplineOrder_space;
        rhs_space_x = zeros(ntest_space-(p.offset_space_test(1)+p.offset_space_test(2)),1);
        
        for j=1+p.offset_space_test(1):ntest_space-p.offset_space_test(2)
            g = @(x) p.f_space_x{i}(x) * Ndiff(p.T_space, p.bSplineOrder_space,0,j,x);
            for step = 0:p.bSplineOrder_space-1
                rhs_space_x(j-p.offset_space_test(1)) = rhs_space_x(j-p.offset_space_test(1)) ...
                    + Gaussq(p.T_space(j+step),p.T_space(j+step+1),g,5);
            end
        end
        
        % second space component
        rhs_space_y = zeros(ntest_space-(p.offset_space_test(1)+p.offset_space_test(2)),1);
        
        for j=1+p.offset_space_test(1):ntest_space-p.offset_space_test(2)
            g = @(y) p.f_space_y{i}(y) * Ndiff(p.T_space, p.bSplineOrder_space,0,j,y);
            for step = 0:p.bSplineOrder_space-1
                rhs_space_y(j-p.offset_space_test(1)) = rhs_space_y(j-p.offset_space_test(1)) ...
                    + Gaussq(p.T_space(j+step),p.T_space(j+step+1),g,5);
            end
        end
        
        % third space component
        rhs_space_z = zeros(ntest_space-(p.offset_space_test(1)+p.offset_space_test(2)),1);
        
        for j=1+p.offset_space_test(1):ntest_space-p.offset_space_test(2)
            g = @(z) p.f_space_z{i}(z) * Ndiff(p.T_space, p.bSplineOrder_space,0,j,z);
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
    
    
    %% Add the inital condtions
    for i=1+p.offset_time_test(1):length(p.T_time) - p.bSplineOrder_time - p.offset_time_test(2)
        
        if Ndiff(p.T_time, p.bSplineOrder_time,0,i,0) == 0 && ...
                Ndiff(p.T_time, p.bSplineOrder_time,1,i,0) == 0
            break
        end
        
        % first space component
        rhs_x_1 = sparse(ntest_space-(p.offset_space_test(1)+p.offset_space_test(2)),1);
        rhs_x_2 = sparse(ntest_space-(p.offset_space_test(1)+p.offset_space_test(2)),1);
        
        for j=1+p.offset_space_test(1):ntest_space-p.offset_space_test(2)
            g1 = @(x) p.u_1_x(x) * Ndiff(p.T_time, p.bSplineOrder_time,0,i,0) ...
                .* Ndiff(p.T_space, p.bSplineOrder_space,0,j,x);
            g2 = @(x) -p.u_0_x(x) * Ndiff(p.T_time, p.bSplineOrder_time,1,i,0) ...
                .* Ndiff(p.T_space, p.bSplineOrder_space,0,j,x);
            for step = 0:p.bSplineOrder_space-1
                rhs_x_1(j-p.offset_space_test(1)) = rhs_x_1(j-p.offset_space_test(1)) ...
                    + Gaussq(p.T_space(j+step),p.T_space(j+step+1),g1,5);
                rhs_x_2(j-p.offset_space_test(1)) = rhs_x_2(j-p.offset_space_test(1)) ...
                    + Gaussq(p.T_space(j+step),p.T_space(j+step+1),g2,5);
            end
        end
        
        % second space component
        rhs_y_1 = sparse(ntest_space-(p.offset_space_test(1)+p.offset_space_test(2)),1);
        rhs_y_2 = sparse(ntest_space-(p.offset_space_test(1)+p.offset_space_test(2)),1);
        
        for j=1+p.offset_space_test(1):ntest_space-p.offset_space_test(2)
            g1 = @(y) p.u_1_y(y) .* Ndiff(p.T_space, p.bSplineOrder_space,0,j,y);
            g2 = @(y) p.u_0_y(y) .* Ndiff(p.T_space, p.bSplineOrder_space,0,j,y);
            for step = 0:p.bSplineOrder_space-1
                rhs_y_1(j-p.offset_space_test(1)) = rhs_y_1(j-p.offset_space_test(1)) ...
                    + Gaussq(p.T_space(j+step),p.T_space(j+step+1),g1,5);
                rhs_y_2(j-p.offset_space_test(1)) = rhs_y_2(j-p.offset_space_test(1)) ...
                    + Gaussq(p.T_space(j+step),p.T_space(j+step+1),g2,5);
            end
        end
        
        % third space component
        rhs_z_1 = sparse(ntest_space-(p.offset_space_test(1)+p.offset_space_test(2)),1);
        rhs_z_2 = sparse(ntest_space-(p.offset_space_test(1)+p.offset_space_test(2)),1);
        
        for j=1+p.offset_space_test(1):ntest_space-p.offset_space_test(2)
            g1 = @(z) p.u_1_z(z) .* Ndiff(p.T_space, p.bSplineOrder_space,0,j,z);
            g2 = @(z) p.u_0_z(z) .* Ndiff(p.T_space, p.bSplineOrder_space,0,j,z);
            for step = 0:p.bSplineOrder_space-1
                rhs_z_1(j-p.offset_space_test(1)) = rhs_z_1(j-p.offset_space_test(1)) ...
                    + Gaussq(p.T_space(j+step),p.T_space(j+step+1),g1,5);
                rhs_z_2(j-p.offset_space_test(1)) = rhs_z_2(j-p.offset_space_test(1)) ...
                    + Gaussq(p.T_space(j+step),p.T_space(j+step+1),g2,5);
            end
        end
        
        p.rhs(:,i) = p.rhs(:,i) + kron(rhs_x_1, kron(rhs_y_1, rhs_z_1)) ...
            + kron(rhs_x_2, kron(rhs_y_2, rhs_z_2));
        
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