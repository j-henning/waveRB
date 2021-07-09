% Gets a problem configuration and creates a problem with all matrices, the
% right hand side and various variables for retrieving the solution
function [p] = create3DWaveProblemNonOptimalImproved(pC)


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
p.A_space = p.c^2 * (kron3(p.A_space_local, p.M_space_local, p.M_space_local) ...
    + kron3(p.M_space_local, p.A_space_local, p.M_space_local) ...
    + kron3(p.M_space_local, p.M_space_local, p.A_space_local));







%% rhs

if ~isempty(p.f_time)
    
    %     ntest_time = length(p.T_time_test) - p.bSplineOrder_time_test;
    %     for i=1:length(p.f_time)
    %         rhs_time = zeros(ntest_time-(p.offset_time_test(1)+p.offset_time_test(2)),1);
    %         % time component
    %
    %
    %         for j=1+p.offset_time_test(1):ntest_time-p.offset_time_test(2)
    %             g = @(t) p.f_time{i}(t) * Ndiff(p.T_time_test, p.bSplineOrder_time_test,0,j,t);
    %             for step = 0:p.bSplineOrder_time_test-1
    %                 rhs_time(j-p.offset_time_test(1)) = rhs_time(j-p.offset_time_test(1)) ...
    %                     + Gaussq(p.T_time_test(j+step),p.T_time_test(j+step+1),g,5);
    %             end
    %         end
    %
    %         % first space component
    %         ntest_space = length(p.T_space_test) - p.bSplineOrder_space_test;
    %         rhs_space_x = zeros(ntest_space-(p.offset_space_test(1)+p.offset_space_test(2)),1);
    %
    %         for j=1+p.offset_space_test(1):ntest_space-p.offset_space_test(2)
    %             g = @(x) p.f_space_x{i}(x) * Ndiff(p.T_space_test, p.bSplineOrder_space_test,0,j,x);
    %             for step = 0:p.bSplineOrder_space_test-1
    %                 rhs_space_x(j-p.offset_space_test(1)) = rhs_space_x(j-p.offset_space_test(1)) ...
    %                     + Gaussq(p.T_space_test(j+step),p.T_space_test(j+step+1),g,5);
    %             end
    %         end
    %
    %         % second space component
    %         rhs_space_y = zeros(ntest_space-(p.offset_space_test(1)+p.offset_space_test(2)),1);
    %
    %         for j=1+p.offset_space_test(1):ntest_space-p.offset_space_test(2)
    %             g = @(y) p.f_space_y{i}(y) * Ndiff(p.T_space_test, p.bSplineOrder_space_test,0,j,y);
    %             for step = 0:p.bSplineOrder_space_test-1
    %                 rhs_space_y(j-p.offset_space_test(1)) = rhs_space_y(j-p.offset_space_test(1)) ...
    %                     + Gaussq(p.T_space_test(j+step),p.T_space_test(j+step+1),g,5);
    %             end
    %         end
    %
    %         % third space component
    %         rhs_space_z = zeros(ntest_space-(p.offset_space_test(1)+p.offset_space_test(2)),1);
    %
    %         for j=1+p.offset_space_test(1):ntest_space-p.offset_space_test(2)
    %             g = @(z) p.f_space_z{i}(z) * Ndiff(p.T_space_test, p.bSplineOrder_space_test,0,j,z);
    %             for step = 0:p.bSplineOrder_space_test-1
    %                 rhs_space_z(j-p.offset_space_test(1)) = rhs_space_z(j-p.offset_space_test(1)) ...
    %                     + Gaussq(p.T_space_test(j+step),p.T_space_test(j+step+1),g,5);
    %             end
    %         end
    %
    %         if i == 1
    %             p.rhs = kron(rhs_time', kron(rhs_space_z, kron(rhs_space_y, rhs_space_x)));
    %         else
    %             p.rhs = p.rhs + kron(rhs_time', kron(rhs_space_z, kron(rhs_space_y, rhs_space_x)));
    %         end
    %     end
    
    p.ntest_time = length(p.T_time_test) - p.bSplineOrder_time_test;
    p.ntest_space = length(p.T_space_test) - p.bSplineOrder_space_test;
    p.rhs = zeros((p.ntest_space-sum(p.offset_space_test))^3,p.ntest_time - sum(p.offset_time_test));
    
    U0 = zeros(p.ntest_space-sum(p.offset_space_test), ...
        p.ntest_space-sum(p.offset_space_test), ...
        p.ntest_space-sum(p.offset_space_test));
    
    
    %     U1 = zeros(ntest-sum(test_offset),ntest-sum(test_offset));
    
    
    
    for j=1+p.offset_space_test(1):p.ntest_space-p.offset_space_test(2)
        for k=1+p.offset_space_test(1):p.ntest_space-p.offset_space_test(2)
            for l=1+p.offset_space_test(1):p.ntest_space-p.offset_space_test(2)
                
                
                g = @(z)@(y)@(x) p.u0(x,y,z) .* ...
                    Ndiff(p.T_space_test, p.bSplineOrder_space_test,0,j,x) ...
                    .* Ndiff(p.T_space_test, p.bSplineOrder_space_test,0,k,y) ...
                    .* Ndiff(p.T_space_test, p.bSplineOrder_space_test,0,l,z); %*u0
                
                %                 g2 = @(x,y,z) p.u0(x,y,z) .* ...
                %                     Ndiff(p.T_space_test, p.bSplineOrder_space_test,0,j,x) ...
                %                     .* Ndiff(p.T_space_test, p.bSplineOrder_space_test,0,k,y) ...
                %                     .* Ndiff(p.T_space_test, p.bSplineOrder_space_test,0,l,z);
                % h = @(y)@(x) u1(x,y) .* Ndiff(T_space_test, splineOrder,0,j,x) .* Ndiff(T_space_test, splineOrder,0,k,y); %*u1
                
                % Integrate
                for stepX = 0:p.bSplineOrder_space_test-1
                    for stepY = 0:p.bSplineOrder_space_test-1
                        for stepZ = 0:p.bSplineOrder_space_test-1
                            
                            a = [p.T_space_test(j+stepX), p.T_space_test(k+stepY), p.T_space_test(l+stepZ)];
                            b = [p.T_space_test(j+stepX+1), p.T_space_test(k+stepY+1), p.T_space_test(l+stepZ+1)];
                            
                            
                            % Check if u0 is zero in the center of the
                            % inteval, as well as on the boundary
                            % if yes assume that the integral is zero.
                            % This is done for performance reasons
                            center = 0.5 * (a+b);
                            if abs(p.u0(center(1),center(2),center(3))) < eps && ...
                                    abs(p.u0(a(1),a(2),a(3))) < eps && ...
                                    abs(p.u0(a(1),a(2),b(3))) < eps && ...
                                    abs(p.u0(a(1),b(2),a(3))) < eps && ...
                                    abs(p.u0(a(1),b(2),b(3))) < eps && ...
                                    abs(p.u0(b(1),b(2),b(3))) < eps && ...
                                    abs(p.u0(b(1),b(2),a(3))) < eps && ...
                                    abs(p.u0(b(1),a(2),b(3))) < eps && ...
                                    abs(p.u0(b(1),a(2),a(3))) < eps
                                continue;
                            end
                            
                            U0(j-p.offset_space_test(1), ...
                                k-p.offset_space_test(1), ...
                                l-p.offset_space_test(1)) = ...
                             U0(j-p.offset_space_test(1), ...
                                k-p.offset_space_test(1), ...
                                l-p.offset_space_test(1)) ...
                             + Gaussq(a,b,g,3);
                            %                                intVal;
                            %
                            
                            %                         U1(j-test_offset(1), k-test_offset(1)) = U1(j-test_offset(1), k-test_offset(1)) ...
                            %                             + Gaussq(a,b,h,5);
                        end
                    end
                end
                
            end
        end
    end
    for i=1+p.offset_time_test(1):length(p.T_time_test) - p.bSplineOrder_time_test - p.offset_time_test(2)
        v_dot0 = Ndiff(p.T_time_test, p.bSplineOrder_time_test,1,i,0);
        
        if v_dot0 == 0
            continue;
        end
        
        p.rhs(:,i) = p.rhs(:,i) + reshape(-v_dot0 * U0, [numel(U0) 1]);
    end
    
    % Calculate some additional values
    p.nsol_time = length(p.T_time_ansatz) - p.bSplineOrder_time_ansatz;
    p.nsol_space = length(p.T_space_ansatz) - p.bSplineOrder_space_ansatz;
    p.dim_time =  p.nsol_time-(p.offset_time_ansatz(1)+p.offset_time_ansatz(2));
end
end

function D = kron3(A,B,C)
D = kron(A,kron(B,C));
end