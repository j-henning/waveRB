% Gets a problem configuration and creates a problem with all matrices, the
% right hand side and various variables for retrieving the solution
function [p] = create2DWaveProblemImproved(pC)


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

p.Q_space_local =StiffMat(p.T_space, ... % Tsol
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
p.A_space_local = StiffMat(p.T_space, ... % Tsol
    p.T_space, ... % Ttest
    p.bSplineOrder_space, ... %ksol
    p.bSplineOrder_space, ... %ktest
    p.offset_space_ansatz,... %offsetsol
    p.offset_space_test, ... %offsettest
    [2, 0] ... %diff
    );



p.Q_space =   p.c^4 *  (kron(p.Q_space_local, p.M_space_local) + kron(p.A_space_local, p.A_space_local') ...
    + kron(p.A_space_local', p.A_space_local) ...
    + kron(p.M_space_local, p.Q_space_local));

p.M_space = kron(p.M_space_local, p.M_space_local);
p.A_space =   p.c^2 *  ( -kron(p.A_space_local, p.M_space_local) - kron(p.M_space_local, p.A_space_local));

%% rhs
% if length(p.f_time) ~= length(p.f_space_x) || length(p.f_time) ~= length(p.f_space_y)
%     error('f_time and f_space should have the same dimensions')
% end

p.ntest_time = length(p.T_time) - p.bSplineOrder_time;
p.ntest_space = length(p.T_space) - p.bSplineOrder_space;

p.rhs = zeros((p.ntest_space-sum(p.offset_space_test))^2, ...
    p.ntest_time-sum(p.offset_time_test));

if ~isempty(p.f_time)
    
    
    
    %     for i=1:length(p.f_time)
    %
    %         % time component
    %
    %         rhs_time = zeros(p.ntest_time-(p.offset_time_test(1)+p.offset_time_test(2)),1);
    %         for j=1+p.offset_time_test(1):p.ntest_time-p.offset_time_test(2)
    %             g = @(t) p.f_time{i}(t) * Ndiff(p.T_time, p.bSplineOrder_time,0,j,(t==0)*(t+eps) + (t==1)*(t-eps) + (t > 0 && t < 1)*t);
    %             for step = 0:p.bSplineOrder_time-1
    %                 rhs_time(j-p.offset_time_test(1)) = rhs_time(j-p.offset_time_test(1)) ...
    %                     + Gaussq(p.T_time(j+step),p.T_time(j+step+1),g,5);
    %             end
    %         end
    %
    %         % first space component
    %         rhs_space_x = zeros(p.ntest_space-(p.offset_space_test(1)+p.offset_space_test(2)),1);
    %
    %         for j=1+p.offset_space_test(1):p.ntest_space-p.offset_space_test(2)
    %             g = @(x) p.f_space_x{i}(x) * Ndiff(p.T_space, p.bSplineOrder_space,0,j,(x==0)*(x+eps) + (x==1)*(x-eps) + (x > 0 && x < 1)*x);
    %             for step = 0:p.bSplineOrder_space-1
    %                 rhs_space_x(j-p.offset_space_test(1)) = rhs_space_x(j-p.offset_space_test(1)) ...
    %                     + Gaussq(p.T_space(j+step),p.T_space(j+step+1),g,5);
    %             end
    %         end
    %
    %         % second space component
    %         p.ntest_space = length(p.T_space) - p.bSplineOrder_space;
    %         rhs_space_y = zeros(p.ntest_space-(p.offset_space_test(1)+p.offset_space_test(2)),1);
    %
    %         for j=1+p.offset_space_test(1):p.ntest_space-p.offset_space_test(2)
    %             g = @(y) p.f_space_y{i}(y) * Ndiff(p.T_space, p.bSplineOrder_space,0,j,(y==0)*(y+eps) + (y==1)*(y-eps) + (y > 0 && y < 1)*y);
    %             for step = 0:p.bSplineOrder_space-1
    %                 rhs_space_y(j-p.offset_space_test(1)) = rhs_space_y(j-p.offset_space_test(1)) ...
    %                     + Gaussq(p.T_space(j+step),p.T_space(j+step+1),g,5);
    %             end
    %         end
    %
    %         if i == 1
    %             p.rhs = kron(rhs_time', kron(rhs_space_x, rhs_space_y));
    %         else
    %             p.rhs = p.rhs + kron(rhs_time', kron(rhs_space_x, rhs_space_y));
    %         end
    %     end
    %
    
    

%             %% Add the inital condtions
%             for i=1+p.offset_time_test(1):length(p.T_time) - p.bSplineOrder_time - p.offset_time_test(2)
%     
%                 if Ndiff(p.T_time, p.bSplineOrder_time,0,i,0) == 0 && ...
%                         Ndiff(p.T_time, p.bSplineOrder_time,1,i,0) == 0
%                     break
%                 end
%     
%                 % first space component
%                 rhs_x_1 = zeros(p.ntest_space-(p.offset_space_test(1)+p.offset_space_test(2)),1);
%                 rhs_x_2 = zeros(p.ntest_space-(p.offset_space_test(1)+p.offset_space_test(2)),1);
%     
%                 for j=1+p.offset_space_test(1):p.ntest_space-p.offset_space_test(2)
%                     g1 = @(x) p.u_1_x(x) * Ndiff(p.T_time, p.bSplineOrder_time,0,i,0) ...
%                         .* Ndiff(p.T_space, p.bSplineOrder_space,0,j,x);
%                     g2 = @(x) -p.u_0_x(x) * Ndiff(p.T_time, p.bSplineOrder_time,1,i,0) ...
%                         .* Ndiff(p.T_space, p.bSplineOrder_space,0,j,x);
%                     for step = 0:p.bSplineOrder_space-1
%                         rhs_x_1(j-p.offset_space_test(1)) = rhs_x_1(j-p.offset_space_test(1)) ...
%                             + Gaussq(p.T_space(j+step),p.T_space(j+step+1),g1,5);
%                         rhs_x_2(j-p.offset_space_test(1)) = rhs_x_2(j-p.offset_space_test(1)) ...
%                             + Gaussq(p.T_space(j+step),p.T_space(j+step+1),g2,5);
%                     end
%                 end
%     
%                 % second space component
%                 rhs_y_1 = zeros(p.ntest_space-(p.offset_space_test(1)+p.offset_space_test(2)),1);
%                 rhs_y_2 = zeros(p.ntest_space-(p.offset_space_test(1)+p.offset_space_test(2)),1);
%     
%                 for j=1+p.offset_space_test(1):p.ntest_space-p.offset_space_test(2)
%                     g1 = @(y) p.u_1_y(y) .* Ndiff(p.T_space, p.bSplineOrder_space,0,j,y);
%                     g2 = @(y) p.u_0_y(y) .* Ndiff(p.T_space, p.bSplineOrder_space,0,j,y);
%                     for step = 0:p.bSplineOrder_space-1
%                         rhs_y_1(j-p.offset_space_test(1)) = rhs_y_1(j-p.offset_space_test(1)) ...
%                             + Gaussq(p.T_space(j+step),p.T_space(j+step+1),g1,5);
%                         rhs_y_2(j-p.offset_space_test(1)) = rhs_y_2(j-p.offset_space_test(1)) ...
%                             + Gaussq(p.T_space(j+step),p.T_space(j+step+1),g2,5);
%                     end
%                 end
%     
%     %             p.rhs(:,i) = p.rhs(:,i) + kron(rhs_x_1, rhs_y_1) + kron(rhs_x_2, rhs_y_2);
%              p.rhs(:,i) = p.rhs(:,i) + kron(rhs_x_1, rhs_y_1) + kron(rhs_x_2, rhs_y_2);
%     
%             end
    
    % TODO: Reimplemente the rhs and u_1
    
    
    p.ntest_time = length(p.T_time) - p.bSplineOrder_time;
    p.ntest_space = length(p.T_space) - p.bSplineOrder_space;
    p.rhs = zeros((p.ntest_space-sum(p.offset_space_test))^2,p.ntest_time - sum(p.offset_time_test));
    
    
    for i=1+p.offset_time_test(1):length(p.T_time) - p.bSplineOrder_time - p.offset_time_test(2)
        v_dot0 = Ndiff(p.T_time, p.bSplineOrder_time,1,i,0);
        
        if v_dot0 == 0
            continue;
        end
        
        
        
        U0 = zeros(p.ntest_space-sum(p.offset_space_test), ...
            p.ntest_space-sum(p.offset_space_test));
        %     U1 = zeros(ntest-sum(test_offset),ntest-sum(test_offset));
        for j=1+p.offset_space_test(1):p.ntest_space-p.offset_space_test(2)
            for k=1+p.offset_space_test(1):p.ntest_space-p.offset_space_test(2)
                g = @(y)@(x) -v_dot0 .* p.u0(x,y) .* ...
                    Ndiff(p.T_space, p.bSplineOrder_space,0,j,x) ...
                    .* Ndiff(p.T_space, p.bSplineOrder_space,0,k,y); %*u0
                % h = @(y)@(x) u1(x,y) .* Ndiff(T_space_test, splineOrder,0,j,x) .* Ndiff(T_space_test, splineOrder,0,k,y); %*u1
                
                for stepX = 0:p.bSplineOrder_space-1
                    for stepY = 0:p.bSplineOrder_space-1
                        % Integrate from a to b
                        a = [p.T_space(j+stepX); p.T_space(k+stepY)];
                        b = [p.T_space(j+stepX+1); p.T_space(k+stepY+1)];
                        U0(j-p.offset_space_test(1), k-p.offset_space_test(1)) = ...
                            U0(j-p.offset_space_test(1), k-p.offset_space_test(1)) ...
                            + Gaussq(a,b,g,5);
                        %                         U1(j-test_offset(1), k-test_offset(1)) = U1(j-test_offset(1), k-test_offset(1)) ...
                        %                             + Gaussq(a,b,h,5);
                    end
                end
            end
        end
        p.rhs(:,i) = p.rhs(:,i) + reshape(U0', [numel(U0) 1]); % Transpose?
    end
    
    
end




% Calculate some additional values
p.nsol_time = length(p.T_time) - p.bSplineOrder_time;
p.nsol_space = length(p.T_space) - p.bSplineOrder_space;
p.dim_time =  p.nsol_time-(p.offset_time_ansatz(1)+p.offset_time_ansatz(2));
end