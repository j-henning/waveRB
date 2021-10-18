# waveRB
A Reduced Basis Method for the wave equation in the very weak space-time setting. The B-Spline implementation is taken from C. Mollet's PhD Thesis *Parabolic PDEs in Space-Time Formulations: Stability for Petrov-Galerkin Discretizations with B-Splines and Existence of Moments for Problems with Random Coefficients* which can be found [here](https://kups.ub.uni-koeln.de/6872/)

## Quick start
First of all, it is always a good idea to add the whole source folder with the command

    addpath(genpath('../source'))
	
This recursively adds all folders which lie inside the source folder.

### Solving the Wave Equation	
We want to solve a one dimensional wave problem. After we added the source folder as described above, we can define our problem:

    problemConfiguration = defineProblem(1, ... % Dimension
        {@(t) t.^2}, ... % f_time
        {@(x) sin(pi*x)}, ... % f_space
        @(x) sin(3*pi*x), ... % u_0
        @(x) x .* (x-1), ... % u_1
        1, ... % Wave speed
        5, ... % Time refinement level
        5); % Space refinement level
Defining a problem only creates a struct where all the needed parameters are stored. To assemble the matrices and the right hand side, we have to invoke 

    problem = createProblem(problemConfiguration);

After that, the linear system can be solved via

    U = solveProblem(problem);	

If we want to evaluate the solution, we first have to specify the resolution in space and time by


    resolution.x = 10; % Space resolution, 2^10 + 1 points in space
    resolution.t = 10; % Time resolution, 2^10 + 1 points in time

and then get the solution by

    solution = getSolution(problem, U, resolution);

If we want to plot the solution, we can use

    plotSolution(solution, resolution); % Surf plot 

The whole code also works with two and three dimensional problems. We just need the adapt the call to `defineProblem` and specify `resolution.y` and `resolution.z`

### Creating a Reduced Basis model
Creating a Reduced Basis model takes a bit more lines of code, but is still easy to create (please refer to the documentation of each method for the input and output parameters)

    problemConfiguration = defineProblem(...); 
    pOne = createProblem(problemConfiguration);
    
    [tree, treeEstimator, pOne, splines] = rbOffline(problemConfiguration, ...
        resolution, muMin, muMax, Nh, Np, hTolerance, pTolerance, height,...
        NhEst, NpEst,hToleranceEst, pToleranceEst, heightEst, Xi, ...
        solver, maxIt, tolerance1, tolerance2);
        
    matrices = {kron(pOne.Q_time, pOne.M_space); ...
        kron(pOne.D_time, pOne.A_space') + kron(pOne.D_time', pOne.A_space); ...
        kron(pOne.M_time, pOne.Q_space)};
    coefficients = {@(mu) 1; @(mu) mu; @(mu) mu.^2};  
    
    tree = computeReducedMatrices(tree, pOne, matrices, coefficients);
    treeEstimator = computeReducedMatrices(treeEstimator, pOne, matrices, coefficients);
    
    theta = thetaCalculation(problemConfiguration, muMin, muMax, 1000, ...
        100, 0.05, tree, treeEstimator, pOne, splines, resolution, solver, ...
        maxIt, tolerance1, tolerance2); 
		
The online stage could look like this

    mu_test = 0.712; 
    u_N_rec = getRBhpSolutionVector(tree, mu_test);
    u_M_rec = getRBhpSolutionVector(treeEstimator, mu_test);
        
    sol_rec_N = getSolution(p, u_N_rec,  resolution, splines);    
    sol_rec_M = getSolution(p, u_M_rec,  resolution, splines);
        
    % Calculate a truth solution
    p = changeWaveSpeed(pOne, mu_test);
    U_test = solveProblem(p, solver, maxIt, tolerance1, tolerance2);
    sol_test{j} = getSolution(p, Utest,  resolution, splines);
    
    % Calculate the errors
    error_test = sqrt(mean((solTest{j}-sol_rec_N).^2, [1 2]));
    errorEst_test = sqrt(mean((sol_rec_M-sol_rec_N).^2, [1 2])) / (1 - theta);    		