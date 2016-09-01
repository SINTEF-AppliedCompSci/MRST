classdef NonLinearSolverSLS < NonLinearSolver

   properties
      lambda_prev
      f_scaling 
      f_prev
      dx_prev
      dx_scaling
      H 
      %previousIncrement
   end
   
   methods
   
      function solver = NonLinearSolverSLS(varargin)
         solver = solver@NonLinearSolver(varargin{:});
         solver.lambda_prev = [];
         %solver.f_prev = 0; 
         solver.dx_prev = {};
         solver.dx_scaling = [];
         
         tmp = merge_options(struct('H', 0.03, 'f_scaling', []), varargin{:});
         solver.H = tmp.H;
         solver.f_scaling = tmp.f_scaling;
         
      end
      % ----------------------------------------------------------------------------
      function [state, report, ministates] = solveTimestep(solver, state0, dT, model, varargin)
         solver.resetSLSfields();
         [state, report, ministates] = solveTimestep@NonLinearSolver(solver, state0, dT, model, varargin{:});
      end
      % ----------------------------------------------------------------------------
      function resetSLSfields(solver)
         solver.lambda_prev = [];
         solver.dx_prev = {};
         solver.dx_scaling = [];
      end
      % ----------------------------------------------------------------------------
      function dx = stabilizeNewtonIncrements(solver, model, problem, dx)
         
         if isempty(solver.dx_scaling)
            % At first iteration, initialize f_scaling (hopefully, the
            % residuals are reasonable representative)
            for i = 1:numel(dx)
               if strcmpi(problem.primaryVariables{i}, 'pressure') ||...
                  strcmpi(problem.primaryVariables{i}, 'bhp')
                  solver.dx_scaling(i) = 1 ./ (5*max(problem.state.pressure));
               else
                  solver.dx_scaling(i) = 1;
               end
               %solver.dx_scaling(i) = 1./max(max(dx{i}), 1);
            end
         end
         
         if isempty(solver.dx_prev)
            solver.dx_prev = dx;
            for i = 1:numel(dx)
               solver.dx_prev{i} = solver.dx_prev{i} * 0;
            end
         end
         if isempty(solver.lambda_prev)
            solver.lambda_prev = ones(numel(dx), 1);
         end
         
                              
         %f = vertcat(dx{:}) .* solver.f_scaling;

         %% Using derivative-based controller
         % for i = 1:numel(dx)
         %    lambda(i) = sqrt( (2 * solver.H * solver.lambda_prev(i)) / ...
         %                      (norm(dx{i} - solver.dx_prev{i}) * solver.dx_scaling(i)));
         %    lambda(i) = min(lambda(i), 1);
         % end
         
         %% Using Richardson controller
         for i = 1:numel(dx)
            lambda(i) = solver.H / (norm(dx{i} - solver.dx_prev{i}) * solver.dx_scaling(i));
            lambda(i) = min(lambda(i), 1);
         end
         
         lambda
         %% updating internal state
         solver.dx_prev = dx;
         solver.lambda_prev = lambda;
         
         for i = 1:numel(dx)
            dx{i} = dx{i} * lambda(i);
         end
      end
      % % ----------------------------------------------------------------------------
      % function dx = stabilizeNewtonIncrements(solver, model, problem, dx)
         
      %    if isempty(solver.f_scaling)
      %       % At first iteration, initialize f_scaling (hopefully, the
      %       % residuals are reasonable representative)
      %       solver.f_scaling = solver.initFScaling(dx);
      %    end
                              
      %    f = vertcat(dx{:}) .* solver.f_scaling;

      %    %% Using derivative-based controller
      %    % lambda = sqrt( (2 * solver.H * solver.lambda_prev) / ...
      %    %                norm(f - solver.f_prev));
      %    % lambda = min(lambda,1);
         
      %    %% Using Richardson controller
      %    lambda = solver.H / norm(f - solver.f_prev);
      %    lambda = min(lambda, 1)
         
      %    %% updating internal state
      %    solver.f_prev = f;
      %    solver.lambda_prev = lambda;
         
      %    for i = 1:numel(dx)
      %       dx{i} = dx{i} * lambda;
      %    end
      % end
      % ----------------------------------------------------------------------------
      function fscaling = initFScaling(solver, dx)
         fscaling = rldecode(1./cellfun(@max, dx), cellfun(@numel, dx));
         fscaling(fscaling > 1) = 1; % @ necessary?
         fscaling(isinf(fscaling)) = 1;

      end
   end % methods
   
end
