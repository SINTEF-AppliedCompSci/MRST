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
         solver.f_prev = 0; 
         solver.dx_prev = {};
         solver.dx_scaling = [];
         
         tmp = merge_options(struct('H', 0.2, 'f_scaling', []), varargin{:});
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
         solver.f_scaling = [];
         solver.lambda_prev = [];
         solver.dx_prev = {};
         solver.dx_scaling = [];
         solver.f_prev = 0;
      end
      % ----------------------------------------------------------------------------
      % function dx = stabilizeNewtonIncrements(solver, model, problem, dx)
         
      %    if isempty(solver.dx_scaling)
      %       % At first iteration, initialize f_scaling (hopefully, the
      %       % residuals are reasonable representative)
      %       for i = 1:numel(dx)
      %          if strcmpi(problem.primaryVariables{i}, 'pressure') ||...
      %             strcmpi(problem.primaryVariables{i}, 'bhp')
      %             solver.dx_scaling(i) = 1 ./ (1*max(problem.state.pressure));
      %          else
      %             solver.dx_scaling(i) = 1;
      %          end
      %          %solver.dx_scaling(i) = 1./max(max(dx{i}), 1);
      %       end
      %    end
         
      %    if isempty(solver.dx_prev)
      %       solver.dx_prev = dx;
      %       for i = 1:numel(dx)
      %          solver.dx_prev{i} = solver.dx_prev{i} * 0;
      %       end
      %    end
      %    if isempty(solver.lambda_prev)
      %       solver.lambda_prev = ones(numel(dx), 1);
      %    end
         
                              
      %    %f = vertcat(dx{:}) .* solver.f_scaling;

      %    %% Using derivative-based controller
      %    % for i = 1:numel(dx)
      %    %    lambda(i) = sqrt( (2 * solver.H * solver.lambda_prev(i)) / ...
      %    %                      (norm(dx{i} - solver.dx_prev{i}) * solver.dx_scaling(i)));
      %    %    lambda(i) = min(lambda(i), 1);
      %    % end
         
      %    % %% Using Odd's modified derivative-based controller
      %    % for i = 1:numel(dx)
      %    %    f      = dx{i} * solver.dx_scaling(i);
      %    %    f_prev = solver.dx_prev{i} * solver.dx_scaling(i);
      %    %    fac = 2 * solver.H * norm(f) / norm(f-f_prev);
      %    %    lambda(i) = solver.lambda_prev(i) * fac;
      %    %    lambda(i) = min(lambda(i), 1);
      %    % end
         
         
      %    %% Using Richardson controller
      %    for i = 1:numel(dx)
      %       f1 = dx{i} / norm(dx{i});
      %       f2 = solver.dx_prev{i} / norm(solver.dx_prev{i});
      %       lambda(i) = solver.H / (norm(f1 - f2));
      %       %lambda(i) = solver.H / (norm(dx{i} - solver.dx_prev{i}) * solver.dx_scaling(i));
      %       lambda(i) = min(lambda(i), 1);
      %    end
         
      %    lambda
      %    %% updating internal state
      %    solver.dx_prev = dx;
      %    solver.lambda_prev = lambda;
         
      %    for i = 1:numel(dx)
      %       dx{i} = dx{i} * lambda(i);
      %    end
      % end
      % ----------------------------------------------------------------------------
      function dx = stabilizeNewtonIncrements(solver, model, problem, dx)
         solver.H = 0.5; % @@
         if isempty(solver.f_scaling)
            % At first iteration, initialize f_scaling (hopefully, the
            % residuals are reasonable representative)
            solver.f_scaling = solver.initFScaling(problem, dx);
         end
         if isempty(solver.lambda_prev)
            solver.lambda_prev = 1;
         end                              

         f = vertcat(dx{:}) .* solver.f_scaling;

         %% Using derivative-based controller
         % lambda = sqrt( (2 * solver.H * solver.lambda_prev) / ...
         %                norm(f - solver.f_prev));
         % lambda = min(lambda,1);
         
         %% Using Odd's modified derivative-based controller
         % DOESN'T WORK
         % lambda = solver.lambda_prev * 2 * solver.H * norm(f) / norm(f-solver.f_prev);
         % lambda = min(lambda,1)
         
         %% Using Richardson controller
         % SEEMS TO WORK WELL WITH THE NEW SCALING I INTRODUCED (NOT
         % REFLECTED IN PAPER)
         f1 = f / norm(f);
         f2 = solver.f_prev / norm(solver.f_prev);
         lambda = solver.H / norm(f1 - f2);
         %lambda = solver.H / norm(f - solver.f_prev)
         lambda = min(lambda, 1)
         
         %% updating internal state
         solver.f_prev = f;
         solver.lambda_prev = lambda;
         
         for i = 1:numel(dx)
            dx{i} = dx{i} * lambda;
         end
      end
      % ----------------------------------------------------------------------------
      function fscaling = initFScaling(solver, problem, dx)

         tmp = dx;
         for i = 1:numel(dx)
            if strcmpi(problem.primaryVariables{i}, 'pressure') || ...
                       strcmpi(problem.primaryVariables{i}, 'bhp')
               tmp{i} = tmp{i} * 0 + 1 ./ (2*max(problem.state.pressure));
            else
               tmp{i} = tmp{i} * 0 + 1;
            end
         end
         fscaling = vertcat(tmp{:});
         %fscaling = rldecode(1./(2*cellfun(@max, dx)), cellfun(@numel, dx));
         %fscaling(fscaling > 1) = 1; % @ necessary?
         %fscaling(isinf(fscaling)) = 1;
      end
   end % methods
   
end
