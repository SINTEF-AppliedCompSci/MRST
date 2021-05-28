classdef AMGCLSolverAD < LinearSolverAD
    % Linear solver that calls external compiled multigrid solver
    %
    % SYNOPSIS:
    %   solver = AMGCLSolverAD()
    %
    % DESCRIPTION:
    %    AD-interface for the AMGCL interface.
    %
    % NOTE:
    %    This solver requires AMGCL to be installed and working.
    %
    % SEE ALSO:
    %   `BackslashSolverAD`

   properties
       amgcl_setup
       reuseMode = 1; % 1 for no reuse, 2 for re-use
   end

   methods
       function solver = AMGCLSolverAD(varargin)
            require linearsolvers
            solver = solver@LinearSolverAD();
            solver.reduceToCell = true;
            solver.tolerance = 1e-6;
            solver.maxIterations = 100;
            
            [solver, extra] = merge_options(solver, varargin{:});
            solver.amgcl_setup = getAMGCLMexStruct(extra{:});
       end
       
       function [result, report] = solveLinearSystem(solver, A, b, varargin)
           [result, report] = solver.callAMGCL_MEX(A, b, 1, varargin{:});
       end
       
       function setCoarsening(solver, varargin)
           solver.setParameterGroup('coarsening', [], varargin{:});
       end
       
       function setRelaxation(solver, varargin)
           solver.setParameterGroup('relaxation', [], varargin{:});
       end
       
       function setPreconditioner(solver, varargin)
           solver.setParameterGroup('preconditioner', [], varargin{:});
       end
       
       function setSolver(solver, varargin)
           solver.setParameterGroup('solver', [], varargin{:});
       end
       
       function setParameterGroup(solver, group, fld, v)
           group = lower(group);
           if isempty(fld)
               fld = group;
           end
           if nargin == 3
               fprintf('No %s argument given. Available options:\n', group)
               [~, choices, descriptions] = translateOptionsAMGCL(group, []);
               for i = 1:numel(choices)
                   fprintf('%10s: ', choices{i});
                   disp(descriptions{i});
               end
           else
               solver.amgcl_setup.(fld) = translateOptionsAMGCL(group, v);
           end
       end
       
       function [n, name, descr, parameters] = getParameterGroup(solver, group, fld)
           group = lower(group);
           if nargin < 3 || isempty(fld)
               fld = group;
           end
           [n, name, descr, parameters] = translateOptionsAMGCL(group, solver.amgcl_setup.(fld));
       end
       
       function [result, report] = callAMGCL_MEX(solver, A, b, id, varargin)
            timer = tic();
            if ~isempty(varargin) && isempty(varargin{1})
                varargin = {};
            end
            [result, res, its] = amgcl_matlab(A', b, solver.amgcl_setup, solver.tolerance, solver.maxIterations, id, solver.reuseMode, varargin{:});
            t_solve = toc(timer);
            if res > solver.tolerance
                warning(['Solver did not converge to specified tolerance of %1.3e in %d iterations. ', ...
                    'Reported residual estimate was %1.3e after %2.2f seconds'], solver.tolerance, its, res, t_solve);
            elseif solver.verbose
                fprintf('AMGCL solver converged to %1.3e in %2d iterations after %2.2f seconds.\n', res, its, t_solve);
            end
            report = solver.getSolveReport(...
                            'Converged',  res <= solver.tolerance, ...
                            'Residual',   res,...
                            'Iterations', its);
       end
       
        function  solver = cleanupSolver(solver, A, b, varargin) %#ok
            if solver.reuseMode > 1
                resetAMGCL();
            end
        end
        
        function [d, sn] = getDescription(solver)
            sn = 'AMGCL';
            if solver.amgcl_setup.block_size > 1
                sn = [sn, '-block'];
            end
            sn = [sn, solver.id];

            prm = {'solver', 'preconditioner', 'relaxation'};
            if solver.amgcl_setup.preconditioner == 1
                prm{end+1} = 'coarsening';
            end
            tmp = cell(1, numel(prm));
            for i = 1:numel(prm)
                s = prm{i};
                tmp{i} = solver.getFeatureDescription(s);
            end
            d = [sprintf('General AMGCL solver. Configuration:\n'), ...
                 sprintf('\t%s\n', tmp{:})];
        end
        
        function out = getFeatureDescription(solver, name, extra)
            if nargin == 2
                extra = name;
            end
            [ix, choice, description, param] = solver.getParameterGroup(name, extra); %#ok
            out = sprintf('%15s: %s', extra, choice);
            out = [out, sprintf(' (%s)', description)];
            if ~isempty(param)
                s = '';
                for j = 1:numel(param)
                     p = param{j};
                     s = [s, sprintf('\n\t\t\t\t\t - %s = %g', p, solver.amgcl_setup.(p))];
                end
                out = [out, s];
            end
        end
   end
end

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
