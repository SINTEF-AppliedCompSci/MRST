function problem = transformProblem(problem, varargin)
    % Apply linear transformation to linearized problem
    % Transforms problem.equations{i} as
    %
    %  - ML*A*MR*x   = ML*b  + B if size(A) = [nc, nc],
    %  - ML*A*MRw*x  = ML*b      if size(A) = [nc, nw],
    %  - MLw*A*MR*x  = MLw*b     if size(A) = [nw, nc],
    %  - MLw*A*MRw*x = MLw*b     if size(A) = [nw, nw],
    %
    % where A*x = b is the linearized form of the equation i, with respect
    % to variable set j, nc is the number of grid cells, and nw the number
    % of wells. Empty input matrices are initialized to 1. B is a function
    % i and j, initilized to 0 if empty.

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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

    opt = struct('ML', [], 'MR', [], 'MLw', [], 'MRw', [], 'B', []);
    opt = merge_options(opt, varargin{:});
    opt = checkInputs(opt);
    % Get subset and problem sizes
    is_cell = strcmpi(problem.types, 'cell');
    ndof = cellfun(@(eq) numel(value(eq)), problem.equations);
    nc   = ndof(is_cell); nc = nc(1);
    ML   = opt.ML;
    MLw  = opt.MLw;
    % Restrict equations
    isAD = isa(problem.equations{1}, 'ADI');
    for i = 1:problem.numeq
        eq = problem.equations{i};
        if ~isAD
            n = numel(eq);
            if n == nc
                % Cell equation (nc)
                eq = ML*eq;
            else
                % Well equation (nw)
                eq = MLw*eq;
            end
        else
            eq = transformEquations(eq, i, nc, opt);
        end
        problem.equations{i} = eq;
    end
end
        
%-------------------------------------------------------------------------%
function eq = transformEquations(eq, i, nc, opt)
    ML  = opt.ML;
    MR  = opt.MR;
    MLw = opt.MLw;
    MRw = opt.MRw;
    B   = opt.B;
    % Get equation and Jacobian block sizes
    n = numel(value(eq));
    m = cellfun(@(jac) sizeJac(jac, 2), eq.jac);
    % Restrict equation value
    if n == nc
        % Cell equation (nc)
        eq.val = ML*eq.val;
    else
        % Well equation (nw)
        eq.val = MLw*eq.val;
    end
    % Restrict equation Jacobians
    for j = 1:numel(eq.jac)
        if n == nc && m(j) >= nc
            % nc x nc
            eq.jac{j} = ML*eq.jac{j}*MR + B(i,j);
        elseif n == nc
            % nc x nw
            eq.jac{j} = ML*eq.jac{j}*MRw;
        elseif m(j) >= nc
            % nw x nc
            eq.jac{j} = MLw*eq.jac{j}*MR;
        elseif n < nc && m(j) < nc
            % nw x nw
            eq.jac{j} = MLw*eq.jac{j}*MRw;
        else
            error('Inconsistent Jacobian size!')
        end
    end
end

%-------------------------------------------------------------------------%
function opt = checkInputs(opt)
    matrices = {'ML', 'MR', 'MLw', 'MRw'};
    for m = matrices
        if isempty(opt.(m{1}))
            opt.(m{1}) = 1;
        end
    end
    if isempty(opt.B)
        opt.B = @(i,j) sparse(1,1);
    end
end
