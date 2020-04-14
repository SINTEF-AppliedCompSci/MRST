function problem = transformProblem(problem, varargin)
    opt = struct('ML', [], 'MR', [], 'MLw', [], 'MRw', [], 'B', [], 'zero', []);
    opt = merge_options(opt, varargin{:});
    opt = checkInputs(opt);
    % Get subset and problem sizes
    is_cell = strcmpi(problem.types, 'cell');
    ndof = cellfun(@(eq) numel(value(eq)), problem.equations);
    nc   = ndof(is_cell); nc = nc(1);
    n    = numel(value(problem.equations{1}));
    ML   = opt.ML;
    MLw  = opt.MLw;
    % Restrict equations
    isAD = isa(problem.equations{1}, 'ADI');
    for i = 1:numel(problem)
        eq = problem.equations{i};
        if ~isAD
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