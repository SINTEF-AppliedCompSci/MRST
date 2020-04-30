function buildLinearSolvers(varargin)
% Build linear solvers
    mrstModule add ad-core
    if mod(nargin, 2) == 1
        rebuild = varargin{1};
        varargin = varargin(2:end);
    else
        rebuild = false;
    end
    opt = struct('names', {{}});
    opt = merge_options(opt, varargin{:});
    if isempty(opt.names)
        names = {'amgcl_matlab', 'amgcl_matlab_block'};
    else
        names = opt.names;
        if ~iscell(names)
            names = {names};
        end
    end
    pth = fullfile(mrstPath('linearsolvers'), 'amgcl', 'utils');

    buildMexExtensions(rebuild, names, pth);
end

