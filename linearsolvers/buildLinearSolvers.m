function buildLinearSolvers(varargin)
% Build linear solvers

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
