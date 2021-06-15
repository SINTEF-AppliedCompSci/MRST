function [suite, names] = getExampleIntegrationTestSuiteMRST(modules, varargin)
%Undocumented Utility Function

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

    opt = struct('seperateModules', false);
    opt = merge_options(opt, varargin{:});
    
    mrstModule add ad-unittest
    import matlab.unittest.TestSuite;
    suite = TestSuite.fromClass(?MRSTExampleTests);
    if nargin > 0
        if ~iscell(modules)
            modules = {modules};
        end
        suite = filter_module(suite, modules);
    end
    
    if opt.seperateModules
        suite0 = suite;
        mods = mrstPath();
        nm = numel(mods);
        suite = cell(nm, 1);
        names = cell(nm, 1);
        for i = 1:nm
            suite{i} = filter_module(suite0, mods(i));
            names{i} = lower(mods{i});
        end
        keep = not(cellfun(@isempty, suite));
        suite = suite(keep);
        names = names(keep);
    else
        names = {'examples'};
    end
end

function suite = filter_module(suite, modules)
        active = arrayfun(@(x) any(strcmpi(x.Parameterization(2).Value, modules)),  suite);
        suite = suite(active);

end
