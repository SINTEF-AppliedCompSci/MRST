function [suite, names] = getExampleIntegrationTestSuiteMRST(modules, varargin)
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