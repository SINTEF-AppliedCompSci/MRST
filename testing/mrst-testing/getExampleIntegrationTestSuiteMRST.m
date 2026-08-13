function [suite, names] = getExampleIntegrationTestSuiteMRST(modules, varargin)
%Get test suite of MRST example smoke tests.
%
% SYNOPSIS:
%   [suite, names] = getExampleIntegrationTestSuiteMRST()
%   [suite, names] = getExampleIntegrationTestSuiteMRST(modules)
%   [suite, names] = getExampleIntegrationTestSuiteMRST(modules, 'key', val)
%
% OPTIONAL PARAMETERS:
%   modules         - Cell array of module names to include. If omitted,
%                     all registered modules are included.
%   'seperateModules' - Return one suite per module (logical, default false).
%   'excludeTags'   - Cell array of tags to exclude. Examples annotated
%                     with any of these tags are omitted from the suite.
%                     Common values: 'slow', 'data-required', 'interactive'.
%                     Default: {'slow', 'data-required'}.
%
% RETURNS:
%   suite - TestSuite object (or cell array when seperateModules is true).
%   names - Cell array of suite names.
%
% EXAMPLE:
%   % Default: excludes slow and data-required examples
%   suite = getExampleIntegrationTestSuiteMRST();
%
%   % Include slow examples too (e.g. for a nightly run)
%   suite = getExampleIntegrationTestSuiteMRST([], 'excludeTags', {});
%
%   % Only data-required filter
%   suite = getExampleIntegrationTestSuiteMRST([], 'excludeTags', {'data-required'});
%
% SEE ALSO:
%   MRSTExampleTests, runPreReleaseTests

%{
Copyright 2009-2026 SINTEF Digital, Mathematics & Cybernetics.

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

    opt = struct('seperateModules', false, ...
                 'excludeTags'    , {{'slow', 'data-required'}});
    opt = merge_options(opt, varargin{:});
    
    mrstModule add mrst-testing

    suite = matlab.unittest.TestSuite.fromClass(?MRSTExampleTests);

    if nargin > 0 && ~isempty(modules)
        if ~iscell(modules)
            modules = {modules};
        end
        suite = filter_module(suite, modules);
    end

    % Filter by tags: remove examples that carry any of the excluded tags
    if ~isempty(opt.excludeTags)
        suite = filter_by_tags(suite, opt.excludeTags);
    end

    if opt.seperateModules
        suite0 = suite;
        mods = mrstPath();
        nm   = numel(mods);
        suite = cell(nm, 1);
        names = cell(nm, 1);
        for i = 1:nm
            suite{i} = filter_module(suite0, mods(i));
            names{i} = lower(mods{i});
        end
        keep  = ~cellfun(@isempty, suite);
        suite = suite(keep);
        names = names(keep);
    else
        names = {'examples'};
    end
end

% -------------------------------------------------------------------------

function suite = filter_module(suite, modules)
    active = arrayfun( ...
        @(x) any(strcmpi(x.Parameterization(2).Value, modules)), suite);
    suite = suite(active);
end

% -------------------------------------------------------------------------

function suite = filter_by_tags(suite, excludeTags)
    % Remove test cases whose example file carries an excluded tag.
    keep = true(numel(suite), 1);
    for i = 1:numel(suite)
        exName = suite(i).Parameterization(1).Value;
        opt    = mrstExampleOptions(exName);
        if any(ismember(lower(opt.tags), lower(excludeTags)))
            keep(i) = false;
        end
    end
    suite = suite(keep);
end
