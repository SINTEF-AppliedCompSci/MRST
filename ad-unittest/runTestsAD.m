function results = runTestsAD(varargin)
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

    opt = struct(...
        'runUnit', true, ...
        'runIntegration', true, ...
        'runExamples', true, ...
        'writeXML',   false, ...
        'writeToDisk',  false ...
        );
    opt = merge_options(opt, varargin{:});
    
    import matlab.unittest.TestRunner;
    if opt.writeXML
        import matlab.unittest.plugins.XMLPlugin;
        out_xml = fullfile(mrstPath('query', 'ad-unittest'), 'output', 'XML');
        if exist(out_xml, 'dir') == 0
            mkdir(out_xml);
        end
    end
    ix = 1;
    [names, suites] = deal(cell(opt.runUnit + opt.runIntegration, 1));
    if opt.runUnit
        suites{ix} = getUnitTestSuiteMRST();
        names{ix} = 'unit';
        ix = ix+1;
    end
    
    if opt.runIntegration
        suites{ix} = getIntegrationTestSuiteMRST();
        names{ix} = 'integration';
        ix = ix+1;
    end
    
    if opt.runExamples
        [suite, nm] = getExampleIntegrationTestSuiteMRST(mrstPath(), 'seperateModules', true);
        suites = [suites; suite];
        names = [names; nm];
    end
    
    results = [];
    for i = 1:numel(suites)
        fprintf('Running test suite %d of %d: %s\n', i, numel(suites), names{i});
        runner = TestRunner.withTextOutput;
        if opt.writeXML
            jFile = fullfile(out_xml, [names{i}, '.xml']);
            p = XMLPlugin.producingJUnitFormat(jFile);
            runner.addPlugin(p);
        end
        res = runner.run(suites{i});
        results = [results, res]; %#ok
        if opt.writeToDisk
            outd = fullfile(mrstPath('query', 'ad-unittest'), 'output', 'TAP');
            if exist(outd, 'dir') == 0
                mkdir(outd);
            end
            mrstModule add ad-unittest
            tapFile = fullfile(outd, [names{i}, '.tap']);
            writeTestsTAP_YAMLISH(res, tapFile)
        end
    end
end
