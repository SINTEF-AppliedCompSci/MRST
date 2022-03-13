classdef RegressionTestGroup
    properties
        name
        tests
        verbose
    end
    
    methods
        %-----------------------------------------------------------------%
        function group = RegressionTestGroup(name, tests, varargin)
            [group, extra] = merge_options(group, varargin{:});
            opt = struct('regTestOpt', {{}});
            opt = merge_options(opt, extra{:});
            group.name = name;
            group.tests = cell(numel(tests), 1);
            
            if isempty(opt.regTestOpt) || ~iscell(opt.regTestOpt{1})
                opt.regTestOpt = repmat({opt.regTestOpt}, group.numTests, 1);
            end
            for i = 1:group.numTests()
                test = tests{i};
                if ischar(test) || isa(test, 'TestCase')
                    test = RegressionTest(test,                 ...
                                          'group', group.name , ...
                                          opt.regTestOpt{i}{:});
                elseif isa(test, 'RegressionTest')
                    test.group = group.name;
                else
                    error(['Each test must be either strings corresponding ', ...
                           'to test case setup functions, or instances ', ...
                           'of the TestCase class'])
                end
                group.tests{i} = test;
            end
        end
        
        %-----------------------------------------------------------------%
        function report = runRegressionTests(group)
            report = struct('passed', 1);
            for i = 1:group.numTests()
                testReport = group.tests{i}.runRegressionTest();
                report.(group.tests{i}.name) = testReport;
                report.passed = min(report.passed, testReport.passed);
            end
            group.printRegressionTestReport(report);
        end

        %-----------------------------------------------------------------%
        function n = numTests(group)
            n = numel(group.tests);
        end
        
        %-----------------------------------------------------------------%
        function printRegressionTestReport(group, report)
            switch report.passed
                case -1
                    str = 'inconclusive';
                case 0
                    str = 'failed';
                case 1
                    str = 'passed';
            end
            strHead = sprintf('Regression test report (overall status: %s)', str);
            names = setdiff(fieldnames(report)', {'passed', 'name'});
            nn = max(cellfun(@numel, [names, 'Test case']));
            npf = max(cellfun(@numel, {'passed', 'failed', 'inconclusive'}));

            headers = {'Test case', 'Status', 'Details'};

            v  = struct2cell(report);
            nd = max(cellfun(@(tr) numel(tr.details), v(cellfun(@isstruct, v))));
            nc = [nn, npf, nd];
            num = max(sum(nc+3) + 1, numel(strHead) + 4);            
            nc(end) = nc(end) + num - sum(nc+2) - 4;
            
            fprintf('\n\n');
            fprintf([repmat('=', 1, num), '\n']);
            fprintf(['| %-', num2str(num-3), 's|\n'], strHead);
            fprintf([repmat('=', 1, num), '\n']);
            fprintf('|');
            for i = 1:numel(headers)
                fprintf([' %-', num2str(nc(i)), 's |'], headers{i});
            end
            fprintf(['\n', repmat('=', 1, num), '\n']);
            for i = 1:numel(names)
                fprintf('|');
                fprintf([' %-', num2str(nc(1)), 's |'], names{i});
                testReport = report.(names{i});
                pf = group.flag2str(testReport.passed);
                fprintf([' %-', num2str(nc(2)), 's |'], pf);
                fprintf([' %-', num2str(nc(3)), 's |'], testReport.details);
                fprintf('\n');
            end
            fprintf([repmat('=', 1, num), '\n\n']);
        end
        
    end
    
    methods (Static)        
        %-----------------------------------------------------------------%
        function str = flag2str(flag)
            switch flag
                case -1
                    str = 'inconclusive';
                case 0
                    str = 'failed';
                case 1
                    str = 'passed';
            end
        end
    end

end

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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
