classdef ScheduleTest < matlab.unittest.TestCase
    properties
        model
        schedule
        rock
        state0
        storeResults
        outputPath
        relativeTolerance
        
        % If the case is very big, some tests will ignore it
        caseIsBig
    end
    
    methods
        function test = ScheduleTest(name, varargin)
            test.storeResults = false;
            test.schedule = [];
            test.state0 = [];
            test.relativeTolerance = 0.01;
            test.caseIsBig = false;
            
            % Attempt to be a clean MRST state
            mrstModule reset ad-unittest ad-blackoil ad-core ad-eor
            gravity reset off

            test.outputPath = fullfile(...
                mrstOutputDirectory(),...
                'TestOutput', ...
                'ScheduleTest'...
                );
            if ~exist(test.outputPath, 'dir')
                mkdir(test.outputPath);
            end
            test = merge_options(test, varargin{:});
            
            if ~(exist(test.outputPath, 'dir') == 7)
                mkdir(test.outputPath);
            end
            % Get model
            [s, m, s0] = getBenchmarkAD(name);
            
            % Avoid outputting a lot of extra data derived from primary
            % variables.
            m.extraStateOutput = false;
            m.extraWellSolOutput = false;
            
            test.state0 = s0;
            test.schedule = s;
            test.model = m;
            test.rock = m.rock;
            
        end
    end
    
    methods
        function [res_states, ref_states] = runSchedule(test, name, varargin)
            [res_states, failure] = runSimulationProblem(test.model, test.state0, ...
                test.schedule, varargin{:});
            test.assertFalse(failure, 'Non-linear solver was unable to complete timestep');
            
            fn = [fullfile(test.outputPath, name), '.mat'];
            if test.storeResults || ~(exist(fn, 'file') == 2)
                disp(['Storing test results for test ''', name, '''']);
                save(fn, 'res_states');
                ref_states = [];
            else
                import matlab.unittest.constraints.RelativeTolerance;
                import matlab.unittest.constraints.AbsoluteTolerance;
                import matlab.unittest.constraints.IsEqualTo;

                % Compare with ref
                tmp = load(fn);
                ref_states = tmp.res_states;
                reltol = RelativeTolerance(test.relativeTolerance);
                abstol = RelativeTolerance(0.1);
                
                for i = 1:numel(ref_states)
                    ref = ref_states{i};
                    res = res_states{i};
                    test.compareStates(ref, res, reltol, abstol, i);
                end
            end
            
            % Store well sols
            if isfield(res_states{1}, 'wellSol')
                d = clock();
                clockstr = [sprintf('%02d', d(4)), sprintf('%02d', d(5))];
                datestr = [date(), '-', clockstr];
                
                fn = [fullfile(test.outputPath, [name, '-wellsols-', datestr]), '.mat'];
                timesteps = test.schedule.step.val;
                wellSols = cellfun(@(x) x.wellSol, res_states, 'UniformOutput', false); %#ok
                save(fn, 'wellSols', 'timesteps');
            end
        end
    
        function compareStates(test, ref, res, reltol, abstol, i)
            import matlab.unittest.constraints.IsEqualTo;
            import matlab.unittest.constraints.RelativeTolerance;
            import matlab.unittest.constraints.AbsoluteTolerance;



            fn_ref = sort(fieldnames(ref));
            fn_res = sort(fieldnames(res));
            test.verifyEqual(fn_ref, fn_res, ...
                ['Field names differ in state ', num2str(i)]);

            for j = 1:numel(fn_ref)
                fn = fn_ref{j};
                e = ['Comperator failed for field ', fn, ' at state number ', num2str(i)];
                
                switch(lower(fn))
                    case 's'
                        s = res.(fn);
                        test.verifyThat(sum(s, 2), ...
                            IsEqualTo(ones(size(s, 1), 1), ...
                            'Within', RelativeTolerance(1e-6)), ...
                            ['Saturations  did not sum to one for state ' num2str(i)]);

                        test.verifyThat(ref.(fn), ...
                            IsEqualTo(res.(fn), 'Within', abstol), ...
                            ['Saturations were not equal for state ' num2str(i)]);
                    case {'flux', 'cmax', 'c'}
                        % Ignore fields intentionally, numerically not that
                        % stable for comparison.
                    case 'wellsol'
                        for k = 1:numel(ref.(fn))
                            % Recursive definition
                            test.compareStates(ref.(fn)(k), res.(fn)(k), reltol, abstol, i);
                        end
                    case {'qws', 'qgs', 'qos', 'qts', 'qs', 'cqs'}
                        % Misc fields that may be close to zero
                        fudge = 1e8;
                        a = round(fudge*res.(fn))/fudge;
                        b = round(fudge*ref.(fn))/fudge;
                        test.verifyThat(a, IsEqualTo(b, 'Within', reltol), e);
                    otherwise
                        test.verifyThat(ref.(fn), ...
                            IsEqualTo(res.(fn), 'Within', reltol), e);
                end
            end
        end
    end

    
    methods (Test)
        
        function states = singleStep(test)
            if test.caseIsBig
                usecpr = true;
            else
                usecpr = false;
            end
            name = test.getIdentifier('singlestep');
            states = test.runSchedule(name, 'stepcount', 1, 'useCPR', usecpr);
        end
        
        function states = baseline(test)
            if test.caseIsBig
                test.assertFail('Test is too big for mldivide direct solver')
            end
            name = test.getIdentifier('baseline');
            states = test.runSchedule(name);
        end
        
        function states = CPR_mldivide(test)
            name = test.getIdentifier('cpr_mldivide');
            states = test.runSchedule(name, 'useCPR', true, 'useAGMG', false);
        end
                
        function states = CPR_AGMG(test)
            name = test.getIdentifier('cpr_agmg');
            mrstModule add agmg
            try
                agmg(speye(3), ones(3, 1));
            catch
                test.verifyFail( ...
                    'AGMG is not installed properly, test cannot proceed')
                return
            end
            states = test.runSchedule(name, 'useCPR', true, 'useAGMG', true);
        end
        function states = selectLinearSolver(test)
            name = test.getIdentifier('select_linear_solver');
            mrstModule add linearsolvers
            states = test.runSchedule(name, 'selectLinearSolver', true);
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
