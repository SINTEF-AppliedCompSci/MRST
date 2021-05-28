classdef TestSimulateScheduleAD < matlab.unittest.TestCase
    properties
       storeResults
       outputPath
    end
    
    methods
        function test = TestSimulateScheduleAD(varargin)
            test.storeResults = false;
            
            % Attempt to be a clean MRST state
            mrstModule reset ad-unittest ad-core ad-blackoil
            gravity reset off
            mrstVerbose off

            test.outputPath = fullfile(...
                mrstPath('query', 'ad-unittest'),...
                'output', ...
                'simulateScheduleAD'...
                );
            
            test = merge_options(test, varargin{:});
            
            if ~(exist(test.outputPath, 'dir') == 7)
                mkdir(test.outputPath);
            end
        end
    end
    
    methods (Static)

    end
    
    methods (Test)
        function testMinistepOutput(test)
            [state0, schedule, model] = setupSimpleOW();
            
            % Disable flux output for comparison purposes
            model.outputFluxes = false;
            model.extraStateOutput = false;
            model.extraWellSolOutput = false;
            % Override schedule with a single timestep
            schedule.step.val = 1*day;
            schedule.step.control = 1;
            
            substepno = 5;
            ts = SimpleTimeStepSelector('maxTimestep', schedule.step.val/substepno, ...
                                        'verbose', true);
            
            nonlinear = NonLinearSolver('TimeStepSelector', ts,...
                                        'verbose', true);
            
                                    
            solve = @(ministepOn)  ...
                simulateScheduleAD(state0, model, schedule,...
                            'OutputMinisteps', ministepOn, ...
                            'NonLinearSolver', nonlinear);

            [ws_mini, s_mini, r_mini] = solve(true);
            test.assertFalse(r_mini.Failure) 
            
            [ws, s, r] = solve(false);
            test.assertFalse(r.Failure) 
            
            [ws_1s, s_1s, r_1s] = simulateScheduleAD(state0, model, schedule);
            
            test.verifyLength(ws_mini, substepno, ...
                'Number of ministep wellSol outputs did not match target ministep count.')
            
            test.verifyLength(s_mini, substepno, ...
                'Number of ministep state outputs did not match target ministep count.')

            test.verifyLength(ws, numel(schedule.step.val),...
                'Number of wellSol outputs did not match schedule.')
            
            test.verifyLength(s, numel(schedule.step.val),...
                'Number of wellSol outputs did not match schedule.')
            
            import matlab.unittest.constraints.RelativeTolerance;
            import matlab.unittest.constraints.IsEqualTo;


            % First, check that the final output of ministeps is equal to
            % the exactly same simulation with only output at control steps
            reltol_strict = RelativeTolerance(1e-6);
            
            test.verifyThat(s{end}, IsEqualTo(s_mini{end}, 'Within', reltol_strict), ...
                'Final reservoir state using ministeps differed from using a single report!')
            
            test.verifyThat(ws{end}, IsEqualTo(ws_mini{end}, 'Within', reltol_strict), ...
                'Final well state using ministeps differed from using a single report!')
            
            % Verify that doing several smaller steps is approximately
            % equal to a single long step
            reltol_loose = RelativeTolerance(0.1);
            
            test.verifyThat(ws{end}, IsEqualTo(ws_1s{end}, 'Within', reltol_loose), ...
                'Final well state using ministeps differed from a single step!')
            
            test.verifyThat(s{end}, IsEqualTo(s_1s{end}, 'Within', reltol_loose), ...
                'Final reservoir state using ministeps differed from a single step!')
        end
        
        
        function testTimestepSelection(test)
            
            
            opt = {'maxRelativeAdjustment', inf,...
                   'minRelativeAdjustment', sqrt(eps), ...
                   'verbose', true};
            
            ts = SimpleTimeStepSelector('maxTimestep', 10*day, ...
                                        opt{:});
            
            sel = {ts};
            
            sel = {sel{:}, IterationCountTimeStepSelector('targetIterationCount', 3, ...
                                                      opt{:})};
            
            mrstVerbose on
            
            ok = test.runTimestepSelector([]);
            test.assertFalse(ok,...
            ['Base case with long timestep converged without timestep', ...
            ' control! Resulting test data is useless.']);
            
            for i = 1:numel(sel)
                ts = sel{i};
                ok = test.runTimestepSelector(ts);
                test.verifyTrue(ok, ...
                    ['Time step selector ', class(ts), ...
                    ' failed to compute valid timesteps']);
            end
            
        end
        
    end
    methods
        function ok = runTimestepSelector(test, ts)
            [state0, schedule, model] = setupSimpleOW();
            
            % Set tolerances to a strict value
            model.nonlinearTolerance = 1e-8;
            
            % Override schedule with a single timestep
            schedule.step.val = [10*day; 10*day; 1*year];
            schedule.step.control = [1; 1; 1];

            
            nonlinear = NonLinearSolver('TimeStepSelector', ts,...
                                        'maxTimestepCuts', 0, ...
                                        'errorOnFailure', false, ...
                                        'maxIterations', 8, ...
                                        'verbose', true);
            
            
            [ws, s, r] = simulateScheduleAD(state0, model, schedule, 'NonLinearSolver', nonlinear);
            ok = ~r.Failure;
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
