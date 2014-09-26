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
        function test = ScheduleTest(varargin)
            test.storeResults = false;
            test.schedule = [];
            test.state0 = [];
            test.relativeTolerance = 0.01;
            test.caseIsBig = false;
            
            % Attempt to be a clean MRST state
            mrstModule reset ad-unittest ad-blackoil ad-core
            gravity reset off

            test.outputPath = fullfile(...
                mrstPath('query', 'ad-unittest'),...
                'output', ...
                'schedule'...
                );
            
            test = merge_options(test, varargin{:});
            
            if ~(exist(test.outputPath, 'dir') == 7)
                mkdir(test.outputPath);
            end            
        end
    end
    
    methods
        function [states, ref] = runSchedule(test, name, varargin)
            [states, failure] = runSimulationProblem(test.model, test.state0, ...
                test.schedule, varargin{:});
            test.assertFalse(failure, 'Non-linear solver was unable to complete timestep');
            
            fn = [fullfile(test.outputPath, name), '.mat'];
            if test.storeResults || ~(exist(fn, 'file') == 2)
                disp(['Storing test results for test ''', name, '''']);
                save(fn, 'states');
                ref = [];
            else
                import matlab.unittest.constraints.RelativeTolerance;
                import matlab.unittest.constraints.IsEqualTo;

                % Compare with ref
                tmp = load(fn);
                ref = tmp.states;
                reltol = RelativeTolerance(test.relativeTolerance);
                
                if isfield(ref{1}, 's')
                    for i = 1:numel(ref)
                        % Validate sum of saturations for solution -
                        % physical constraint.
                        s = states{i}.s;
                        test.verifyThat(sum(s, 2), ...
                            IsEqualTo(ones(size(s, 1), 1), ...
                            'Within', RelativeTolerance(1e-6)), ...
                            'Saturations  did not sum to one.');
                        
                        % Avoid zero comparison for saturations by adding a
                        % reasonable constant before comparison
                        states{i}.s = states{i}.s + .1;
                        ref{i}.s = ref{i}.s + .1;
                    end
                end
                test.verifyThat(states, IsEqualTo(ref,...
                                    'Within', reltol));
            end
        end
        
    end
    
    methods (Test)
        
        function singleStep(test)
            if test.caseIsBig
                usecpr = true;
            else
                usecpr = false;
            end
            name = test.getIdentifier('singlestep');
            test.runSchedule(name, 'stepcount', 1, 'useCPR', usecpr);
        end
        
        function baseline(test)
            if test.caseIsBig
                test.assertFail('Test is too big for mldivide direct solver')
            end
            name = test.getIdentifier('baseline');
            test.runSchedule(name);
        end
        
        function CPR_mldivide(test)
            name = test.getIdentifier('cpr_mldivide');
            test.runSchedule(name, 'useCPR', true, 'useAGMG', false);
        end
                
        function CPR_AGMG(test)

            name = test.getIdentifier('cpr_agmg');
            mrstModule add agmg
            try
                agmg(speye(3), ones(3, 1));
            catch
                test.verifyFail( ...
                    'AGMG is not installed properly, test cannot proceed')
                return
            end
            test.runSchedule(name, 'useCPR', true, 'useAGMG', true);
        end
    end

end