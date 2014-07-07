classdef ScheduleTest < matlab.unittest.TestCase
    properties
        model
        schedule
        state0
        storeResults
        outputPath
        relativeTolerance
    end
    
    methods
        function test = ScheduleTest(varargin)
            test.storeResults = false;
            test.schedule = [];
            test.state0 = [];
            test.relativeTolerance = 0.01;
            
            % Attempt to be a clean MRST state
            mrstModule reset ad-unittest
            gravity reset off
            mrstVerbose off

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
                ref = load(fn);
                reltol = RelativeTolerance(test.relativeTolerance);
                
                test.verifyThat(states, IsEqualTo(ref.states,...
                                    'Within', reltol));
            end
        end
        
        function [G, rock, fluid, deck, schedule] = setupADcase(test, fn)
            moddir = mrstPath('query', 'ad-testdata');
            fn = fullfile(moddir, fn);
            if ~exist(fn, 'file')
                error(['Did not find dataset at expected location: (', fn , ')'])
            end
            
            deck = readEclipseDeck(fn);
            deck = convertDeckUnits(deck);
            
            G = initEclipseGrid(deck);
            G = computeGeometry(G);
            
            rock  = initEclipseRock(deck);
            rock  = compressRock(rock, G.cells.indexMap);
            
            fluid = initDeckADIFluid(deck);
            schedule = convertDeckScheduleToMRST(G, rock, deck);
        end
    end
end