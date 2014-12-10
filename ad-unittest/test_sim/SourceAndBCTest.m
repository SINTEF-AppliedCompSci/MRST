classdef SourceAndBCTest < matlab.unittest.TestCase
    properties
        
    end
    
    methods
        function test = SourceAndBCTest(varargin)
            
            
        end
    end
    
    methods
        
        function [G, rock] = simpleGridRock(test)
            G = cartGrid([10, 1, 1], [10 1 1]*meter);
            G = computeGeometry(G);
            
            rock = struct('perm', 1000*ones(G.cells.num, 1)*milli*darcy, ...
                'poro', ones(G.cells.num, 1));
        end
        
        function [state0, model] = getSimpleOW(test)
            [G, rock] = test.simpleGridRock();
            
            mrstModule add ad-props
            fluid = initSimpleADIFluid();
            model = TwoPhaseOilWaterModel(G, rock, fluid);
            state0 = initResSol(G, 0, [.5, .5]);
            state0.wellSol = initWellSolAD([], model, state0);
        end
        
        function [state0, model] = getSimpleBO(test)
            [G, rock] = test.simpleGridRock();
            
            mrstModule add ad-props
            fluid = initSimpleADIFluid();
            model = ThreePhaseBlackOilModel(G, rock, fluid);
            state0 = initResSol(G, 0, [1 1 1]/3);
            state0.rs = 0;
            state0.rv = 0;
            state0.wellSol = initWellSolAD([], model, state0);
        end
        
        function bc = linearPressureDrop(test, G, sat)
            bc = [];
            bc = pside(bc, G, 'xmin', 1*barsa, 'sat', sat);
            bc = pside(bc, G, 'xmax', 0, 'sat', sat);
        end
        
        
        function VerifyDirichletBC(test, state0, model, cases)
            import matlab.unittest.constraints.RelativeTolerance;
            import matlab.unittest.constraints.AbsoluteTolerance;
            import matlab.unittest.constraints.IsEqualTo;
            
            G = model.G;
            for i = 1:size(cases, 1)
                sat = cases(i, :);
                bc = test.linearPressureDrop(model.G, sat);
                
                state = state0;
                solver = NonLinearSolver();
                
                dT = 10*year;
                n = 10;
                for j = 1:n
                    state = solver.solveTimestep(state, dT, model, 'bc', bc);
                end
                % We assume stationary flow, i.e. the bc completely
                % determine the saturations at the endpoint
                expected = repmat(sat, G.cells.num, 1);
                test.verifyThat(state.s, ...
                    matlab.unittest.constraints.IsEqualTo(...
                    expected, 'Within', AbsoluteTolerance(1e-2)) );
                
                % Next, check the linear pressure drop
                expected = interp1(G.faces.centroids(bc.face, 1), bc.value, G.cells.centroids(:, 1));
                
                test.verifyThat(state.pressure, ...
                    matlab.unittest.constraints.IsEqualTo(...
                    expected, 'Within', RelativeTolerance(1e-2)) );
            end
        end
    end
    
    
    methods (Test)
        
        function TestDirichletBCSimpleBO(test)
            [state0, model] = test.getSimpleBO();
            
            cases = [eye(3); [1 1 1]/3];
            test.VerifyDirichletBC(state0, model, cases);
        end
        
        
        function TestDirichletBCSimpleOW(test)
            [state0, model] = test.getSimpleOW();
            
            cases = [0,   1; ...
                0.5, 0.5;...
                1,   0];
            test.VerifyDirichletBC(state0, model, cases);
        end
        
        function states = TestMixedBC(test)
            
        end
        
        function states = TestSourceTerms(test)
            
        end
    end
end

