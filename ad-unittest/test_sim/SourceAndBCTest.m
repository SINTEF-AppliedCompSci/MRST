classdef SourceAndBCTest < matlab.unittest.TestCase
    properties
        
    end
    
    methods
        function test = SourceAndBCTest(varargin)
            mrstModule reset
            mrstModule add ad-core ad-blackoil ad-fi ad-unittest
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
            model = model.validateModel();
            state0 = validateState(model, state0);
        end
        
        function [state0, model] = getSimpleBO(test)
            [G, rock] = test.simpleGridRock();
            
            mrstModule add ad-props
            fluid = initSimpleADIFluid();
            model = ThreePhaseBlackOilModel(G, rock, fluid);
            state0 = initResSol(G, 0, [1 1 1]/3);
            model = model.validateModel();
            state0 = validateState(model, state0);
        end
        
        function bc = linearPressureDrop(test, G, sat)
            bc = [];
            bc = pside(bc, G, 'xmin', 1*barsa, 'sat', sat);
            bc = pside(bc, G, 'xmax', 0, 'sat', sat);
        end
        
        function bc = mixedBC(test, G, sat, time, sgn)
            injRate = 5*sgn*sum(G.cells.volumes)/time;
            bc = [];
            bc = fluxside(bc, G, 'xmin', -injRate, 'sat', sat);
            bc = pside(bc, G, 'xmax', 0, 'sat', sat);
        end
        
        function [bc, src] = srcWithBC(test, G, sat, time)
            [src, bc] = deal([]);
            injRate = 5*sum(G.cells.volumes)/time;
            src = addSource(src, ceil(G.cells.num/2), injRate, 'sat', sat);
            src = addSource(src, G.cells.num, -injRate, 'sat', sat);
            % BC to fixate pressure - avoid singular values
            bc = pside(bc, G, 'xmin', 0*barsa, 'sat', sat);

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
        
        function VerifyMixedBC(test, state0, model, cases, sgn)
            import matlab.unittest.constraints.RelativeTolerance;
            import matlab.unittest.constraints.AbsoluteTolerance;
            import matlab.unittest.constraints.IsEqualTo;
            
            G = model.G;
            for i = 1:size(cases, 1)
                sat = cases(i, :);
                                
                dT = 10*year;
                n = 10;
                
                bc = test.mixedBC(model.G, sat, dT*n, sgn);
                
                state = state0;
                solver = NonLinearSolver();

                for j = 1:n
                    state = solver.solveTimestep(state, dT, model, 'bc', bc);
                end
                % We assume stationary flow, i.e. the bc completely
                % determine the saturations at the endpoint
                expected = repmat(sat, G.cells.num, 1);
                test.verifyThat(state.s, ...
                    matlab.unittest.constraints.IsEqualTo(...
                    expected, 'Within', AbsoluteTolerance(1e-2)) );
            end
        end
        
        function VerifySRC(test, state0, model, cases)
            import matlab.unittest.constraints.RelativeTolerance;
            import matlab.unittest.constraints.AbsoluteTolerance;
            import matlab.unittest.constraints.IsEqualTo;
            
            G = model.G;
            for i = 1:size(cases, 1)
                sat = cases(i, :);
                                
                dT = 10*year;
                n = 10;
                
                [bc, src] = test.srcWithBC(G, sat, dT*n);
                
                state = state0;
                solver = NonLinearSolver();

                for j = 1:n
                    state = solver.solveTimestep(state, dT, model, 'bc', bc, 'src', src);
                end
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
        
        function TestMixedBC_BO(test)
            [state0, model] = test.getSimpleBO();
            
            cases = [eye(3); [1 1 1]/3];
            
            test.VerifyMixedBC(state0, model, cases, -1);
            test.VerifyMixedBC(state0, model, cases, +1);
        end

        function TestMixedBC_OW(test)
            [state0, model] = test.getSimpleOW();
            
            cases = [0,   1; ...
                    0.5, 0.5;...
                    1,   0];
            test.VerifyMixedBC(state0, model, cases, -1);
            test.VerifyMixedBC(state0, model, cases, +1);
        end
        
        function TestSourceTermsBO(test)
            [state0, model] = test.getSimpleBO();
            cases = [eye(3); [1 1 1]/3];
            test.VerifySRC(state0, model, cases);
        end
        
        function TestSourceTermsOW(test)
            [state0, model] = test.getSimpleOW();
            cases = [eye(2); [1 1]/2];
            test.VerifySRC(state0, model, cases);
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
