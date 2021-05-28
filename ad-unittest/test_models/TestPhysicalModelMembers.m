classdef TestPhysicalModelMembers < matlab.unittest.TestCase
    methods
        function test = TestPhysicalModelMembers()
            mrstModule reset
            mrstModule add ad-unittest ad-core
        end
    end
    
    methods (Test)
        % Add your own, test specific tests here
        function badVariableInput(test)
            % Try to retrieve an variable not known by the model, check for
            % failure.
            import matlab.unittest.qualifications.Assertable;
            model = MockPhysicalModel();
            state = struct('unknownfield', 1);
            fn = @() model.getProp(state, 'unknownfield');
            test.assertError(fn, 'PhysicalModel:UnknownVariable');
        end
        
        function firstColumnRetrieval(test)
            model = MockPhysicalModel();
            d = rand(10, 2);
            state = struct('x', d);
            % x1 is the first column according to the mock model definition
            test.assertEqual(d(:,1), model.getProp(state, 'x1'));
        end
        
        function secondColumnRetrieval(test)
            model = MockPhysicalModel();
            d = rand(10, 2);
            state = struct('x', d);
            % x1 is the second column according to the mock model definition
            test.assertEqual(d(:,2), model.getProp(state, 'x2'));
        end
        function allColumnsRetrieval(test)
            model = MockPhysicalModel();
            d = rand(10, 2);
            state = struct('x', d);
            % All columns
            test.assertEqual(d, model.getProp(state, 'x'));
        end
        
        function updateState(test)
            n = 10;
            state = struct('x', rand(n, 2), 'y', rand(n,1));
            
            model = MockPhysicalModel();
            
            % Make mock newton increments for second column of x and y
            problem = LinearizedProblem({'x2', 'y'}, state, 0);
            dx = {ones(n, 1), 2*ones(n, 1)};
            statenew = model.updateState(state, problem, dx, []);
            
            % Check that resulting update is equal to the manually
            % calculated values
            test.assertEqual(statenew.x(:,2), model.getProp(state, 'x2') + dx{1})
            test.assertEqual(statenew.y, model.getProp(state, 'y') + dx{2})
        end
        
        function limitUpdateRelative(test)
            val = [100; 10; 1];
            dv  = [1; 1; 1];
            model = MockPhysicalModel();
            % 10% maximum relative update
            maxRelativeUpdate = .1;
            [dvnew, change] = model.limitUpdateRelative(dv, val, maxRelativeUpdate);
            
            test.assertTrue(all(abs((dvnew + val)./abs(val)) <= 1.1));
            test.assertTrue(all(abs((dvnew + val)./abs(val)) >= 0.9));
            
            relUp = (dvnew + val)./val;
            % Check that relative change is either the relative update or
            % the unlimited update.
            ok = relUp == (1 + maxRelativeUpdate);
            ok = ok | relUp == (dv + val)./val;
            test.assertTrue(all(ok))
        end
        
        function limitUpdateAbsolute(test)
            for sgn = [-1, 1]
                model = MockPhysicalModel();
                dv = sgn.*[0; 1; 5; 10];
                maxAbsUp = 5;
                dvnew = model.limitUpdateAbsolute(dv, maxAbsUp);
                test.assertTrue(all(abs(dvnew) <= maxAbsUp));
            end
        end
        
        function getPropsSingle(test)
            model = MockPhysicalModel();
            d = rand(10, 2);
            state = struct('x', d);
            % x1 is the first column according to the mock model definition
            test.assertEqual(d(:,1), model.getProps(state, 'x1'));
        end
        
        function getPropsMultiple(test)
            model = MockPhysicalModel();
            d = rand(10, 2);
            state = struct('x', d, 'y', 2*d(:,1));
            % x1 is the first column according to the mock model definition
            [x, y] = model.getProps(state, 'x', 'y');
            
            test.assertEqual(state.x, x);
            test.assertEqual(state.y, y);
        end
        
        function testInsertionAD(test)
            model = MockPhysicalModel();
            d = rand(10, 2);
            state = struct('x', d, 'y', 2*d(:,1));
            
            tmp = initVariablesADI(d(:, 1));
            state = model.setProp(state, 'x2', tmp);
            
            x2 = model.getProp(state, 'x2');
            test.assertEqual(x2, tmp)
        end
        
        function testSingleAD(test)
            model = MockPhysicalModel();
            d = rand(10, 2);
            state = struct('x', d, 'y', 2*d(:,1));
            
            tmp = initVariablesADI(d(:, 1));
            state = model.setProp(state, 'y', tmp);
            
            y = model.getProp(state, 'y');
            test.assertEqual(y, tmp)
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
