classdef TestValue < matlab.unittest.TestCase
    methods
        function test = TestValue()
            mrstModule reset
            mrstModule add ad-unittest ad-core
        end
    end
    
    methods (Test)
        function testTrivial(test)
            y = rand(10, 1);
            v = value(y);
            test.assertEqual(y, v)
        end
        
        function testAD(test)
            y = rand(10, 1);
            x = initVariablesADI(y);
            
            v = value(x);
            test.assertEqual(y, v)
        end
        
        function testCellArrayHorizontal(test)
            y = rand(10, 1);
            x = rand(10, 1);
            
            c = {x, y};
            
            v = value(c);
            answer = [x, y];
            test.assertEqual(answer, v)
        end
        
        function testCellArrayVertical(test)
            y = rand(10, 1);
            x = rand(10, 1);
            
            c = {x; y};
            
            v = value(c);
            test.assertEqual(c, v)
        end
        
        function testCellAD(test)
            y = rand(10, 1);
            x = rand(10, 1);
            z = initVariablesADI(y);
            
            c = {x, z};
            
            v = value(c);
            answer = [x, y];
            test.assertEqual(answer, v)
        end
        
        function testStruct(test)
            y = rand(10, 1);
            x = rand(10, 1);
            z = initVariablesADI(y);
            
            s = struct();
            s.c = {x, z};
            s.y = z;
            s.x = x;
            
            v = value(s);
            test.assertEqual(v.c, [x, y])
            test.assertEqual(v.y, y);
            test.assertEqual(v.x, x);
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
