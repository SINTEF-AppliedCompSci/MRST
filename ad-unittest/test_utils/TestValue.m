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

