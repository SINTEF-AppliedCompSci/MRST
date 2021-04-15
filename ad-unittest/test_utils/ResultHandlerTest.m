classdef ResultHandlerTest < matlab.unittest.TestCase
    properties
        
    end
    
    methods
        function test = ResultHandlerTest()
            mrstModule reset
            mrstModule add ad-unittest ad-core
        end
    end
    
    methods
        function d = makeTestData(test, id)
            d = struct('identifier', id);
        end
        
        function runRetrievalTests(test, handler)
            A = test.makeTestData(1);
            B = test.makeTestData('hello');
            
            
            index = {1, [3, 5, 6], [3, 5, 6], [3, 5, 6], []};
            input = {A, {A, B, A}, A        , {A, B},    []};
            valid = [1,  1       , 1        , 0,         1];
            
            for i = 1:numel(index)
                handler.resetData();
                
                ind = index{i};
                in = input{i};
                try
                    handler{ind} = in;
                catch e
                    if ~valid(i)
                        test.verifyError(@() rethrow(e), ...
                            'ResultHandler:dimMismatch', ...
                            'Dimension assignment mismatch');
                        continue
                    end
                end
                switch(numel(ind))
                    case 1 
                        out = handler{ind};
                    case 0
                        continue
                    otherwise
                        out = handler(ind);
                end
                test.verifyEqual(out, in);
            end
        end
            
        function testInput(test, index, input, memory)
            handler = ResultHandler('storeInMemory', memory, ...
                                    'writeToDisk',  ~memory);
            if ~memory
                handler.resetData();
            end
            tmp = {};
            
            if iscell(input)
                % Interpret as multiple value to be stored
                handler(index) = input;
                tmp(index) = input;
            else
                % Interpret as single value to be stored
                handler{index} = input;
                tmp{index} = input;
            end

            

            % Test curly braces
            out = {handler{index}};
            test.verifyEqual(out, {tmp{index}});
            if numel(index) > 0
                % Test paranthesis
                out = handler(index);
                test.verifyEqual(out, tmp(index));
            end
        end
        
        function testGetValidIds(test, memory)
            handler = ResultHandler('storeInMemory', memory, ...
                                    'writeToDisk',  ~memory);
            
            if ~memory
                % May exist remnants of old runs
                handler.resetData();
            end
            A = test.makeTestData(1);
            
            td = [1 32 103];
            for i = td
                handler{i} = A;
            end
            
            ids = handler.getValidIds();

            test.assertEqual(ids, td)
            
            handler.resetData();
            ids = handler.getValidIds();
            test.assertTrue(isempty(ids));
        end
    end
    
    
    methods (Test)
        
        function SingleInputMemory(test)
            d = test.makeTestData(1);
            test.testInput(1, d, true)
        end
        
        function MultipleInputMemory(test)
            A = test.makeTestData(1);
            B = test.makeTestData('hello!');
            
            test.testInput([1,2], {A, B}, true)
        end
        
        function BadInputMemory(test)
            % Verify that assignment mismatch goes wrong
            A = test.makeTestData(1);
            B = test.makeTestData('hello!');
            try
                test.testInput([1, 2], {A, B, A}, true)
            catch e
                test.verifyError(@() rethrow(e), ...
                                'ResultHandler:dimMismatch', ...
                                'Dimension assignment mismatch');
            end
        end
        
        function SingleInputDisk(test)
            d = test.makeTestData(1);
            test.testInput(1, d, false)
        end
        
        function MultipleInputDisk(test)
            A = test.makeTestData(1);
            B = test.makeTestData('hello!');
            
            test.testInput([1,2], {A, B}, false)
        end
        
        function BadInputDisk(test)
            % Verify that assignment mismatch goes wrong
            A = test.makeTestData(1);
            B = test.makeTestData('hello!');
            try
                test.testInput([1, 2], {A, B, A}, false)
            catch e
                test.verifyError(@() rethrow(e), ...
                                'ResultHandler:dimMismatch', ...
                                'Dimension assignment mismatch');
            end
        end
        
        function GetValidIdsMemory(test)
            test.testGetValidIds(true);
        end
        
        function GetValidIdsDisk(test)
            test.testGetValidIds(false);
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
