classdef ResultHandlerTest < matlab.unittest.TestCase
    properties
        
    end
    
    methods
        function test = ResultHandlerTest(varargin)
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
    end
    
    methods (Test)
        function StoreInMemoryTest(test)
            handler = ResultHandler('storeInMemory', true);
            test.runRetrievalTests(handler);
        end

    end

end