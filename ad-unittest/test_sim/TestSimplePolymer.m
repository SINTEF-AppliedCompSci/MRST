classdef TestSimplePolymer < ScheduleTest
    properties
        
    end
    
    methods
        function test = TestSimplePolymer(varargin)
            test = test@ScheduleTest('simplepolymer');
            test = merge_options(test, varargin{:});
        end
        function s = getIdentifier(test, name) %#ok
            s = [mfilename('class'), '_', name];
        end
    end
    
    methods (Test)
        % Add your own, test specific tests here
    end
    
end

