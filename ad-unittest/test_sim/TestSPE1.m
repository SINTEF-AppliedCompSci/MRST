classdef TestSPE1 < ScheduleTest
    properties
        
    end
    
    methods
        function test = TestSPE1(varargin)
            test = test@ScheduleTest('spe1');
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

