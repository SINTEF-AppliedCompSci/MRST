classdef TestSPE3 < ScheduleTest
    properties
        
    end
    
    methods
        function test = TestSPE3(varargin)
            test = test@ScheduleTest('spe3');
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

