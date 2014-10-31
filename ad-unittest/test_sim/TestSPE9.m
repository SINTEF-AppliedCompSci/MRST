classdef TestSPE9 < ScheduleTest
    properties
        
    end
    
    methods
        function test = TestSPE9(varargin)
            test = test@ScheduleTest('spe9');
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

