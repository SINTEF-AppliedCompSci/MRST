classdef TestEGG < ScheduleTest
    properties
        
    end
    
    methods
        function test = TestEGG(varargin)
            test = test@ScheduleTest('egg');
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

