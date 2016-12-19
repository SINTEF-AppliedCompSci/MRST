function [x, isAD] = getSampleAD(varargin)
% Utility for getting a AD value if it exists from a list of possible
% AD-values
    x = varargin{1};
    isAD = false;
    for i = 1:numel(varargin)
        if isa(varargin{i}, 'ADI')
            x = varargin{i};
            isAD = true;
            return
        end
    end
end