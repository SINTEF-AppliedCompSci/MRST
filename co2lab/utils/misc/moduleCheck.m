function moduleCheck(varargin)
%% Load modules whose names are in the argument list, unless they are loaded already.
    for i = 1:nargin
        try
           require(varargin{i});
        catch
            fprintf('Loading module %s\n', varargin{i});
            mrstModule('add',varargin{i});
        end
    end
end

