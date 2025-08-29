function mstruct = removeStructFields(mstruct, varargin)

    for iarg = 1 : numel(varargin)
        
        mstruct = removeStructField(mstruct, varargin{iarg});
        
    end
    
end

