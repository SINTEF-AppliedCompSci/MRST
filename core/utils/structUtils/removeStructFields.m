function jsonstruct = removeStructFields(jsonstruct, varargin)

    for iarg = 1 : numel(varargin)
        
        jsonstruct = removeStructField(jsonstruct, varargin{iarg});
        
    end
    
end

