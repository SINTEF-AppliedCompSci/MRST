function jsonstruct = removeJsonStructFields(jsonstruct, varargin)

    for iarg = 1 : numel(varargin)
        
        jsonstruct = removeJsonStructField(jsonstruct, varargin{iarg});
        
    end
    
end

