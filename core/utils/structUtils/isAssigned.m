function status = isAssigned(value, fieldnamelist)

    if nargin < 2
        fieldnamelist = {};
    end

    value = getStructField(value, fieldnamelist);
    
    if isa(value, 'UnAssigned')
        
        status = false;
        
    else
        
        status = true;
        
    end
    
end
