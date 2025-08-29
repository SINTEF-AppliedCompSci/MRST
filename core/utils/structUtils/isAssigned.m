function status = isAssigned(value, fieldnamelist)

    if nargin > 1

        % the value input is then a jsonstruct, we extract the value and return the isAssign flag
        value = getJsonStructField(value, fieldnamelist);

        status = isAssigned(value);
        
        return
    end
    
    if isa(value, 'UnAssigned')
        
        status = false;
        
    else
        
        status = true;
        
    end
    
end
