function status = isUnAssigned(value, fieldnamelist)
    
    if nargin < 2
        fieldnamelist = {};
    end
    
    status = isAssigned(value, fieldnamelist);
    status = ~status;
    
end
