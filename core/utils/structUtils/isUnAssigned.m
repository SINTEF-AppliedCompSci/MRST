function status = isUnAssigned(value, fieldnamelist)
    
    if nargin > 1
        status = isAssigned(value, fieldnamelist);
    else
        status = isAssigned(value);
    end

    status = ~status;
    
end
