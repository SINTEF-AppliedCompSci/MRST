function status = isUnAssignedOrEmpty(value, fieldnamelist)

    if isempty(getStructField(value, fieldnamelist))
        status = true;
        return
    end
    
    status = isAssigned(value, fieldnamelist);
    status = ~status;
    
end
