function status = isUnAssignedOrEmpty(value, fieldnamelist)

    if nargin < 2
        fieldnamelist = {};
    end
    
    value = getStructField(value, fieldnamelist);

    if isUnAssigned(value) || isempty(value)
        status = true;
    else
        status = false;
    end
    
end
