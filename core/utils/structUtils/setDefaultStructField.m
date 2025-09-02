function mstruct = setDefaultStructField(mstruct, fieldnamelist, defaultValue)

    if isUnAssigned(mstruct)
        mstruct = setStructField([], fieldnamelist, defaultValue);
        return
    end
    
    value = getStructField(mstruct, fieldnamelist);

    if isUnAssigned(value)
        mstruct = setStructField(mstruct, fieldnamelist, defaultValue);
    end
    
end
