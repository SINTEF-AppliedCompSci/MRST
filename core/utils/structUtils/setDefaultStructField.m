function mstruct = setDefaultStructField(mstruct, fieldnamelist, defaultValue)

    value = getStructField(mstruct, fieldnamelist);

    if isUnAssigned(value)
        mstruct = setStructField(mstruct, fieldnamelist, defaultValue);
    end
    
end
