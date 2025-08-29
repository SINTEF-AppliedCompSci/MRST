function jsonstruct = setDefaultStructField(jsonstruct, fieldnamelist, defaultValue)

    value = getStructField(jsonstruct, fieldnamelist);

    if isUnAssigned(value)
        jsonstruct = setStructField(jsonstruct, fieldnamelist, defaultValue);
    end
    
end
