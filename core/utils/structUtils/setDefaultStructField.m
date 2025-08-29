function jsonstruct = setDefaultJsonStructField(jsonstruct, fieldnamelist, defaultValue)

    value = getJsonStructField(jsonstruct, fieldnamelist);

    if isUnAssigned(value)
        jsonstruct = setJsonStructField(jsonstruct, fieldnamelist, defaultValue);
    end
    
end
