function value = getStructField(jsonstruct, fieldnamelist, defaultValue)

    if ischar(fieldnamelist)
        % handle case where fieldnamelist is just a char
        fieldnamelist = {fieldnamelist};
        if nargin > 2
            value = getStructField(jsonstruct, fieldnamelist, defaultValue);
        else
            value = getStructField(jsonstruct, fieldnamelist);
        end
        
        return
    end

    fieldname = fieldnamelist{1};

    if numel(fieldnamelist) > 1

        fieldnamelist = fieldnamelist(2 : end);
        getValue = false;
        
    else
        
        getValue = true;
        
    end

    if isempty(jsonstruct) || (~isfield(jsonstruct, fieldname) && ~isprop(jsonstruct, fieldname))

        value = UnAssigned();

    else

        if getValue

            value = jsonstruct.(fieldname);

        else

            value = getStructField(jsonstruct.(fieldname), fieldnamelist);

        end
        
    end

    if isUnAssigned(value) && nargin > 2
        
        value = defaultValue;
        
    end
        
end
