function value = getStructField(mstruct, fieldnamelist, defaultValue)

    if ischar(fieldnamelist)
        % handle case where fieldnamelist is just a char
        fieldnamelist = {fieldnamelist};
        if nargin > 2
            value = getStructField(mstruct, fieldnamelist, defaultValue);
        else
            value = getStructField(mstruct, fieldnamelist);
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

    if isempty(mstruct) || (~isfield(mstruct, fieldname) && any(~isprop(mstruct, fieldname)))

        value = UnAssigned();

    else

        if getValue

            value = mstruct.(fieldname);

        else

            value = getStructField(mstruct.(fieldname), fieldnamelist);

        end
        
    end

    if isUnAssigned(value) && nargin > 2
        
        value = defaultValue;
        
    end
        
end
