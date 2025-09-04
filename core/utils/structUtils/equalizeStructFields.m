function [mstruct, value] = equalizeStructFields(mstruct, fieldnamelists, defaultvalue)

    if nargin < 3
        defaultvalue = UnAssigned();
    end
    
    values = cell(numel(fieldnamelists), 1);

    for ifd = 1 : numel(fieldnamelists)
        fieldnamelist = fieldnamelists{ifd};
        values{ifd} = getStructField(mstruct, fieldnamelist);
    end

    assigned_inds = cellfun(@(v) isAssigned(v), values);
    assigned_inds = find(assigned_inds);
    
    if isempty(assigned_inds)

        value = defaultvalue;

    else

        fieldnamelist = fieldnamelists{assigned_inds(1)};
        value = getStructField(mstruct, fieldnamelist);

        for iass = 2 : numel(assigned_inds)
            
            fieldnamelist = fieldnamelists{assigned_inds(iass)};
            v = getStructField(mstruct, fieldnamelist);
            if ~isequal(v, value)
                error('The values of the specified fields are not equal.');
            end
            
        end

    end

    for ifd = 1 : numel(fieldnamelists)
        fieldnamelist = fieldnamelists{ifd};
        mstruct = setStructField(mstruct, fieldnamelist, value);
    end

end
