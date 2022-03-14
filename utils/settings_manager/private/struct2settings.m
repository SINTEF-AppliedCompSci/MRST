function settings = struct2settings(oldsettings)

    settings = settingsStruct();
    fnames = fieldnames(oldsettings);
    
    for i = 1:numel(fnames)
        assert(isfield(settings,fnames{i}),...
            sprintf('Unknown field %s in mrstsettings.mat, cannot convert to settingsStruct',fnames{i}));
        settings.(fnames{i}) = oldsettings.(fnames{i});
    end
    
end

