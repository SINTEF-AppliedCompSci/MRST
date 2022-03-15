function settings = struct2settings(oldsettings)

    settings = settingsStruct();
    fnames = fieldnames(oldsettings);

    supported = isfield(settings, fnames);
    if ~all(supported)
       unsupported = sprintf(' * %s\n', fnames{~supported});
       pl = ''; if sum(~supported) ~= 1, pl = 's'; end
       error('Unsupported setting name%s\n%s\nCannot convert to ''settingsStruct''.', pl, unsupported);
    end

    for setting = reshape(fnames, 1, [])
       settings.(setting{1}) = oldsettings.(setting{1});
    end    

end

