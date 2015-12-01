function Seff = getStorageEfficiency(fmName)
% returns storage efficiency (%) of formation. Values come from NPD's CO2
% Storage Atlas

    if strcmpi(fmName,'G_gsf') || strcmpi(fmName,'G_albatross') ...
            || strcmpi(fmName,'G_askeladd')
        Seff = 3;
        
    elseif strcmpi(fmName,'G_prospectC') || strcmpi(fmName,'G_prospectD') ...
            || strcmpi(fmName,'G_prospectE') || strcmpi(fmName,'G_prospectF')
        Seff = 10;
        
    elseif strcmpi(fmName,'G_prospectG') || strcmpi(fmName,'G_prospectH')
        Seff = 5;
        
    else
        warning('\n No storage efficiency value available for %s. \n', fmName)
        Seff = [];
    end


end

