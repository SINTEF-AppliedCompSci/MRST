function [horizons, ix] = sortHorizons(horizons)
    z = cellfun(@(x) min(x.z(:)), horizons);
    [z, ix] = sort(z);
    
    horizons = horizons(ix);
end