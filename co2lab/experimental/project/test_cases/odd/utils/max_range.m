function range = max_range(old_range, values)
    if isempty(old_range)
        old_range = [values(1), values(1)];
    end
    
    valmin = min(values(:));
    valmax = max(values(:));
        
    valspan = valmax - valmin;
    
    range = [min(old_range(1), valmin - 0.1 * valspan),...
             max(old_range(2), valmax + 0.1 * valspan)];
end
