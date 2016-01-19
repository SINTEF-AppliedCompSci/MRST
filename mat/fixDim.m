function vec = fixDim(vec)

    if size(vec, 2) == 1
        vec = vec';
    end

end