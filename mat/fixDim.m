function vec = fixDim(vec)

    if size(vec, 1) == 1
        vec = vec';
    end

end