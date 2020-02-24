function ip = innerProduct(u, v, cubature, elements, type)
    if isa(u, 'SpatialVector') || isa(v, 'SpatialVector')
        uv = dot(u, v);
    else
        uv = u.*v;
    end
    switch type
        case 'cell'
            ip = cellIntegral(uv, cubature, elements);
        case 'face'
            ip = faceIntegral(uv, cubature, elements);
    end
end