function p = convertTo3DPlane(p, T, R)
% Convert the points p from horizontal plane to the fully 3D plane. T and R
% are transformation matrix, can be obtained by `convertToXYPlane`
    p_extend = [p, ones(size(p,1),1)];
    p = p_extend / (T*R);
    p = p(:, 1:3);
end