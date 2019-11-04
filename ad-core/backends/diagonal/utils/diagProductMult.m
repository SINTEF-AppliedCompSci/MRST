function [x, D1, D2] = diagProductMult(v1, v2, x, y, D1, D2)
    [x, D2] = diagMult(v2, x, D2);
    [y, D1] = diagMult(v1, y, D1);
    x = x + y;
end