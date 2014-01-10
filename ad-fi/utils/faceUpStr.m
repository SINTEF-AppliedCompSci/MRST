function upsExpr = faceUpStr(flag, expr)
N = getGridNeighbors;
inx = N(:,2);
inx(flag) = N(flag,1);
upsExpr = expr(inx);
end
