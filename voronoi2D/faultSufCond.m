function [p, removed] = faultSufCond(p, F)
  TOL = 50*eps;
  nc = size(F.c.CC,1);
  np = size(p,1);
  
  CRSqr = F.c.R.^2;
  removed = zeros(np,1);
  for i = 1:nc
    distSqr = sum(bsxfun(@minus, F.c.CC(i,:), p).^2,2);
    removed = removed + (distSqr<CRSqr(i)-TOL);
  end
  removed = logical(removed);
  p = p(~removed,:);
end