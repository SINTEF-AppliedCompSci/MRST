function d=minPdist2(X, Y)
% Calculate the distance from one point set to the nearest point in another
%
% SYNOPSIS
%     d = minPdist2(x, x)
% Arguments:
%    x        n x d array holding n points
%    y        n x d array holding m points
% Returns:
%    d        d(i) is the distance from x(i,:) to the nearest point y
%
% Written by Knut-Andreas Lie
d = sqrt(min(abs(bsxfun(@plus, sum(X.^2,2),sum(Y.^2,2)') - 2*(X*Y')),[],2));
end