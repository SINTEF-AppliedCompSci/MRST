function var = computeVariance(varA, varB, meanA, meanB, nA, nB, normfn)
    if nargin < 7
        normfn = @(u) abs(u);
    end
    var = (nA-1)*varA + (nB-1)*varB + nA*nB*normfn(meanB - meanA)^2/(nA + nB);
    var = var/max(nA + nB - 1,1);
end