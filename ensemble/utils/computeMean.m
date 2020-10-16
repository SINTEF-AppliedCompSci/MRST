function mean = computeMean(meanA, meanB, nA, nB)
    % Compute from two means of different sizes
    mean = (nA.*meanA + nB.*meanB)./(nA + nB);
end