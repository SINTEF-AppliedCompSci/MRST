function [lumping,subset] = reorganizeIndices(indices)

A = cellfun(@(C)numel(C),indices);

sizeAllIndices = sum(A);

lumping = zeros(sizeAllIndices,1);
subset  = zeros(sizeAllIndices,1);

reel = 1;
for i=1:numel(indices)
    n = numel(indices{i});
    subset(reel:n+reel-1) = indices{i};
    lumping(reel:n+reel-1)= i;
    reel = reel+ n;
end