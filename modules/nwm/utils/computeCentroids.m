function pmid = computeCentroids(p)
% Compute centroids of the 2D polygon specified by points p
    p0 = sum(p,1)/size(p,1);
    areas = zeros(size(p,1),1);
    pmids = zeros(size(p,1),2);
    p = [p; p(1,:)];
    for i = 1 : length(areas)
        pts = [p(i:i+1,:); p0];
        areas(i) = tri_area(pts(1,:), pts(2,:), pts(3,:));
        pmids(i,:) = sum(pts,1)/3;
    end
    pmid = [0, 0];
    pmid(1) = sum( pmids(:,1) .* areas ) / sum(areas);
    pmid(2) = sum( pmids(:,2) .* areas ) / sum(areas);
end