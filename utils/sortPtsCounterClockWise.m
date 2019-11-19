function t = sortPtsCounterClockWise(p, t)
% Sort the points in counter clockwise order for each element specified by
% the connectivity list t
%   p - 2D point set
%   t - Connectivity list, cell x 1
%
    fTheta = @(x,y)2*pi*double(sign(atan2(y,x))<0) + atan2(y,x);
    for k = 1 : length(t)
        t0 = t{k};
        xy = p(t0,:);
        xy = bsxfun(@minus, xy, mean(xy,1));
        theta = fTheta(xy(:,1), xy(:,2));
        [~, i] = sort(theta);
        t{k} = t0(i);
    end
end