function [pdis, fd] = passToDistmesh(pIB, pOB, multiplier, maxIter, varargin)
% Generate parameters passed to 'distmesh_2d_nwm'
    opt = struct('pIBRadius', []);
    opt = merge_options(opt, varargin{:});

    pIB = [pIB; pIB(1,:)];
    pOB = [pOB; pOB(1,:)];
    
    assert( all(inpolygon(pIB(:,1), pIB(:,2), pOB(:,1), pOB(:,2))), ...
        ['The well region boundary is outside the VOI region, ' ...
        'please enlarge the VOI boundary']) 
    
    fdI = @(p)dpoly(p, pIB);
    fdO = @(p)dpoly(p, pOB);
    fd  = @(p)ddiff(fdO(p), fdI(p));

    lIB = pIB(2:end, :) - pIB(1:end-1, :);
    lIB = sqrt( sum(lIB.^2, 2) );
    lIBave = mean(lIB);
    
    lOB = pOB(2:end, :) - pOB(1:end-1, :);
    lOB = sqrt( sum(lOB.^2, 2) );
    lOBave = mean(lOB);
    
    a = multiplier;
    fh = @(p)min(a*fdI(p)+lIBave, lOBave);
    h0 = lIBave;
 
    pfix = [pIB(1:end-1, :); pOB(1:end-1, :)];

    x_min = min(pOB(:,1));
    x_max = max(pOB(:,1));
    y_min = min(pOB(:,2));
    y_max = max(pOB(:,2));
    bbox = [ [x_min, y_min]; [x_max, y_max] ];

    fprintf('    * Dist Mesh iteration information: \n')
    pdis = distmesh_2d_nwm(fd, fh, h0, bbox, maxIter, pfix, false);
%     close(gcf)

    % Remove points too close to the boundary
    lIBmax = max(lIB);
    lOBmax = max(lOB);

    pdisG = pdis(size(pfix,1)+1:end, :);
    idx = abs(fdI(pdisG)) <= abs(fdO(pdisG));
    
    pdisI = pdisG(idx, :);
    if isempty(opt.pIBRadius)
        idxI = abs(fdI(pdisI)) > lIBmax/2;
    else
        idxI1 = abs(fdI(pdisI)) > 0.1;
        in = nan(size(pdisI,1), size(pIB, 1)-1);
        for i = 1 : size(pIB, 1)-1
            in(:, i) = pointsInCircle(pdisI, pIB(i,:), opt.pIBRadius(i));
        end
        idxI2 = all(~in, 2);
        idxI = idxI1 & idxI2;
    end
    pdisI = pdisI(idxI, :);
    
    pdisO = pdisG(~idx,:);
    idxO = abs(fdO(pdisO)) > lOBmax/2;
    pdisO = pdisO(idxO, :);
    
    pdis = [pfix; pdisI; pdisO];
end

function in = pointsInCircle(p, pc, r)
    d = bsxfun(@minus, p, pc);
    d = sqrt(sum(d.^2, 2));
    in = d <= r;
end