function [W, W_ve, bc] = wellSetup(Gt, rock, wellPts, wellamounts, rho, varargin)
    opt = struct('radius', 0.15*meter);
    opt = merge_options(opt, varargin{:});
    
    rock2D = averageRock(rock, Gt);
    
    if numel(wellamounts) == 1
        wellamounts = repmat(wellamounts, size(wellPts, 1), 1);
    end
    
    W = [];
    for i=1:size(wellPts,1)
        
        wpos = wellPts(i,:);
        % Amount in kilograms !!!
        rates = wellamounts(i)*kilogram./(year*rho(1)*kilogram*meter^3);
        dist = sqrt(sum(bsxfun(@minus,Gt.cells.centroids(:,1:2), double(wpos)).^2,2));
        [dd, cellnum] = min(dist); %#ok
        [ix,iy]=ind2sub(Gt.cartDims,Gt.cells.indexMap(cellnum));
        wellIx = double([ix iy]); 
        W      = verticalWell(W, Gt.parent, rock, wellIx(1), wellIx(2), ...
                          1, 'Type', 'rate', 'Val', rates(1), ...
                          'Radius', opt.radius, 'comp_i', [1,0], 'name', ...
                          'I','InnerProduct','ip_tpf');
    end

    % Well in 2D model
    W_ve = convertwellsVE(W, Gt.parent, Gt, rock2D);

    % BC in 2D model
    i = any(Gt.faces.neighbors==0, 2);  % find all outer faces
    I = i(Gt.cells.faces(:,1));         % vector of all faces of all cells, true if outer
    j = false(6,1);                     % mask, cells can at most have 6 faces,
    %j(1:4)=true;
    j(1)=true;
    %   extract east, west, north, south
    J = j(Gt.cells.faces(:,2));         % vector of faces per cell, true if E,W,N,S
    bcIxVE = Gt.cells.faces(I & J, 1);

    bc = addBC([], bcIxVE, 'pressure', ...
                Gt.faces.z(bcIxVE)*rho(2)*norm(gravity));
    bc.sat = zeros(size(bc.face));
    bc.h = zeros(size(bc.face));
end

