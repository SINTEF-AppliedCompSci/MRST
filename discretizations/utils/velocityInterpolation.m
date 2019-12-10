function vi = velocityInterpolation(G, type)
    % Construct velocity vector from face fluxes
    switch type
        case 'mimetic'

            cellNo    = rldecode((1:G.cells.num)', diff(G.cells.facePos), 1);
            faceNo = G.cells.faces(:,1);
            C = G.faces.centroids(faceNo, :) - G.cells.centroids(cellNo, :);
            sgn = 1 - 2*(cellNo ~= G.faces.neighbors(faceNo,1));

            D = cell(1,G.griddim);
            for dNo = 1:G.griddim
                D{dNo} = sparse(cellNo, faceNo, C(:,dNo).*sgn, G.cells.num, G.faces.num);
            end

            f2c = @(v) faceFlux2cellVelocity(D,v);
            vi = struct('D', {D}, 'faceFlux2cellVelocity', f2c);
            
        otherwise
                error('Unknown velocity interpolation type')
    end

end

function vc = faceFlux2cellVelocity(D, v)
    
    vc = cell(1,numel(D));
    for dNo = 1:numel(D)
        vc{dNo} = D{dNo}*v;
    end
    vc = SpatialVector(vc{:});

end