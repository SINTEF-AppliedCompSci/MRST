function vi = velocityInterpolation(G, type)

    switch type
        case 'mimetic'

            cellNo    = rldecode((1:G.cells.num)', diff(G.cells.facePos), 1);
            faceNo = G.cells.faces(:,1);
            X = G.faces.centroids(faceNo, :) - ...
            G.cells.centroids(cellNo   , :);
            sgn = 1 - 2*(cellNo ~= G.faces.neighbors(faceNo,1));

            D = cell(1,G.griddim);
            for dNo = 1:G.griddim
                D{dNo} = sparse(cellNo, faceNo, X(:,dNo).*sgn, G.cells.num, G.faces.num)./G.cells.volumes;
            end

            if G.griddim == 2
                faceFlux2CellVelocity = @(v) [D{1}*v, D{2}*v];
            else
                faceFlux2CellVelocity = @(v) [D{1}*v, D{2}*v, D{3}*v];
            end

            vi = struct('D', {D}, 'faceFlux2cellVelocity', faceFlux2CellVelocity);
            
        otherwise
                error('Unknown velocity interpolation type')
    end

end