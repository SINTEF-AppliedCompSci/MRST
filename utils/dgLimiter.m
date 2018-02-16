function limiter = dgLimiter(model, type, varargin)

    opt = struct('threshold', 0, 'innerType', 'minmod');
    opt = merge_options(opt, varargin{:});

    G = model.G;

    switch type
        case 'tvb'
            lim = getLimiter(opt.innerType);
            
            h     = G.cells.diameters;
            
            faces = G.cells.faces(:,1);
            isbf  = any(G.faces.neighbors(faces,:) == 0,2);
            faces = faces(~isbf);
            
            sW     = @(x,dof,c) getSatFromDof(x, c, dof, model);
            
            cells = rldecode((1:G.cells.num)', diff(G.cells.facePos), 1);
            nbf = accumarray(cells, isbf);
            cells = rldecode((1:G.cells.num)', diff(G.cells.facePos) - nbf, 1);
            
            
            xf = (G.faces.centroids(faces,:) - G.cells.centroids(cells,:))./h(cells);
            
            sWjump = @(dof) abs(sW(xf, dof, G.faces.neighbors(faces,1)) ...
                              - sW(xf, dof, G.faces.neighbors(faces,2))  );
            indicator = @(dof) accumarray(cells, sWjump(dof) > opt.threshold) > 0;

            limiter = @(dof) tvbLimiter(model, dof, indicator, opt);
            
            op = model.operators;
            
            dx = G.cells.centroids(op.N(:,2), :) - ...
                 G.cells.centroids(op.N(:,1), :);
             
            [C, pts, grad_basis, supports, linear_weights, scaling] = getTriangulation(model);
             
            dMat = cell(G.griddim,1);
            for dNo = 1:G.griddim
                dMat{dNo} = sparse(op.N(:,1), op.N(:,2), -abs(dx(:,dNo)), G.cells.num, G.cells.num);
            end
             
            tol = 1e-3;
            connections = cell(G.griddim,1);
            for dNo = 1:G.griddim
                
                [cr, cc] = find(dMat{dNo} == min(dMat{dNo},[], 2));
                
                connections{dNo} = repmat((1:G.cells.num)', 1, 2);
                connections{dNo}(op.N(:,1),1) = op.N(:, 2);
                
%                 .*(dx(:,dNo) > 0);
                
                
                connections{dNo}(op.N(:,2),2) = op.N(: ,2).*(dx(:,dNo) < 0);
            end

            ll = getLimiter(opt.innerType);
            
%             model.operators
            
%             limiter = @(dof,c) 2/sqrt(G.*lim(dof) > M;
    end

end

function limiter = tvbLimiter(model, dof, indicator, opt)
    
    dx = G.cells.centroids(G.faces.neighbors(:,1), :) - G.cells.centroids(G.faces.neighbors(:,2), :);

    ll = getLimiter(opt.innerType);
    
    
    
    
%     dof0 = dof(1:nDof:G.cells.num*nDof);
%     for dNo = 1:G.griddim
%         dofx = dof((dNo+1):nDof:G.cells.num*nDof);
%         for cNo = 1:G.cells.num
%             dofx(c)
%             ll([dofx, theta*dof]);
%         end
%     end
    limiter = indicator(dof);

end

function limiter = getLimiter(type)
    switch type
        case 'minmod'
            limiter = @(val) max(val)*all(val < 0) + min(val)*all(val > 0);
    end
end