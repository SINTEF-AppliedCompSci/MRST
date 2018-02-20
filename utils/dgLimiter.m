function limiter = dgLimiter(disc, type, varargin)

    opt = struct('threshold', 0, 'innerType', 'minmod');
    opt = merge_options(opt, varargin{:});

    G = disc.G;

    switch type
        case 'tvb'
            lim = getLimiter(opt.innerType);
            
            h     = G.cells.diameters;
            
            faces = G.cells.faces(:,1);
            isbf  = any(G.faces.neighbors(faces,:) == 0,2);
            faces = faces(~isbf);
            
            sW     = @(x,dof,c) getSatFromDof(x, c, dof, disc);
            
            cells = rldecode((1:G.cells.num)', diff(G.cells.facePos), 1);
            nbf   = accumarray(cells, isbf);
            cells = rldecode((1:G.cells.num)', diff(G.cells.facePos) - nbf, 1);
            
            xf = (G.faces.centroids(faces,:) - G.cells.centroids(cells,:))./h(cells);
            
            sWjump = @(dof) abs(sW(xf, dof, G.faces.neighbors(faces,1)) ...
                              - sW(xf, dof, G.faces.neighbors(faces,2))  );
            indicator = @(dof) accumarray(cells, sWjump(dof) > opt.threshold) > 0;

            limiter = @(dof) approx_grad(dof, disc);
            
%             limiter = @(dof) tvbLimiter(disc, dof, indicator, opt);
            

            ll = getLimiter(opt.innerType);
            
%             model.operators
            
%             limiter = @(dof,c) 2/sqrt(G.*lim(dof) > M;
    end

end

function sigma = approx_grad(dof, disc)
    ind = 1:disc.basis.nDof:disc.G.cells.num*disc.basis.nDof;
    q = dof(ind);
    sigma = cell(1, disc.dim);
    [sigma{:}] = deal(0);
    for d = 1:disc.dim
        b = disc.interp_setup.tri_basis{d};
        for l = 1:disc.dim+1
            loc_cells = disc.interp_setup.C(:, l);
            ccl = disc.interp_setup.tri_cells(loc_cells);
            ds = b(:, l).*q(ccl);

            sigma{d} = sigma{d} + ds;
        end
    end
end

function [newdof, flag] = tvbLimiter(model, dof, indicator, opt)
    
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
    flag   = indicator(dof);
    newdof = 1;

end

function limiter = getLimiter(type)
    switch type
        case 'minmod'
            limiter = @(val) max(val)*all(val < 0) + min(val)*all(val > 0);
    end
end