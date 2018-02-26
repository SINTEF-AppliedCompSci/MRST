function limiter = dgLimiter(disc, type, varargin)

    opt = struct('threshold', 0, 'innerType', 'minmod');
    opt = merge_options(opt, varargin{:});

    G = disc.G;

    switch type
        case 'tvb'
            
            h     = G.cells.diameters;
            
            faces = G.cells.faces(:,1);
            isbf  = any(G.faces.neighbors(faces,:) == 0,2);
            faces = faces(~isbf);
            
            sW     = @(x,dof,c) disc.evaluateSaturation(x, c, dof);
            
            cells = rldecode((1:G.cells.num)', diff(G.cells.facePos), 1);
            nbf   = accumarray(cells, isbf);
            cells = rldecode((1:G.cells.num)', diff(G.cells.facePos) - nbf, 1);
            
            xf = G.faces.centroids(faces,:);
            
            c_l  = G.faces.neighbors(faces,1);
            xf_l = disc.transformCoords(xf, c_l);
            c_r  = G.faces.neighbors(faces,2);
            xf_r = disc.transformCoords(xf, c_r);
            
%             xf = (G.faces.centroids(faces,:) - G.cells.centroids(cells,:))./(h(cells)/(2*sqrt(G.griddim)));
            
            sWjump = @(dof) abs(sW(xf_l, dof, c_l) ...
                              - sW(xf_r, dof, c_r)  );
            indicator = @(dof) accumarray(cells, sWjump(dof) > opt.threshold) > 0;
            
            limiter = @(dof) tvbLimiter(disc, dof, indicator);
            
    end

end

function dofbar = approx_grad(dof, disc)

    ind  = 1:disc.basis.nDof:disc.G.cells.num*disc.basis.nDof;
    q    = dof(ind);
    G    = disc.G;
    nDof = disc.basis.nDof;
    
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
    
    dofbar = zeros(G.cells.num, disc.dim);
    for cNo = 1:G.cells.num
        
%         gradix = gradpos(cNo):gradpos(cNo+1)-1;
        gradix = disc.interp_setup.cell_support{cNo};
        for dNo = 1:disc.dim
            dofix = (cNo-1)*nDof + 1 + dNo;
            ss = sigma{dNo}(gradix);
            ss = ss(ss~=0);
            val = [dof(dofix); ss];
            dofbar(cNo, dNo) = minmod(val);
        end

    end
    
end

function newdof = tvbLimiter(disc, dof, indicator)

    flag   = indicator(dof);
    
    dofbar = approx_grad(dof, disc)';
    dofbar = dofbar(:);
    
    newdof = dof;
    
    nDof = disc.basis.nDof;
    
    ix1 = mcolon((find(flag)-1)*nDof + 1 + 1, (find(flag)-1)*nDof + 1 + disc.dim);
    ix2 = mcolon((find(flag)-1)*disc.dim + 1, (find(flag)-1)*disc.dim + disc.dim);
    newdof(ix1) = dofbar(ix2);
    

%     dx = G.cells.centroids(G.faces.neighbors(:,1), :) - G.cells.centroids(G.faces.neighbors(:,2), :);
% 
%     ll = getLimiter(opt.innerType);
%     
%     
%     
%     
% %     dof0 = dof(1:nDof:G.cells.num*nDof);
% %     for dNo = 1:G.griddim
% %         dofx = dof((dNo+1):nDof:G.cells.num*nDof);
% %         for cNo = 1:G.cells.num
% %             dofx(c)
% %             ll([dofx, theta*dof]);
% %         end
%     end


end

function v = minmod(val)
    v = max(val)*all(val < 0) + min(val)*all(val > 0);
end
        
% 
% function limiter = getLimiter(type)
%     switch type
%         case 'minmod'
%             limiter = @(val) max(val)*all(val < 0) + min(val)*all(val > 0);
%     end
% end