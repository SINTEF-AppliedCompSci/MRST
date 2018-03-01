function [limiter, sWjump, indicator] = dgLimiter(disc, type, varargin)

    opt = struct('threshold', 0.5, 'innerType', 'minmod');
    opt = merge_options(opt, varargin{:});

    G = disc.G;

    switch type
        case 'tvb'
            
            faces = G.cells.faces(:,1);
            isbf  = any(G.faces.neighbors(faces,:) == 0,2);
            faces = faces(~isbf);
            
            sW     = @(x, c, dof) disc.evaluateSaturation(x, c, dof);
            
            cells = rldecode((1:G.cells.num)', diff(G.cells.facePos), 1);
            nbf   = accumarray(cells, isbf);
            cells = rldecode((1:G.cells.num)', diff(G.cells.facePos) - nbf, 1);
            
            xf = G.faces.centroids(faces,:);
            
            c_l  = G.faces.neighbors(faces,1);
            xf_l = disc.transformCoords(xf, c_l);
            c_r  = G.faces.neighbors(faces,2);
            xf_r = disc.transformCoords(xf, c_r);
            
%             xf = (G.faces.centroids(faces,:) - G.cells.centroids(cells,:))./(h(cells)/(2*sqrt(G.griddim)));
            
            sWjump = @(dof) abs(sW(xf_l, c_l, dof) ...
                              - sW(xf_r, c_r, dof)  );
            indicator = @(dof) accumarray(cells, sWjump(dof) > opt.threshold) > 0;
            
            limiter = @(dof) tvbLimiter(disc, dof, indicator);
            
        case 'cap'
            
            faces = G.cells.faces(:,1);
            nodes = G.faces.nodes(mcolon(G.faces.nodePos(faces), G.faces.nodePos(faces+1)-1));
%             nodes = reshape(nodes, 2, [])';
%             swap  = G.faces.neighors(faces,1) ~= rldecode((1:G.cells.num)', diff(G.cells.facePos), 1);
%             nodes(swap,:) = nodes(swap,[2,1]); nodes = nodes(:,1);


            
            nfn = diff(G.faces.nodePos);
            ncn = accumarray(rldecode((1:G.cells.num)', diff(G.cells.facePos), 1), nfn(faces));
            
            x = G.nodes.coords(nodes,:);
            
            cells = rldecode((1:G.cells.num)', ncn, 1);
            x = disc.transformCoords(x, cells);
            
            s = @(dof) disc.evaluateSaturation(x, cells, dof);
            indicator = @(dof) accumarray(cells, s(dof) < 0 | s(dof) > 1) > 0;
            
            sWjump = @(dof) 0;
            limiter = @(dof) false;
            
    end

end

function dofbar = approx_grad(dof, disc)

    ind  = 1:disc.basis.nDof:disc.G.cells.num*disc.basis.nDof;
    q    = dof(ind);
    G    = disc.G;
    nDof = disc.basis.nDof;
    
    sigma = cell(1, disc.dim);
    [sigma{:}] = deal(0);
    
    if disc.dim == 1
        left  = disc.interp_setup.cellNeighbors(:,1);
        right = disc.interp_setup.cellNeighbors(:,2);
        left(left == 0) = 1;
        right(right == 0) = G.cells.num;
        cells = 1:G.cells.num;
        dx    = G.cells.dx(:,1);
        for d = 1:disc.G.griddim
            if d == 1
%                 dd = [(q(cells) - q(left) )./(dx(cells)/2 + dx(left)/2 ), ...
%                             (q(right) - q(cells))./(dx(cells)/2 + dx(right)/2)];
                dd = [(q(cells) - q(left) ), (q(right) - q(cells))];
                sigma{d} = reshape(dd', [], 1);
                sigma{d} = sigma{d}(2:end-1);
            else
                sigma{d} = zeros(2*G.cells.num,1);
                sigma{d} = sigma{d}(2:end-1);
            end
        end
        
        disc.interp_setup.cell_support = cell(G.cells.num,1);
        for cNo = 1:G.cells.num
            if cNo == 1
                cs = 1;
            elseif cNo < G.cells.num
                cs = (((cNo-1)*2+1):cNo*2) - 1;
            else
                cs = G.cells.num*2-2;
            end
            disc.interp_setup.cell_support{cNo} = cs;
        end
    
    else
    
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
    
    dofbar = zeros(G.cells.num, disc.G.griddim);
    for cNo = 1:G.cells.num
        
%         gradix = gradpos(cNo):gradpos(cNo+1)-1;
        gradix = disc.interp_setup.cell_support{cNo};
        for dNo = 1:G.griddim
            dofix = (cNo-1)*nDof + 1 + dNo;
            ss = sigma{dNo}(gradix);
%             ss = ss(ss~=0);
            val = [dof(dofix); 2*ss];
            db = minmod(val);
            dofbar(cNo, dNo) = db;
        end

    end
    
end

function newdof = tvbLimiter(disc, dof, indicator)

    flag   = indicator(dof);
    newdof = dof;
    
    if any(flag)
    
        dofbar = approx_grad(dof, disc)';
        dofbar = dofbar(:);

        nDof = disc.basis.nDof;


        cells = find(flag);

        ix1 = disc.getDofIx(2:disc.G.griddim+1, cells);

    %     ix1 = mcolon((find(flag)-1)*nDof + 1 + 1, (find(flag)-1)*nDof + 1 + disc.G.griddim);
        ix2 = mcolon((find(flag)-1)*disc.G.griddim + 1, (find(flag)-1)*disc.G.griddim + disc.G.griddim);
        newdof(ix1) = dofbar(ix2);
        if disc.degree > 1
            ix = disc.getDofIx((disc.G.griddim + 2):nDof, cells);
    %         ix = mcolon((find(flag)-1)*nDof + 1 + disc.G.griddim + 1, (find(flag)-1)*nDof + 1 + disc.G.griddim + nDof);
            newdof(ix) = 0;
        end
    
    end

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