function [div] = VEM2D_div(G, varargin)
% For now define all the primary operators with out units i.e without
% volum, area, length


%{ 
Copyright 2009 - 2014 SINTEF ICT, Applied Mathematics
%} 
    opt = struct('addface_dof', false); 
    opt = merge_options(opt, varargin{:}); 
    div = VEM2d_loc_int(G, opt); 
end

% ---------------------------------------------------------------------------- 

function div = VEM2d_loc_int(G, opt)
    qf     = calculateQF(G); 
    qf     = reshape(qf', [], 1); 
    dofs   = mcolon(G.griddim * (G.cells.nodes - 1) + 1, ...
                    G.griddim * (G.cells.nodes - 1) + G.griddim); 
    ndofs  = G.nodes.num * G.griddim; 
    cellno = rldecode([1:G.cells.num]', diff(G.cells.nodePos) * G.griddim); 
    
    if (opt.addface_dof)
        
        ncf        = diff(G.cells.facePos); 
        fcellno    = rldecode([1:G.cells.num]', ncf); 
        faces      = G.cells.faces(:, 1); 
        fasgn      = (2 * (G.faces.neighbors(faces, 2) == fcellno) - 1) .* ...
            G.faces.areas(faces); 
        cellno     = [cellno; fcellno]; dofs = [dofs'; ndofs + faces]; 
        qf         = [qf; fasgn]; 
        ndofs      = ndofs + G.faces.num; 
        
    end    
    
    div = sparse(cellno, dofs, qf, G.cells.num, ndofs); 
    
end

