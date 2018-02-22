classdef DGDiscretization < WENODiscretization
    
    properties
        degree
        basis
        limiter
    end
    
    methods
        %-----------------------------------------------------------------%
        function disc = DGDiscretization(model, dim, varargin)
            
            disc = disc@WENODiscretization(model, dim);
            
            disc.degree  = 1;
            disc.basis   = 'legendre';
            disc.limiter = 'tvb';
            disc         = merge_options(disc, varargin{:});
            
            disc.basis   = dgBasis(disc.degree, dim, disc.basis);
            disc.limiter = dgLimiter(disc     , disc.limiter);
            
        end
        
        %-----------------------------------------------------------------%
        function xhat = scaling(disc, x, cells)
            
            G = disc.G;
            xhat = (x - G.cells.centroids(cells))...
                   ./(G.cells.diameters(cells)/(2*sqrt(G.griddim)));
               
        end
            
    end
    
end