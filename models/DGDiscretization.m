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
            
            disc.basis = DGBasisFunctions(disc.G, disc.degree);
            
%             disc.basis   = dgBasis(disc.degree, dim, disc.basis);
            disc.limiter = dgLimiter(disc     , disc.limiter);
            
        end
        
        %-----------------------------------------------------------------%
        function [xhat, translation, scaling] = transformCoords(disc, x, cells)
            
            G           = disc.G;
            translation = -G.cells.centroids(cells);
            scaling     = 1./(G.cells.diameters(cells)/(2*sqrt(G.griddim)));
            xhat        = (x + translation)./scaling;
               
        end
        
%         function accumulation = getAccumulationTerm(disc, sdof, sdof0)
%             
%             
%             sW = getSatFromDof(x, c, sdof
%             
%         end
            
    end
    
end