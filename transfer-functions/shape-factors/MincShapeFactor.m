classdef MincShapeFactor < ShapeFactor
    %MINC_SHAPE_FACTOR
    %

    properties
        opt
    end

    methods

        function shapefactor = MincShapeFactor(block_dimension, varargin)
            opt = struct('volume_fractions', [], ...
                         'domains', []);
            has_permeability = 1;
            shapefactor = shapefactor@ShapeFactor(block_dimension,...
                                                  has_permeability);
            shapefactor.opt = merge_options(opt, varargin{:});
        end

        function [sigma] = calculateShapeFactor(sf, rock)

            [ dist, area] = shell_properties( sf.opt.volume_fractions, sf.block_dimension );
            
            %TODO anisotropic permeability 
            kx = rock.perm(:,1);
            
            
            if sf.opt.domains(1)==0
                sigma = area(sf.opt.domains(2))*kx/dist(1)/prod(sf.block_dimension);
            else
                sigma = area(sf.opt.domains(2))*kx...
                    /(dist(sf.opt.domains(1)) + dist(sf.opt.domains(2)))...
                    /prod(sf.block_dimension);
            end
        end
    end
end
