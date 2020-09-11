classdef CoatsShapeFactor < ShapeFactor
    %COATS_SHAPE_FACTOR
    %
    % calculates the shape factor handling anisotropic permeability
    % dimensions are k.L^-2
    % when the permeability is isotropic this formula collapses to
    % km.sigma as quoted in the literature

    properties

    end

    methods

        function shapefactor = CoatsShapeFactor(block_dimension)
            has_permeability = 1;
            shapefactor = shapefactor@ShapeFactor(block_dimension,...
                                                  has_permeability);
        end

        function [sigma] = calculateShapeFactor(sf, rock)
            
            kx = rock.perm(:,1);
            try
                ky = rock.perm(:,2);
                kz = rock.perm(:,3);
            catch
                ky = rock.perm(:,1);
                kz = rock.perm(:,1);
            end
            lx = sf.block_dimension(:,1);
            ly = sf.block_dimension(:,2);
            lz = sf.block_dimension(:,3);
            assert(lx==ly && lx==lz,...
                'lx!=ly!=lz Coats shape factor assumes equal matrix block lengths')
            
            kavg = (1/3)*(kx + ky + kz);

            sigma = 24*kavg./lx.^2;

        end
    end
end
