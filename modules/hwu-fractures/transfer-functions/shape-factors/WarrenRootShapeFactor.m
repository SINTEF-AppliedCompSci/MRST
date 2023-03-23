classdef WarrenRootShapeFactor < ShapeFactor
    %WARREN_ROOT_SHAPE_FACTOR
    %
    % calculates the shape factor handling anisotropic permeability
    % dimensions are k.L^-2
    % when the permeability is isotropic this formula collapses to
    % km.sigma as quoted in the literature
    %
    %Warren and Root 1963
    % sigma = 4n(n+1)/l^2
    % where n is the sets of normal fractures
    % here n is assumed to be 1

    properties

    end

    methods
        function shapefactor = WarrenRootShapeFactor(block_dimension)
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
                'lx!=ly!=lz Warren&Root shape factor assumes equal matrix block lengths')
            
            kavg = (1/3)*(kx + ky + kz);

            sigma = 60*kavg./lx.^2;
        end
    end
end
