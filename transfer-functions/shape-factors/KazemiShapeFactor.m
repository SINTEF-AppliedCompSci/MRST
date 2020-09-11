classdef KazemiShapeFactor < ShapeFactor
    %KAZEMI_SHAPE_FACTOR
    %
    % calculates the shape factor handling anisotropic permeability
    % dimensions are k.L^-2
    % when the permeability is isotropic this formula collapses to
    % km.sigma as quoted in the literature

    properties

    end

    methods

        function shapefactor = KazemiShapeFactor(block_dimension)
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

            sigma = 4*(kx./(lx.^2)+ky./(ly.^2)+kz./(lz.^2));
        end
    end
end
