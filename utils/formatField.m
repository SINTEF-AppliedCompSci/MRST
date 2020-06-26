function u = formatField(u, dim, type)
    switch type
      case {'displacement', 'u'}
        u = reshape(u, dim, [])';
      case {'stress'}        
        uu = reshape(u, dim*dim, [])';
        nc = size(uu, 1);
        switch dim
          case 2
            vdim = 3; % voigt dimension
            u = NaN(nc, vdim);
            u(:, 1) = uu(:, 1);
            u(:, 2) = uu(:, 4);
            u(:, 3) = 0.5*(uu(:, 2) + uu(:, 3));
          case 3
            vdim = 6; % voigt dimension
            u = NaN(nc, vdim);
            u(:, 1) = uu(:, 1);
            u(:, 2) = uu(:, 5);
            u(:, 3) = uu(:, 9);
            u(:, 4) = 0.5*(uu(:, 6) + uu(:, 8));
            u(:, 5) = 0.5*(uu(:, 3) + uu(:, 7));
            u(:, 6) = 0.5*(uu(:, 2) + uu(:, 4));
          otherwise
            error('dimension not valid (accepted : 2D or 3D)');
        end
      otherwise
        error('type not recognized');
    end
end
