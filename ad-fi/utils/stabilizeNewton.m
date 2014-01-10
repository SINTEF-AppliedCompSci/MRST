function [dx, meta] = stabilizeNewton(dx, meta, system)
    if ~isfield(meta, 'dx')
        meta.dx = [];
    end
    dxold = meta.dx;
    meta.dx = dx;

    omega = meta.relax;
    switch lower(system.nonlinear.relaxType)
        case 'dampen'
            if omega == 1; return; end;
            dx = cellfun(@(x) x*omega, dx, 'UniformOutput', false);
        case 'sor'
            if omega == 1; return; end;
            for i = 1:numel(dx)
                dx{i} = dx{i}*omega + (1-omega)*dxold{i};
            end
        otherwise
            warning(['Unknown relaxation type ''' system.nonlinear.relaxType '''']);
    end
end
