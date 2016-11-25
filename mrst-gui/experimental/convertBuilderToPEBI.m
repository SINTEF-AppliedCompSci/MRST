function [G,pts,F] = convertBuilderToPEBI(out, n, varargin)
    require upr
    fault = out.faults;
    bdr = out.outline;
    wells = cell(1, size(out.points, 1));
    for i = 1:numel(wells)
        wells{i} = out.points(i, :);
    end
    wells = [wells, out.wells];
		varargin = [varargin, {'polyBdr', bdr}];
		[G,pts,F] = pebiGrid(1/sqrt(n), [1, 1], ...
        'faultLines', fault, 'wellLines', wells, ...
        varargin{:});
end