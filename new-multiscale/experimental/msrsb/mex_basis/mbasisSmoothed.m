function [I, I_compressed] = mbasisSmoothed(sg, varargin)
    opt = struct('tolerance', 5e-3, ...
                 'maxiter',   1000, ...
                 'omega',     2/3, ...
                 'Interpolator',        [], ...
                 'active',    true(sg.nf, 1));
    opt = merge_options(opt, varargin{:});
    
    
    nc = sg.nc;
    nf = sg.nf;
    
    [cells, cellPos, cellNum] = indirmap(sg.cells);
    % Same dims as cells
    isBnd = vertcat(sg.isBnd{:});
    
    
    [iiAll, iiPosAll, vvAll, diagAll] = deal(cell(nc, 1));
    for i = 1:nc
        [ii, jj, vv] = find(sg.subsys{i});
        
        [~, tmp] = rlencode(jj);
        iiPos = [0; cumsum(tmp)];
        clear tmp
        
        iiAll{i}    = int32(ii - 1);
        iiPosAll{i} = int32(iiPos);
        vvAll{i}    = vv;
        diagAll{i} = diag(sg.subsys{i});
    end
    
    if isempty(opt.Interpolator)
        % Initial guess
        I_compressed = zeros(size(cells));
        for i = 1:nc
            % Add one for zero indices
            pos = (cellPos(i):cellPos(i+1)-1) + 1;
            I_compressed(pos) = ismember(cells(pos), find(sg.partition == i) - 1);
        end
    else
        I_compressed = opt.Interpolator;
    end
    
    
    interactCoarse = rldecode((1:nc)', cellNum);
    
    % Global boundary used for renormalization
    GB_globalCells = int32(sg.globalBoundary);
    [GBLocal, GBLocalMap] = indirmap(sg.globalBoundaryToLocal);
    
    
    GB_globalCells = int32(GB_globalCells);
    GBLocal        = int32(GBLocal);
    GBLocalMap     = int32(GBLocalMap);
    
    nel = int32(numel(cells));
    
    nc = int32(nc);
    nf = int32(nf);
    cells = int32(cells);
    isBnd = int32(isBnd);
    cellPos = int32(cellPos);
    
    tic()
    I_compressed = ...
    mex_iteratedJacobiBasisFaster(...
        nf, nc, nel, I_compressed,...
        cells, isBnd, cellPos, ...
        GB_globalCells, GBLocal, GBLocalMap, ...
        iiAll, iiPosAll, vvAll, diagAll, ...
        int32(opt.active), opt.tolerance, opt.maxiter, opt.omega);
    toc();
    
    bad = isnan(I_compressed);
    if any(bad)
        c = cells+1;
        fb = find(bad);
        warning('NaN found in operator, attempting fix...');
        if mrstVerbose
            disp('Affected cells:')
            disp(unique(c(fb)))
        end
        isSelf = sg.partition(c(bad)) == interactCoarse(bad);
        % Assign nan cells unitary value if they are inside their own
        % coarse block, otherwise 0.
        I_compressed(fb(isSelf)) = 1;
        I_compressed(fb(~isSelf)) = 0;
    end
    
    I = sparse(double(cells) + 1, interactCoarse, I_compressed, double(nf), double(nc));
    
    % This part is a hack...
%     centers = sg.centers + 1;
%     for i = 1:nc
%         if I(centers(i), i) < 0.9
%             I(centers(i), :) = 0;
%             I(centers(i), i) = 1;
%         end
%     end
end

function [v, vpos, num] = indirmap(data)
    % Zero indexing indirection map
    v = vertcat(data{:});
    num = cellfun(@numel, data);
    vpos = [0; cumsum(num)];
end

