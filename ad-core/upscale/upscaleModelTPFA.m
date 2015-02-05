function model = upscaleModelTPFA(model, partition, varargin)
    opt = struct('validatePartition', true,...
                 'transCoarse',       [], ...
                 'permCoarse',        [], ...
                 'neighborship',      [], ...
                 'poroCoarse',        []);
     
    opt = merge_options(opt, varargin{:});
    
    % Handle grid
    mrstModule add coarsegrid
    
    G = model.G;
    rock = model.rock;
    
    CG = getGrid(G, partition, opt);
    
    rock_c = getRock(rock, CG, opt);
    
    [Tc, Nc] = getTransmissibility(CG, rock_c, opt);
    
    model.G = CG;
    model.rock = rock_c;
    model = model.setupOperators(CG, rock_c, 'neighbors', Nc, 'trans', Tc);
end

function CG = getGrid(G, partition, opt)
    if opt.validatePartition
        partition = processPartition(G, partition);
        partition = compressPartition(partition);
    end
    CG = generateCoarseGrid(G, partition);
    CG = coarsenGeometry(CG);
end

function rock_c = getRock(rock, CG, opt)
    % Handle rock
    poro_c = opt.poroCoarse;
    perm_c = opt.permCoarse;
    
    p = CG.partition;
    cellcount  = accumarray(p, 1);
    if isfield(rock, 'poro') && isempty(poro_c)
        poro = rock.poro;
        % Include NTG directly into porosity for the time being...
        if isfield(rock, 'ntg');
            poro = poro.*rock.ntg;
        end
        % Sum up pore volumes - we want upscaled model to have same porosity.
        poro_c = accumarray(p, poro)./cellcount;
    end
    
    if isfield(rock, 'perm') && isempty(perm_c)
        nK = size(rock.perm, 2);
        perm_c = zeros(CG.cells.num, nK);
        for i = 1:nK
            perm_c = accumarray(p, 1./rock.perm(:, i));
        end
        perm_c = 1./bsxfun(@rdivide, perm_c, cellcount);
    end
    rock_c = makeRock(CG, perm_c, poro_c);
end

function [Tc, N_int] = getTransmissibility(CG, rock_c, opt)
    N = opt.neighborship;
    if isempty(N)
        N = CG.faces.neighbors;
    end
    
    % Handle transmissibility
    Tc = opt.transCoarse;
    nT = numel(Tc);
    
    % No. half faces
    nHf = size(CG.cells.faces, 1);
    % Number of faces
    nF  = CG.faces.num;
    % No. interfaces
    intx = all(N ~= 0, 2);
    nIF = sum(intx);
    
    if nT == 0 || nT == nF
        % We either have no transmissibility given, or it is in the format
        % one entry per face. Either case is fine for setupOperators.
    elseif nT == nHf
        Tc  = 1 ./ accumarray(CG.cells.faces(:,1), 1./Tc, [nF, 1]);
    elseif nT == nIF
        % Use upscaled perm to compute trans for boundary and the provided
        % trans for internal faces.
        Tc_dummy = computeTrans(CG, rock_c);
        Tc_intf = Tc;
        
        Tc = zeros(nF, 1);
        Tc( intx) = Tc_intf(intx);
        Tc(~intx) = Tc_dummy(~intx);
    else
        msg = ['Number of input transmissibility entries (', num2str(nT), ...
            ') does not match either half face count (' num2str(nHf) '),', ...
            ' interface count (', num2str(nIF), ') or facecount (', num2str(nF), ').'];
        error(msg);
    end
    N_int = N(intx, :);
end
