function [G, rock, fluid, state0, schedule] = setupSaigupBC(varargin)
    
    opt = struct('layers', [], ...
                 'time', 20*year, ...
                 'dt', 30*day);
    opt = merge_options(opt, varargin{:});

    try
       grdecl = readGRDECL(fullfile(ROOTDIR, 'examples', 'data', ...
                                    'SAIGUP', 'SAIGUP.GRDECL'));
       grdecl = convertInputUnits(grdecl, getUnitSystem('METRIC'));
       G      = processGRDECL(grdecl);
       G      = computeGeometry(G);
       rock   = grdecl2Rock(grdecl, G.cells.indexMap);
    catch me
       error('SAIGUP model data is not available.')
    end

    if ~isempty(opt.layers)
        ijk = gridLogicalIndices(G);
        c = ismember(ijk{3}, opt.layers);
        G = extractSubgrid(G, c);
        G = computeGeometry(G);
        rock = makeRock(G, rock.perm(c,:), rock.poro(c));
    end
    
    G = computeCellDimensions2(G);
    fluid = initSimpleADIFluid('phases', 'WO'             , ...
                               'n'     , [1,1]            , ...
                               'mu'    , [1,1]*centi*poise, ...
                               'rho'   , [1,1]* kilogram/meter^3, ...
                               'c'     , [1e-6,1e-6]/barsa);
    
    bc = [];
    pInj = 1000*barsa;
    fw = find(G.faces.centroids(:,2) == 9000);
    bc = addBC(bc, fw, 'pressure', pInj, 'sat', [1,0]);
    pProd = 10*barsa;
    fe = find(G.faces.centroids(:,2) == 0);
    bc = addBC(bc, fe, 'pressure', pProd, 'sat', [1,0]);

    dt = rampupTimesteps(opt.time, opt.dt);
    schedule = simpleSchedule(dt, 'bc', bc);
    state0 = initResSol(G, pProd, [0,1]);
end