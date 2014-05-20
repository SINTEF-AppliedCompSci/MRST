function result = makeResultStructure(states, G, rock, EOSCO2, rhoW, tfun)

    if G.griddim == 2
        result = case2D(states, G);
    else 
        assert(G.griddim == 3);
        result = case3D(states, G, rock, EOSCO2, rhoW, tfun);
    end
end

% ----------------------------------------------------------------------------

function result = case2D(states, Gt)
    for i = 1:numel(states)
        result(i).h         = states{i}.h;
        result(i).pI        = states{i}.pressure;
        result(i).tI        = states{i}.extras.tI;
        result(i).rhoI      = states{i}.extras.rhoI;
        result(i).pTop      = states{i}.extras.pTop;
        result(i).tTop      = states{i}.extras.tTop;
        result(i).rhoTop    = states{i}.extras.rhoTop;
        result(i).fluxCO2   = states{i}.extras.fluxCO2;
        result(i).fluxBrine = states{i}.extras.fluxBrine;
    end        
end

% ----------------------------------------------------------------------------

function result = case3D(states, G, rock, EOSCO2, rhoW, tfun)
    
    Gt = topSurfaceGrid(G);

    for i = 1:numel(states)
        result(i).h         = sat2height(states{i}.s(:,2), Gt, rock);
        result(i).pI        = interfacePressure(result(i).h, states{i}, Gt, ...
                                                EOSCO2, rhoW, tfun);
        result(i).tI        = tfun(Gt.cells.z + result(i).h);
        result(i).rhoI      = EOSCO2.rho(result(i).pI, result(i).tI);
        result(i).pTop      = interfacePressure(zeros(Gt.cells.num,1), states{i}, Gt, ... 
                                                EOSCO2, rhoW, tfun);
        result(i).tTop      = tfun(Gt.cells.z);
        result(i).rhoTop    = EOSCO2.rho(result(i).pTop, result(i).tTop);
        result(i).fluxCO2   = []; % atm, flux not supported 
        result(i).fluxBrine = []; 
    end    
end

% ----------------------------------------------------------------------------

function pI = interfacePressure(h, state, Gt, EOSCO2, rhoW, tfun)

    icells = interfaceCells(h, state, Gt);
    idepth = Gt.cells.z + h;
    rhoCO2 = EOSCO2.rho(state.pressure(icells), tfun(idepth));
    
    dz = idepth - Gt.parent.cells.centroids(icells, 3);
    
    pI = state.pressure(icells);
    below_ix = dz>0;
    pI( below_ix) = pI( below_ix) + rhoCO2(below_ix) * norm(gravity) .* dz(below_ix);
    pI(~below_ix) = pI(~below_ix) + rhoW * norm(gravity) .* dz(~below_ix);
    
end

% ----------------------------------------------------------------------------

function icells = interfaceCells(h, state, Gt)
% Determine the cells in the 3D grid that contain the liquid/gas interface
% (when saturation has been transformed into depth, assuming sharp interface)    
    cells3D = Gt.columns.cells;
    z_vals  = Gt.columns.z;
    colpos  = Gt.cells.columnPos;
    
    depths = arrayfun(@(i) z_vals(colpos(i): colpos(i+1)-1), ...
                      1:Gt.cells.num, 'uniformoutput', false);
    
    colcell_ix = arrayfun(@(i) find(depths{i} >= h(i), 1, 'first'), 1:Gt.cells.num);
    
    icells = Gt.columns.cells(colpos(1:end-1) + colcell_ix' - 1);
    
end

% ----------------------------------------------------------------------------

% function val = interfacePressure(state, Gt)
%     tol = 1e-4;
%     cells = Gt.columns.cells;
%     colpos = Gt.cells.columnPos;
%     % identify the cells in the 3D grid that contain the gas/liquid interface
%     iface_colcell = ...
%         arrayfun(@(i) sum(state.s(cells(colpos(i):colpos(i+1)-1),2) > tol, 1), ...
%                  1:Gt.cells.num);
%     iface_colcell(iface_colcell == 0) = 1;  % we don't want zero indices
    
%     iface_cells = cells(colpos(1:end-1) + iface_colcell' - 1);
    
% end
