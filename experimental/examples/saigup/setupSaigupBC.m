function [G, rock, fluid, state0, schedule] = setupSaigupBC(varargin)
%Undocumented Utility Function

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
    
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
