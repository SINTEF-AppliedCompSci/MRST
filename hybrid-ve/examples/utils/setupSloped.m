function [state0, model, schedule, nearWell, setToZero] = setupSloped()
% This is a utility function for setting up the base grid and simulation
% model for the slopedMigration example in the hybrid-ve module.
%
% SYNOPSIS:
%   [state0, model, schedule, nearWell, setToZero] = setupSloped()
%
% REQUIRED PARAMETERS:
%    none
%
% RETURNS:
%
%   state0      - Initial state for sloped simulation.
%
%   model       - Model for sloped simulation.
%
%   schedule    - Schedule for sloped simulation.
%
%   nearWell    - Index of cells near well. 
%
%   setToZero   - Index of faces where transmissibility is set to zero
%                   (sealing faces).
%
% EXAMPLE:
%   slopedMigration
%

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
    uniformRock = true;
    if nargin == 0
        nx = 400;
        ny = 1;
        nz = 50;
    else
        nx = dims(1);
        ny = dims(2);
        nz = dims(3);
    end

    whdist = ceil(nx*0.1);
    wvdist = 1;

    L = 1000;
    G = cartGrid([nx, ny, nz], [L, 1, 100]);
    G0 = computeGeometry(G);

    x = G.nodes.coords(:, 1);
    G.nodes.coords(:, 3) = G.nodes.coords(:, 3) + 1*sin(10*x/50);
    G.nodes.coords(:, 3) = G.nodes.coords(:, 3) + 4*sin(3.5*x/50);
    G.nodes.coords(:, 3) = G.nodes.coords(:, 3) + 8*sin(x/30);

    if 1
        theta = (-2.5/360)*2*pi;
        x = G.nodes.coords(:, 1);
        z = G.nodes.coords(:, 3);
        G.nodes.coords(:, 1) = x*cos(theta) - z*sin(theta);
        G.nodes.coords(:, 3) = x*sin(theta) + z*cos(theta);
    end

    G = computeGeometry(G);

    if uniformRock
        rock = makeRock(G, 100*milli*darcy, 0.3);
        permname = 'uniform';
    else
        rng(0);
        K = logNormLayers(G.cartDims);
        K = K./mean(K);
        K = K.*300*milli*darcy;
        rock = makeRock(G, K, 0.3);

        permname = 'layers';
    end


    [ii, jj, kk] = gridLogicalIndices(G);
% 
%     figure(1);
%     clf
%     plotGrid(G, 'facec', 'none');
%     view(0, 0)


    K = G.cartDims(3);

    setToZero = false(G.faces.num, 1);


        ranges = [200, 800, ceil(5/10*K); ...
                    0, 250, ceil(7.5/10*K); ...
                  750, 1000, ceil(2.5/10*K); ...
              ];
        for i = 1:size(ranges, 1)
            faces = addSealingFaces(G0, 'x_range', ranges(i, 1:2), 'k_range', [ranges(i, 3)-1, ranges(i, 3)]);
            setToZero(faces) = true;
        end

    G = computeGeometry(G);
    c = [1e-4/barsa, 1e-7/barsa];
    fluid = initSimpleADIFluid('phases', 'GW', ...
                               'mu', [6e-2*milli 8e-4]*Pascal*second,...
                               'rho', [760 1200].* kilogram/meter^3, ...
                               'c', c, ...
                               'pRef', 100*barsa);

    time = 1000*year;
    injFrac = 0.03;
    pvFrac = 0.2;

    pv = poreVolume(G, rock);
    inj_rate = pvFrac*sum(pv)/(injFrac*time);
    nt = 100;

    % Specify well information
    W = [];
    bc = [];
    nearWell = false(G.cells.num, 1);

    
    W = verticalWell(W, G, rock, 1, 1, nz, ...
        'type', 'rate', ...  % inject at constant rate
        'val', inj_rate, ... % volumetric injection rate
        'comp_i', [0 1]);    % inject CO2, not water
    
    
    
    bc = pside(bc, G, 'xmax', 100*barsa, 'sat', [1, 0]);
    bc.value = bc.value + fluid.rhoWS.*G.faces.centroids(bc.face, 3)*norm(gravity());
    
    for i = 1:numel(W)
        c = W(i).cells(1);
        hdist = abs(ii - ii(c));
        vdist = abs(kk - kk(c));
        
        nearWell(hdist < whdist & vdist < wvdist) = true;
    end
    if ~isempty(bc)
        nearWell(sum(G.faces.neighbors(bc.face, :), 2)) = true;
    end
    s0 = [1, 0];

   
    dt = [rampupTimesteps(injFrac*time, injFrac*time/nt); repmat((1-injFrac)*time/(2*nt), 2*nt, 1)];
    schedule = simpleSchedule(dt, 'W', W, 'bc', bc);
    schedule.control = [schedule.control; schedule.control];
    schedule.control(2).W(1).status = false;

    schedule.step.control = 1 + (cumsum(schedule.step.val) > injFrac*time);


    T = getFaceTransmissibility(G, rock);
    T(setToZero) = 0;
    model = TwoPhaseWaterGasModel(G, rock, fluid, 1, 1, 'useCNVConvergence', true);
    model.operators.T = T(model.operators.internalConn);
    model.operators.T_all(model.operators.internalConn) = model.operators.T;
    model.extraStateOutput = true;

    state0 = initResSol(G, 0, s0);
end