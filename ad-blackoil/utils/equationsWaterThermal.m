function [problem, state] = equationsWaterThermal(state0, state, model, dt, drivingForces, varargin)
% Get linearized problem for single phase water system with black oil-style
% properties

opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'resOnly', false,...
             'iteration', -1);  % Compatibility only

opt = merge_options(opt, varargin{:});

% W = drivingForces.W;
%assert(isempty(drivingForces.bc) && isempty(drivingForces.src))
assert(isempty(drivingForces.src))


s = model.operators;
G = model.G;
f = model.fluid;

[p, T, wellSol] = model.getProps(state, 'pressure','temperature', 'wellsol');

[p0 , T0, wellSol0] = model.getProps(state0, 'pressure','temperature','wellsol');

%Initialization of independent variables ----------------------------------
[wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);

if ~opt.resOnly
    if ~opt.reverseMode
        % define primary varible x and initialize
        [p, T, wellVars{:}] = initVariablesADI(p, T, wellVars{:});
    else
        wellVars0 = model.FacilityModel.getAllPrimaryVariables(wellSol0);
        [p0, T0, wellVars0{:}] = initVariablesADI(p0, T0, wellVars0{:}); %#ok
    end
end
gdz   = s.Grad(G.cells.centroids) * model.gravity';
%--------------------
%check for p-dependent tran mult:
trMult = 1;
if isfield(f, 'tranMultR'), trMult = f.tranMultR(p); end

%check for p-dependent porv mult:
pvMult = 1; pvMult0 = 1;
if isfield(f, 'pvMultR')
    pvMult =  f.pvMultR(p);
    pvMult0 = f.pvMultR(p0);
end
transMult=1;
if isfield(f, 'transMult')
   transMult=f.transMult(p); 
end

trans=s.T.*transMult;
% -------------------------------------------------------------------------
% water props (calculated at oil pressure OK?)
bW     = f.bW(p, T);bW0 = f.bW(p0,T0);
rhoW   = bW.*f.rhoWS;
% rhoW on face, avarge of neighboring cells (E100, not E300)
rhoWf  = s.faceAvg(rhoW);
mobW   = trMult./f.muW(p,T);
dpW     = s.Grad(p) - rhoWf.*gdz;
% water upstream-index
upcw = (value(dpW)<=0);
vW = - s.faceUpstr(upcw, mobW).*trans.*dpW;
bWvW = s.faceUpstr(upcw, bW).*vW;

% get themal propeties
uR = T.*model.rock.cR.*model.rock.rhoR;uR0 = T0.*model.rock.cR.*model.rock.rhoR;
uW = f.uW(p,T);uW0 = f.uW(p0,T0);
hW = f.hW(p,T);

if model.outputFluxes
    state = model.storeFluxes(state, vW, [], []);
end

if model.extraStateOutput
    state = model.storebfactors(state, bW, [], []);
    state = model.storeMobilities(state, mobW, [], []);
    state = model.storeUpstreamIndices(state, upcw, [], []);
end
% EQUATIONS ---------------------------------------------------------------
% water:
primaryVars = {'pressure', 'T', wellVarNames{:}};


eqs{1} = (s.pv/dt).*( pvMult.*bW - pvMult0.*bW0 ) + s.Div(bWvW);

vol=model.G.cells.volumes;
vQ = -s.T_r .* s.Grad(T);
eqs{2} = (1./dt).*((vol-pvMult.*s.pv).*uR-(vol-pvMult0.*s.pv).*uR0) + s.Div( vQ);
% accumulation of energy and advection of enthalpy in fluid, conduction of
% heat is neglected
eqs{2}  =  eqs{2} + ((s.pv/dt).*( pvMult.*uW.*f.rhoWS.*bW - pvMult0.*uW0.*f.rhoWS.*bW0 )...
        +  s.Div( s.faceUpstr(bWvW>0, f.rhoWS.*hW) .* bWvW));
names = {'water', 'temperature'};
types = {'cell', 'cell'};


sW = ones(model.G.cells.num, 1);
[eqs, state, src] = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
                                                                 {p}, {sW}, {mobW}, {rhoW}, ...
                                                                 {}, {}, ...
                                                                 drivingForces);

% [eqs, qBC, qRes, BCtocellMap, qSRC, srcCells, bcCells] = addFluxesFromSourcesAndBC(model, eqs, ...
%                                           {p},...
%                                         {rhoW},...
%                                         {mobW}, ...
%                                        {ones(G.cells.num,1)}, ...
%                                        drivingForces);
                                   
% accumulation of energy and conduction of heat in rock
% add energy from pressure boundaries assume outside has sime enthapy as
%
% if(~isempty(bcCells))
%     qBC=qBC{1};
%     hWbc=hW(bcCells);%.*BCtocellMap;
%     hWbc(qBC>0)=drivingForces.bc.hW(qBC>0);
%     eqs{2} = eqs{2} - BCtocellMap*(hWbc.*f.rhoWS.*qBC);
% end

bcCells = src.bc.sourceCells;
if(~isempty(bcCells))
    qBC=src.bc.phaseMass{1};
    hWbc=hW(bcCells);%.*BCtocellMap;
    hWbc(qBC>0)=drivingForces.bc.hW(qBC>0);
    
    qtbc = src.bc.mapping*(hWbc.*qBC);
    if isempty(src.bc.mapping)
        qtbc = src.bc.mapping*qtbc;
    end
    eqs{2}(bcCells) = eqs{2}(bcCells) - qtbc;
end
srcCells = src.src.sourceCells;
if(~isempty(srcCells))
    qSRC=src.bc.phaseMass{1};
    hWsrc=hW(srcCells);
    hWsrc(srcCells(qSRC>0))=drivingForces.src.hW;
    eqs{2}(srcCells) = eqs{2}(srcCells) - qSRC.*hWsrc;
end

% add energy from conduction from boundary
if(~isempty(drivingForces.bcT))
    assert(all(strcmp(drivingForces.bcT.type,'temperature')));
    bc_t=strcmp(drivingForces.bcT.type,'temperature');
    if(any(bc_t))
        %assert(isempty(drivingForces.bc),'can only have temprature boundary with nowflow');
        bc=drivingForces.bcT;
        [bQqQbc,bcT_cell] = temperatureBCContrib(G,s,T,bc);
        eqs{2}(bcT_cell)  = eqs{2}(bcT_cell) + bQqQbc;
    end
end
                                   
                                   
% well equations
%sat = {sW, sO, sG};
components = {};
[eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, names, types, wellSol0, wellSol, wellVars, wellMap,...
    p, {mobW}, {rhoW}, hW, components, dt, opt);
problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
end
%--------------------------------------------------------------------------

function  [bQqQbc,bc_cell] = temperatureBCContrib(G,s,T,bc)        
        assert(all(strcmp(bc.type,'temperature')),'only pressure bc allowed');
        Tr_bc=s.T_r_all(bc.face);
        assert(all(sum(G.faces.neighbors(bc.face,:)>0,2)==1),'bc on internal boundary');
        bc_cell=sum(G.faces.neighbors(bc.face,:),2);        
        Tbc=T(bc_cell);
        dTbc=bc.value-Tbc;
        bQqQbc  = -Tr_bc.*(dTbc);
end

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




