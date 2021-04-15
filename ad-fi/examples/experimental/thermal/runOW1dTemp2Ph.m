% runLayered2D
%require ad-fi deckformat
T_res=300;
p_res=300*barsa;
nx=50;
Lx=2*40;Ly=4*1;Lz=80*1;
G=cartGrid([nx 1 1],[Lx, Ly, Lz]);
G.nodes.coords(:,3)=G.nodes.coords(:,3)+1500;
G=computeGeometry(G);
rock.perm=ones(G.cells.num,1)*100*milli*darcy;
rock.poro=ones(G.cells.num,1)*0.3;
G = computeGeometry(G);

gravity off
%rock  = compressRock(rock, G.cells.indexMap);
fluid_o=initSimpleADIFluid('mu',[0.5 5 1e-1]*centi*poise,'rho',[1000 1000 1000],'n',[2 1 1]);
fluid=fluid_o;
mycase='water_co2_simple'
switch mycase
case 'water_co2'
            ff{1}=h2oProps();
            ff{2}=co2Props();
            ff{3}=h2oProps_simple();
        case 'water_co2_simple'
            ff{1}=h2oProps_simple();
            ff{2}=co2Props_simple();
            ff{3}=h2oProps_simple();
    otherwise    
        error();
end
p_ref=300*barsa;T_ref=300;%kelvin
linear_method='adv'
switch linear_method
        case 'ladv'
            for i=1:3
                fluid_pvt{i} = pvt2Linear(ff{i},p_ref,T_ref);%h2oProps; %water
            end
            fluid=fluid2BlackOil(fluid_o,fluid_pvt,p_ref,T_ref);
        case 'adv'
            for i=1:3;
                fluid_pvt{i} = ff{i};
            end
            fluid=fluid2BlackOilTemp(fluid_o,fluid_pvt,p_ref,T_ref);
        otherwise
            error();
end  
 cR       = 2.17*10.^6;% energy density per volume
fluid.uR = @(T) cR.*T;
if(~isfield(fluid,'uW') || strcmp(linear_method,'ladv'))
    fluid.uW=@(p,T) fluid.hW(p,T)-p./(fluid.rhoWS.*fluid.bW(p));
    fluid.uO=@(p,T) fluid.hO(p,T)-p./(fluid.rhoOS.*fluid.bO(p));
end  
%% presure dependent viscosity
tfac=1;


fluid.relPerm =@(Sw) deal(fluid.krW(Sw),fluid.krW(1-Sw));
extra_fields={};%'BW','BO'}
for k=1:numel(extra_fields)
   fluid=rmfield(fluid,extra_fields{k}); 
end

W = verticalWell([], G, rock,  1,   1, (1:G.cartDims(3)),     ...
                     'Type', 'bhp', 'Val', 200*barsa, ...
                     'Radius', 0.4, 'Name', 'P1','Comp_i',[1 0 0],'sign',-1);

W = verticalWell(W, G,rock,  G.cartDims(1),  G.cartDims(2), (1:G.cartDims(3)),     ...
                     'Type', 'bhp', 'Val', p_res, ...
                     'Radius', 0.4, 'Name', 'I1','Comp_i',[1 0 0],'sign',1); 
               %}
W(1).WI=W(1).WI;
W(2).WI=W(2).WI;

% temprature boundary

%set termeratrue in wells
W(1).T=T_res;
W(2).T=280;%T_res;



W_c={W};
% define time steps
dt1=linspace(0,5,20)*day;dt1=diff(dt1);% time scale of drawdown
dt2=linspace(0.5,20,20)*day;dt2=diff(dt2)
dt=[repmat(1e-3*day,1,5),dt1,dt2]*tfac;

n_steps=numel(dt);
step=struct('control',ones(n_steps,1),'val',dt);
schedule=struct('control',struct('W',W_c),'step',step);


system = initADISystem({'Oil','Water'}, G, rock, fluid, 'cpr', false);
system.nonlinear.cprRelTol = 1e-8;
% setting equations explicitely to be temperature
if(false)
    system.getEquations =@ eqsfiOWT;
    system.stepFunction =@(state0, state, meta, dt, W, G,system)...
        stepOWT(state0, state, meta, dt, G, W, system,fluid);
else
    system.getEquations =@ eqsfiOWFullT;
    system.stepFunction =@(state0, state, meta, dt, W, G,system)...
        stepOWFullT(state0, state, meta, dt, G, W, system,fluid); 
end

% computeing rock conductivity
fake_rock.perm=4.0*ones(G.cells.num,1);
T = computeTrans(G,fake_rock);
Trans=1./accumarray(G.cells.faces(:,1),1./T,[G.faces.num,1]);
internal=all(G.faces.neighbors>0,2);
system.s.T_r=Trans(internal);


% initialize
clear state;
state.pressure = ones(G.cells.num,1)*p_res;
state.s        = ones(G.cells.num,1)*[0.0 1.0 ];
state.T        = ones(G.cells.num,1)*T_res;

[wellSols, states] = runScheduleADI(state, G, rock, system, schedule);
% for runing with dynamic timesteps
%schedule.W=W_c;
%[wellSols, states] = runMrstADI(state, G, system, schedule,'dt_min', 0.01*day,'force_step',false,'targetIts',6);
 
s = cellfun(@(x)x(1).cqs, wellSols, 'UniformOutput', false);
ss = vertcat(s{:})*day;
t  = [cumsum(schedule.step.val)/day ];
figure(1),clf
plot(t,ss)
%
figure(2)
xx=1:G.cells.num;
for k=1:numel(states)
    subplot(3,1,1),cla
    plot(states{k}.pressure/barsa)
    subplot(3,1,2),cla
    plot(states{k}.T)
     subplot(3,1,3),cla
    plot(xx,states{k}.s(:,1),'b',xx,states{k}.s(:,2),'r')
    pause(0.1)
end
%
figure(33)
Tw=[];Pw=[];
for k=2:numel(states)
    Tw=[Tw;states{k}.T(W(1).cells)];
    Pw=[Pw;states{k}.pressure(W(1).cells)]   
end
%%
figure(3)
subplot(2,1,1)
plot(t,(Tw-T_res),'r',t,bsxfun(@rdivide,ss,sum(ss,2))')
subplot(2,1,2)
plot(t,sum(abs(ss),2));%t,bsxfun(@rdivide,ss,sum(ss,2))')

% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
