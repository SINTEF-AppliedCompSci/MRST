function dynamic = UpdateDynamicParams(model)
%
% DESCRIPTION: save the simulation results from each schedule row during
%              the simulation
%
% SYNOPSIS:
%   dynamic = UpdateDynamicParams(model)
%
% PARAMETERS:
%   model - struct containing the following fields:
%   - bc: boundary conditions
%   - static: static simulation information (e.g. gauge placement)
%   - dynamic: information saved during the simulation of last schedule
%   rows
%   - experiment: saturation functions used for forward modeling
%   - rock: rock properties like porosity and absolute permeability
%   - grid
%
%
% RETURNS:
%   model - struct containing the following fields:
%   - dynamic: information saved during the new schedule row
%
% ----------------------------------
% (c) 2020-2022
% Siroos Azizmohammadi
% Omidreza Amrollahinasab
% Montanuniversität Leoben, Austria
% Chair of Reservoir Engineering
% https://dpe.ac.at/
% ----------------------------------
%
%%
bc          = model.bc;
static      = model.static;
dynamic     = model.dynamic;
params      = dynamic.params;  
Swi         = model.experiment.rock.Swi.value;
simulation  = model.simulation;
process     = model.experiment.process;
processName = process.name;
processType = lower(process.type);    
bcl         = bc.left;            
if (isfield(simulation,'bCells'))
    total_pv = sum(model.rock.pv(2:end-1));
    pv_array = model.rock.pv(2:end-1);
else
    total_pv = sum(model.rock.pv);
    pv_array = model.rock.pv;
end

for i = 1:length(params.scheduleSteps)
    state  = model.state{i, 1};
    if (isfield(simulation,'bCells'))
        Sw = state.s(2:end-1,1);
    else
        Sw = state.s(:,1);
    end               
    
    % injected fluids
    % water_injection_rate
    WIR = params.qinj(:,1);
    % oil_injection_rate
    OIR = params.qinj(:,2);
    % water_injection_volume
    WIV = params.Qinj(:,1);
    % oil_injection_volume
    OIV = params.Qinj(:,2);            
    
    if(processType == "cent")
        WIR = [WIR; 0]; injection_rate(:,1) = WIR;
        OIR = [OIR; 0]; injection_rate(:,2) = OIR;
        WIV = [WIV; 0]; injection_volume(:,1) = WIV;
        OIV = [OIV; 0]; injection_volume(:,2) = OIV;
    else                       
        WIR = [WIR; bcl.qw_inj]; injection_rate(:,1) = WIR;
        OIR = [OIR; bcl.qo_inj]; injection_rate(:,2) = OIR;        
        WIV = [WIV; WIV(end) + WIR(end) * params.scheduleSteps(i)];            
        OIV = [OIV; OIV(end) + OIR(end) * params.scheduleSteps(i)];
        injection_volume(:,1) = WIV; injection_volume(:,2) = OIV;
    end
    params.qinj = injection_rate; 
    params.Qinj = injection_volume;
    clear injection_rate injection_volume

    PVI = params.PVI;
    params.PVI = [PVI; (WIV(end) + OIV(end)) / total_pv];

    % produced fluids  
    % water_production_rate
    WPR = params.qp_net(:,1);
    % oil_production_rate
    OPR = params.qp_net(:,2); 
    % water_production_volume
    WPV = params.Qp_net(:,1);
    % oil_production_volume
    OPV = params.Qp_net(:,2);

    % water volume in the media
    q = sum(pv_array.* Sw);
    % intial water volume
    qi = total_pv * Swi;

    switch lower(processName)
        case "imbibition"
            % we do not report the water production for imbibition case
            WPV = [WPV; nan]; 
            % net_production_volume
            NPV(:,1) = WPV;
            WPR = [WPR; WPV(end) / params.scheduleSteps(i)]; 
            %net_production_rate
            NPR(:,1) = WPR;
            OPV = [OPV; q-qi]; NPV(:,2) = OPV;
            OPR = [OPR; OPV(end) / params.scheduleSteps(i)]; NPR(:,2) = OPR;
        case "drainage"
            WPV = [WPV; qi-q]; NPV(:,1) = WPV;
            WPR = [WPR; WPV(end) / params.scheduleSteps(i)]; NPR(:,1) = WPR;
            % we do not report oil production for drainage case
            OPV = [OPV; nan]; NPV(:,2) = OPV;
            OPR = [OPR; OPV(end) / params.scheduleSteps(i)]; NPR(:,2) = OPR;
    end 

    % Net production
    params.qp_net = NPR;
    params.Qp_net = NPV;
    
    % Cumulative production
    params.qprod = NPR + [WIR, OIR];
    params.Qprod = NPV  + [WIV, OIV];
    clear NPR NPV
    
    % This gives oil pressure
    pDiff = params.pDiff;
    params.pDiff = [pDiff;state.pressure(static.gauge.left)-...
                          state.pressure(static.gauge.right)];            
    params.gradp = gradient(pDiff);

    % This gives oil pressure in the sample middle (in case of 2
    % pressure at the same time)
    if isfield(model.experiment.observation, "pressure_mid")
        pDiff_mid = params.pDiff_mid;
        params.pDiff_mid = [pDiff_mid ;state.pressure(static.gauge_Pmid.left)-...
                              state.pressure(static.gauge_Pmid.right)];            
        params.gradp_mid = gradient(pDiff_mid);
    end


    SwAvg = params.SwAvg;
    params.SwAvg = [SwAvg;mean(Sw)];
    if min(Sw) < params.Sw_min
        params.Sw_min = min(Sw);
    end
    if max(Sw) > params.Sw_max
        params.Sw_max = max(Sw);
    end

    dynamic.params = params;
    model.dynamic = dynamic;        
end