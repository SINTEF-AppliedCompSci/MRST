%%Post-processing of old code
clc; clear;
load('EagleFord_SDEOR_Primary_rate_1e5_DOE.mat')

wsOld = ws;
% load('OMOEOR_64NF_EOR_60years.mat')
% wsNew = ws;
%% Calculate RF
tinSecs = cumsum(schedule.step.val);
tinDays = tinSecs./86400;
 
%%Post-processing of old code
casename = 'eagleford';%'bakken_light_5comps';%'eagleford'; %'oil_1'
[fluid, info] = getShaleCompFluidCase(casename);
eosname = 'prcorr';  %'srk','rk','prcorr'
EOSModel = EquationOfStateModel([], fluid, eosname);
 
%Surface Conditions
% R = 8.31446261815324;
p_sc = 101325;%1*atm;%101325; %atmospheric pressure
T_sc = 288.706;%300; %288.706;% 60 Farenheit

wsOldcorr = wsOld;
qProdLiq = zeros(numel(wsOld),1);
qProdGas = zeros(numel(wsOld),1);
for ts=1:numel(wsOld)
    if(isempty(wsOld{ts}(1).components))
        qiProd = zeros(1,numel(fluid.molarMass));       % in kg/s  (for each component)
    else
        qiProd = wsOld{ts}.components; %Prod then inj
%         [~, qiProd] = wsOld{ts}.components; %Inj then producer for HnP mat files
        qiProd = sum(qiProd,1);
    end
    qiProd = qiProd ./ fluid.molarMass;  % convert to mol/s (for each component)
    qProd = sum(qiProd);                 % in mol/s for entire fluid mixture
    
    zi = qiProd ./ qProd;                % overall mole-fraction (matched info.initial b4 gas breaks through
    [Lsc, xsc, ysc, Z_Lsc, Z_Vsc, rhoO_Ssc, rhoG_Ssc] = standaloneFlash(p_sc, T_sc, zi, EOSModel);
    qiL = Lsc*qProd*xsc;                 % n^i_L = x_i*L*n bcos n^i_L = x_i*n_L, and n_L = L*n; qiL=n^i_L/time
    qiV = (1-Lsc)*qProd*ysc;             % similar to qiL
    qL = sum(qiL);                       % liquid rate in mol/s
    qV = sum(qiV);                       % gas rate in mol/s
    cL = rhoO_Ssc/sum(xsc.*fluid.molarMass); %liquid molar density in mol/m^3  %Same as p_sc/(Z_Lsc*R*T_sc);
    cV = rhoG_Ssc/sum(ysc.*fluid.molarMass); %gas molar density in mol/m^3     %Same as p_sc/(Z_Vsc*R*T_sc);
    qProdLiq(ts) = qL/cL;                % in m^3/s
    qProdGas(ts) = qV/cV;                % in m^3/s 
    wsOldcorr{ts}(1).qOs = qProdLiq(ts);    % override qOs using a correction to old approach
    wsOldcorr{ts}(1).qGs = qProdGas(ts);    % override qOs using a correction to old approach
    
%     qProdLiq(ts) = qL/cL*1*day/stb;  %stb/d of liq     %qL/cL in m^3/s
%     qProdGas(ts) = qV/cV*1*day*35.314666721489;  %scf/d of gas   % qV/cV % in m^3/s 
end
%% Extracting cumulative gas production (field units)
q_gs_corrected = zeros(numel(wsOldcorr),1);
for i = 1:numel(wsOldcorr)
    q_gs_corrected(i)=-3.051e+6.*wsOldcorr{i}(1).qGs;
end
Q_gs_corrected = cumtrapz(tinDays,q_gs_corrected)/1e6; %cum gas produced in MMscf/D
%% Plot corrected well solution
% plotWellSols(wsOldcorr,cumsum(schedule.step.val))
% data1(isnan(data1))=0;
% Qo = cumtrapz(tinDays,data1)./1e3;
% data11(isnan(data11))=0;
% Qg = cumtrapz(tinDays,data11)./1e6;
%%
% wsCombo = {wsOld,wsOldcorr,wsNew};
% names = {'wsOld','wsOldflashed','wsNew'};
% plotWellSols(wsCombo, cumsum(schedule.step.val), 'datasetnames', names,'field','qOs','linestyles',{'r','b'})
% 
% s0 = [0.23, 0.70, 0.07];
% Boi = wsOld{1,1}(1).qOr/ws{1,1}(1).qOs; %calculate the initial oil formation volume factor from production rates
% STOIIP = 6.289811*sum((G_matrix.cells.volumes .* G_matrix.rock.poro))*(1-s0(1))/Boi %in STB
% 
% Bwi = wsOld{1,1}(1).qWr/ws{1,1}(1).qWs; %calculate the initial oil formation volume factor from production rates
% STWIIP = 6.289811*sum((G_matrix.cells.volumes .* G_matrix.rock.poro))*s0(1)/Bwi %in STB

% ws{5} = wsOldcorr;
% T{5} = T{2};
% names{5} = 'Overall-Oldcorr';
% plotWellSols(ws, T, 'datasetnames', names)