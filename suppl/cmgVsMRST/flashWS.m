function [wsOldcorr,qProdTimes] = flashWS(wsOld, varargin)
%FLASHWS Summary of this function goes here
%   Detailed explanation goes here

    opt = struct('casename', 'bakken_light_5comps',...
                 'eosname', 'prcorr',...
                 'p_sc', 101325,...
                 'T_sc', 288.706,...
                 'recoveryMethod', 'CGEOR');

    opt = merge_options(opt, varargin{:});
    
    [fluid, ~] = getShaleCompFluidCase(opt.casename);
    EOSModel = EquationOfStateModel([], fluid, opt.eosname);

    wsOldcorr = wsOld;
    qProdLiq = zeros(numel(wsOld),1);
    qProdGas = zeros(numel(wsOld),1);
    qProdTimes = zeros(numel(wsOld),1);
    for ts=1:numel(wsOld)
        switch(lower(opt.recoveryMethod))
            case {'sdeor', 'primary'}
                qiProd = wsOld{ts}(1).components; % Producer comes first in SDEOR or primary production
            case 'cgeor'
                qiProd = wsOld{ts}(2).components; % Producer comes second in CGEOR
            otherwise
                error('Unknown recovery method')
        end
        if(isempty(qiProd))
            qiProd = zeros(1,numel(fluid.molarMass)); % in kg/s  (for each component)
        else
            qiProd = sum(qiProd,1);          % in kg/s  (for each component)
        end
        qProdTimes(ts) = sum(qiProd);        % in kg/s for total fluid mixture
        qiProd = qiProd ./ fluid.molarMass;  % convert to mol/s (for each component)
        qProd = sum(qiProd);                 % in mol/s for entire fluid mixture

        zi = qiProd ./ qProd;                % overall mole-fraction (matched info.initial b4 gas breaks through
        [Lsc, xsc, ysc, ~, ~, rhoO_Ssc, rhoG_Ssc] = standaloneFlash(opt.p_sc, opt.T_sc, zi, EOSModel);
        qiL = Lsc*qProd*xsc;                 % n^i_L = x_i*L*n bcos n^i_L = x_i*n_L, and n_L = L*n; qiL=n^i_L/time
        qiV = (1-Lsc)*qProd*ysc;             % similar to qiL
        qL = sum(qiL);                       % liquid rate in mol/s
        qV = sum(qiV);                       % gas rate in mol/s
        cL = rhoO_Ssc/sum(xsc.*fluid.molarMass); %liquid molar density in mol/m^3  %Same as p_sc/(Z_Lsc*R*T_sc);
        cV = rhoG_Ssc/sum(ysc.*fluid.molarMass); %gas molar density in mol/m^3     %Same as p_sc/(Z_Vsc*R*T_sc);
        qProdLiq(ts) = qL/cL;                % in m^3/s
        qProdGas(ts) = qV/cV;                % in m^3/s 
        
        
        switch(lower(opt.recoveryMethod))
            case {'sdeor', 'primary'}
                wsOldcorr{ts}(1).qOs = qProdLiq(ts);    % override qOs using a correction to old approach
                wsOldcorr{ts}(1).qGs = qProdGas(ts);    % override qOs using a correction to old approach
            case 'cgeor'
                wsOldcorr{ts}(2).qOs = qProdLiq(ts);    % override qOs using a correction to old approach
                wsOldcorr{ts}(2).qGs = qProdGas(ts);    % override qOs using a correction to old approach
            otherwise
                error('Unknown recovery method')
        end
        
    end


end

