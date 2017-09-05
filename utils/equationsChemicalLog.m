function [eqs, names, types] = equationsChemicalLog(logporo, logcomps, logmasterComps, logGasComps, logSolidComps, state, model)

    T = model.getProp(state, 'temp');

    if model.nG > 0
        partialPressures = cell(1,model.nG);
        [partialPressures{:}] = deal(model.getProps(state, model.partialPressureNames{:}));
        logPartialPressures = cellfun(@(x) log(x), partialPressures, 'UniformOutput', false);    
    end

    if model.nS > 0
        solidDensities = cell(1,model.nS);
        [solidDensities{:}] = deal(model.getProps(state,  model.solidDensityNames{:}));
    end
    
    
    An  = 6.0221413*10^23;       	% avagadros number [#/mol]
    F   = 9.64853399e4;             % Faraday's Constant [C/mol]
    R   = 8.3144621;             	% Gas Constant [J/(K mol)]
    e_o = 8.854187817620e-12;       % permitivity of free space [C/Vm]
    e_w = 87.740 - 0.4008.*(T-273.15) + 9.398e-4.*(T-273.15).^2 - 1.410e-6*(T-273.15).^3;% Dielectric constant of water
    A   = 1.82e6*(e_w.*T).^(-3/2);

    RM = model.ReactionMatrix;
    CM = model.CompositionMatrix;

    comps = cellfun(@(x) exp(x), logcomps, 'UniformOutput', false);
    gasComps = cellfun(@(x) exp(x), logGasComps, 'UniformOutput', false);
    solidComps = cellfun(@(x) exp(x), logSolidComps, 'UniformOutput', false);
    masterComps = cellfun(@(x) exp(x), logmasterComps, 'UniformOutput', false);
    poro = exp(logporo);
    
    logK = model.LogReactionConstants;

    eqs   = cell(1, model.nR + model.nMC);
    names = cell(1, model.nR + model.nMC);
    types = cell(1, model.nR + model.nMC);

    %% calculate ionic strength
    ionDum = 0;
    nP = sum(cellfun(@(x) ~isempty(x), regexpi(model.CompNames, 'psi')));
    model.ChargeVector = [model.ChargeVector, zeros(1,nP)];
    for i = 1 : model.nC
        ionDum = ionDum + (model.ChargeVector(1,i).^2.*comps{i}).*litre/mol;
    end
    ion = cell(1,model.nC);
    [ion{:}] = deal((1/2)*abs(ionDum));
    
    %% calculate acitivity coefficient by davies equation
    pg = cell(1,model.nC);
    for i = 1 : model.nC
        pg{i} = log(10).*-A.*model.ChargeVector(1,i).^2 .* (ion{i}.^(1/2)./(1 + ion{i}.^(1/2)) - 0.3.*ion{i});
    end
    
    %% active fraction for ion exchange surfaces
    af = cell(1,model.nC);
    [af{:}] = deal(0);
    if ~isempty(model.surfInfo)
        for i = 1 : numel(model.surfInfo.master)
            mcName = model.surfInfo.master{i};
            if strcmpi(model.surfInfo.scm{i}, 'ie')
                Sind = strcmpi(mcName, model.MasterCompNames);
                for j = 1 : numel(model.surfInfo.species{i})
                    p = model.surfInfo.species{i}(j);
                    ind = strcmpi(p, model.CompNames);
                    af{ind} = log(comps{ind}/exp(logmasterComps{Sind}));
                end
            end  
        end
    end
    
    %% reaction matrix
    for i = 1 : model.nR  
        
        eqs{i} = -logK(i);
        
        %%%%%%%%%%%%%%%%%%%%%%%%% put switch here %%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % solid and gas phase contributions have been removed from RM, but
        % will still be enforced as if at local chemical equilibrium. need
        % a vector of alphas to switch the solid phase reactions off
        %
        for k = 1 : model.nC
            eqs{i} = eqs{i} + RM(i, k).*(pg{k} + af{k} + logcomps{k});
        end
        
        for k = 1 : model.nG
            eqs{i} = eqs{i} + model.GasReactionMatrix(i,k).*logPartialPressures{k};
        end
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        names{i} = model.rxns{i};
    end

    assert(all(all(CM>=0)), ['this implementation only supports positive ' ...
                        'master components values']);
                    
    %% composition matrix
    for i = 1 : model.nMC
        j = model.nR + i;
        masssum = 0;
        
        % now using units of moles
        for k = 1 : model.nC
            masssum = masssum + CM(i,k).*comps{k}.*poro;
        end
        
        % solidComps and gasComps has units of volume (m^3)
        % need to implement ideal gas law here for the gasses
        for k = 1 : model.nG
            masssum = masssum + model.GasCompMatrix(i,k).*gasComps{k};
        end
        for k = 1 : model.nS
            masssum = masssum + model.SolidCompMatrix(i,k).*solidComps{k}.*solidDensities{k};
        end
        
        eqs{j} = log(masssum) - (logmasterComps{i} + logporo);

        names{j} = ['Conservation of ', model.MasterCompNames{i}] ;
    end
    
    
    %% conservation of volume
    
    vol = poro;
    for i = 1 : model.nS
        vol = vol + solidComps{i};
    end
    for i = 1 : model.nG
        vol = vol + gasComps{i};
    end    
    
    eqs{end+1} = log(1) - log(vol);
    names{end+1} = 'Conservation of volume';
    types{end+1} = [];
    
    %% surface potentials
    if ~isempty(model.surfInfo)
        for i = 1 : numel(model.surfaces.groupNames)

            surfName = model.surfaces.groupNames{i};

            if strcmpi(model.surfaces.scm{i},'langmuir')
                call = call + 1;
                continue
            end

            mNames = model.surfaces.speciesNames{i};

            sig_0 = 0;
            sig_1 = 0;
            sig_2 = 0;

            for j = 1 : numel(mNames)
                mInd = strcmpi(mNames{j}, model.surfInfo.master);

                % grab the correct info
                S = model.surfInfo.s{mInd}*gram/(meter)^2;
                a = model.surfInfo.a{mInd}*litre/gram;
                C = model.surfaces.c{i};

                % number of species associated with surface
                nSp = numel(model.surfInfo.species{mInd});
                SpNames = model.surfInfo.species{mInd};
                charge = model.surfInfo.charge{mInd};

                switch model.surfaces.scm{i}
                    case 'tlm'

                        % calculate surface charges

                        for k = 1 : nSp
                            SpInd = strcmpi(SpNames{k}, model.CompNames);
                            sig_0 = sig_0 + (F./(S.*a)).*charge{k}(1).*comps{SpInd}.*litre/mol;
                            sig_1 = sig_1 + (F./(S.*a)).*charge{k}(2).*comps{SpInd}.*litre/mol;
                            sig_2 = sig_2 + (F./(S.*a)).*charge{k}(3).*comps{SpInd}.*litre/mol;
                        end

                    case 'ccm'

                        % calculate surface charge
                        for k = 1 : nSp
                            SpInd = strcmpi(SpNames{k}, model.CompNames);
                            sig_0 = sig_0 + (F./(S.*a)).*charge{k}(1).*comps{SpInd}.*litre/mol;
                        end

                end
            end

            switch model.surfaces.scm{i}
                case 'tlm'
                    mysinh = @(x) exp(x)./2 - exp(-x)./2;
                    
                    P2Ind = strcmpi([surfName '_ePsi_2'], model.CompNames);
                    P1Ind = strcmpi([surfName '_ePsi_1'], model.CompNames);
                    P0Ind = strcmpi([surfName '_ePsi_0'], model.CompNames);

                    sig_2 = sig_2 + -(8.*10.^3.*R.*T.*ion{end}.*e_o.*e_w).^(0.5).*mysinh(logcomps{P2Ind}./2);
                    
                    
                    eqs{end+1} = sig_0 + sig_1 + sig_2;
                    names{end+1} = ['charge balance of ' surfName];
                    types{end+1} = [];
                    
                    
                    eqs{end+1} = -sig_0 + C(:,1).*(R*T)./F.*(logcomps{P0Ind} - logcomps{P1Ind});
                    names{end+1} = ['-s0 + C1*(P0 - P1), ' surfName];
                    types{end+1} = [];
                    
                    
                    eqs{end+1} = -sig_2 - C(:,2).*(R*T)./F.*(logcomps{P1Ind} - logcomps{P2Ind});
                    names{end+1} = ['-s2 - C2*(P1 - P2), ' surfName];
                    types{end+1} = [];

                case 'ccm'
                    
                    % explicitly calculate what the potential should be
                    Pind = cellfun(@(x) ~isempty(x), regexpi(model.CompNames, [surfName '_']));
                    eqs{end+1} = -sig_0 + (R*T/F).*logcomps{Pind}.*C(:,1);
                    names{end+1} = ['-s + Psi*C ,' surfName];
                    types{end+1} = [];
            end

        end

    end

    [types{:}] = deal('cell');

end
