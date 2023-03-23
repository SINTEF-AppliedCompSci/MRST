function [eqs, names, types] = equationsChemicalLog(model, state)

    chemsys = model.chemicalSystem;

    logSpecies  = model.getPropAsCell(state, 'logSpecies');
    logElements = model.getPropAsCell(state, 'logElements');
    if chemsys.nLC > 0
        combinationComponents = model.getPropAsCell(state, 'combinationComponents');
    end
    if chemsys.nG > 0
        logPartialPressures = model.getPropAsCell(state, 'logPartialPressures');
    end
    if chemsys.nS > 0    
        logSaturationIndicies = model.getPropAsCell(state, 'logSaturationIndicies');
    end
    if chemsys.nP > 0    
        logSurfaceActivityCoefficients = model.getPropAsCell(state, 'logSurfaceActivityCoefficients');
        logSurfAct = logSurfaceActivityCoefficients; % shortcut
    end
     
    T = model.getProp(state, 'temperature');
    
    An  = 6.0221413*10^23;    % avagadros number [#/mol]
    F   = 9.64853399e4;       % Faraday's Constant [C/mol]
    R   = 8.3144621;          % Gas Constant [J/(K mol)]
    e_o = 8.854187817620e-12; % permitivity of free space [C/Vm]
    e_w = 87.740 - 0.4008.*(T-273.15) + 9.398e-4.*(T-273.15).^2 - 1.410e-6*(T-273.15).^3;% Dielectric constant of water
    A   = 1.82e6*(e_w.*T).^(-3/2);

    CM  = chemsys.compositionMatrix;
    RM  = chemsys.reactionMatrix;
    GM  = chemsys.gasReactionMatrix;
    SM  = chemsys.solidReactionMatrix;
    SPM = chemsys.surfacePotentialMatrix;
    
    species = cellfun(@(x) exp(x), logSpecies, 'UniformOutput', false);
    elements = cellfun(@(x) exp(x), logElements,'UniformOutput', false);

    logK = chemsys.logReactionConstants;

    eqs   = cell(1, chemsys.nR + chemsys.nMC);
    names = cell(1, chemsys.nR + chemsys.nMC);
    types = cell(1, chemsys.nR + chemsys.nMC);

    %% calculate ionic strength
    ionDum = 0;
    nP = chemsys.nP;

    CV = chemsys.chargeVector;
    eInd = strcmpi('e-', chemsys.speciesNames);
    CV(1,eInd) = 0;
    
    for i = 1 : chemsys.nC
        ionDum = ionDum + (CV(1,i).^2.*species{i}).*litre/mol;
    end
    ion = cell(1, chemsys.nC);
    [ion{:}] = deal((1/2)*abs(ionDum));
    
    %% calculate acitivity coefficient by davies equation
    pg = cell(1,chemsys.nC);
    for i = 1 : chemsys.nC
        pg{i} = log(10).*-A.*CV(1,i).^2 .* (ion{i}.^(1/2)./(1 + ion{i}.^(1/2)) - 0.3.*ion{i});
        if CV(1,i) == 0
            pg{i} = ion{i}*0.1;
        end
    end
    
    %% mol fractions
    surfMat = repmat(chemsys.surfMaster, 1, chemsys.nC).*CM;
    surfTest = logical(sum(surfMat));
    
    moleFraction = species;
    
    for i = 1 : chemsys.nC
        if surfTest(i)
            surfDen = 0;
            surfNum = 0;
            for j = 1 : chemsys.nMC
                surfNum = surfNum + CM(j,i).*chemsys.surfMaster(j);
                surfDen = surfDen + double(logical(CM(j,i).*chemsys.surfMaster(j)))*elements{j};
            end
            moleFraction{i} = (surfNum./surfDen).*species{i};
        end

    end
    logMoleFraction = cellfun(@(x) log(x), moleFraction, 'UniformOutput', false);
    
    
    %% reaction matrix
    for i = 1 : chemsys.nR  
        
        eqs{i} = -logK{i}(:);

        % component contribution
        for k = 1 : chemsys.nC
            eqs{i} = eqs{i} + RM(i, k).*(pg{k} + logMoleFraction{k});
        end
        
        % potential contribution
        for k = 1 : chemsys.nP
            eqs{i} = eqs{i} + SPM(i, k).*logSurfaceActivityCoefficients{k};
        end
        
        % gas reactions
        for k = 1 : chemsys.nG
            eqs{i} = eqs{i} + GM(i,k).*logPartialPressures{k};
        end
        
        % solid reactions
        for k = 1 : chemsys.nS
            eqs{i} = eqs{i} + SM(i,k).*logSaturationIndicies{k};
        end
        
        names{i} = chemsys.rxns{i};
    end

    assert(all(all(CM>=0)), ['this implementation only supports positive ' ...
                        'master components values']);
                    
    %% composition matrix
    for i = 1 : chemsys.nMC
        j = chemsys.nR + i;
        masssum = 0;
        
        % now using units of moles
        for k = 1 : chemsys.nC
            masssum = masssum + CM(i,k).*species{k};
        end
     
        
        eqs{j} = log(masssum) - logElements{i};

        names{j} = ['Conservation of ', chemsys.elementNames{i}] ;
    end
    
    
    %% surface potentials
    if ~isempty(chemsys.surfInfo)
        call = 0;
        for i = 1 : numel(chemsys.surfaces.groupNames)

            groupNames = chemsys.surfaces.groupNames{i};

            if ismember(chemsys.surfaces.scm{i},{'langmuir','ie'})
                call = call + 1;
                continue
            end

            funcNames = chemsys.surfaces.masterNames{i};

            sig_0 = 0;
            sig_1 = 0;
            sig_2 = 0;

            for j = 1 : numel(funcNames)
                mInd = strcmpi(funcNames{j}, chemsys.surfInfo.master);

                % grab the correct info
                S = chemsys.surfInfo.s{mInd};
                a = chemsys.surfInfo.a{mInd};
                C = chemsys.surfaces.c{i-call};

                % number of species associated with surface
                nSp = numel(chemsys.surfInfo.species{mInd});
                SpNames = chemsys.surfInfo.species{mInd};
                charge = chemsys.surfInfo.charge{mInd};

                switch chemsys.surfaces.scm{i}
                    case 'tlm'
                        % calculate surface charges
                        for k = 1 : nSp
                            SpInd = strcmpi(SpNames{k}, chemsys.speciesNames);
                            sig_0 = sig_0 + (F./(S.*a)).*charge{k}(1).*species{SpInd};
                            sig_1 = sig_1 + (F./(S.*a)).*charge{k}(2).*species{SpInd};
                            sig_2 = sig_2 + (F./(S.*a)).*charge{k}(3).*species{SpInd};
                        end

                    case 'ccm'
                        % calculate surface charge
                        for k = 1 : nSp
                            SpInd = strcmpi(SpNames{k}, chemsys.speciesNames);
                            sig_0 = sig_0 + (F./(S.*a)).*charge{k}(1).*species{SpInd};
                        end
                end
            end

            switch chemsys.surfaces.scm{i}
                case 'tlm'
                    mysinh = @(x) exp(x)./2 - exp(-x)./2;
                    
                    P2Ind = strcmpi([groupNames '_ePsi_2'], chemsys.surfaceActivityCoefficientNames);
                    P1Ind = strcmpi([groupNames '_ePsi_1'], chemsys.surfaceActivityCoefficientNames);
                    P0Ind = strcmpi([groupNames '_ePsi_0'], chemsys.surfaceActivityCoefficientNames);

                    sig_2 = sig_2 + -(8.*10.^3.*R.*T.*ion{end}.*e_o.*e_w).^(0.5).*mysinh(logSurfAct{P2Ind}./2);
                    
                    eqs{end+1} = sig_0 + sig_1 + sig_2;
                    names{end+1} = ['charge balance of ' groupNames];
                    types{end+1} = [];
                    
                    eqs{end+1} = -sig_0 + C(:,1).*(R*T)./F.*(logSurfAct{P0Ind} - logSurfAct{P1Ind});
                    names{end+1} = ['-s0 + C1*(P0 - P1), ' groupNames];
                    types{end+1} = [];
                    
                    eqs{end+1} = -sig_2 - C(:,2).*(R*T)./F.*(logSurfAct{P1Ind} - logSurfAct{P2Ind});
                    names{end+1} = ['-s2 - C2*(P1 - P2), ' groupNames];
                    types{end+1} = [];

                case 'ccm'
                    % explicitly calculate what the potential should be
                    Pind = cellfun(@(x) ~isempty(x), regexpi(chemsys.surfaceActivityCoefficientNames, groupNames));
                    eqs{end+1} = -sig_0 + (R*T/F).*logSurfAct{Pind}.*C(:,1);
                    names{end+1} = ['-s + Psi*C ,' groupNames];
                    types{end+1} = [];
            end

        end

    end

    [types{:}] = deal('cell');
    
    %% combination matrix
    for i = 1 : chemsys.nLC
        combSum = 0;
        for k = 1 : chemsys.nC
            combSum = combSum + chemsys.combinationMatrix(i,k).*species{k};
        end
        if any(combSum<=0) || any(combinationComponents{i}<=0)
            eqs{end+1} = (combSum - combinationComponents{i});
        else
            eqs{end+1} = (log(combSum) - log(combinationComponents{i}));
        end
        
        names{end + 1} = [chemsys.combinationNames{i}] ;
        types{end + 1} = 'cell';
    end

end
