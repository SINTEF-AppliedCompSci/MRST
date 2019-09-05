function [eqs, names, types] = equationsChemicalLog(model, state, logComponents, logMasterComponents, combinationComponents, ...
                                                       logPartialPressures, logSaturationIndicies,logSurfaceAcitivityCoefficients)


    T = model.getProp(state, 'temperature');
    
    
    An  = 6.0221413*10^23;       	% avagadros number [#/mol]
    F   = 9.64853399e4;             % Faraday's Constant [C/mol]
    R   = 8.3144621;             	% Gas Constant [J/(K mol)]
    e_o = 8.854187817620e-12;       % permitivity of free space [C/Vm]
    e_w = 87.740 - 0.4008.*(T-273.15) + 9.398e-4.*(T-273.15).^2 - 1.410e-6*(T-273.15).^3;% Dielectric constant of water
    A   = 1.82e6*(e_w.*T).^(-3/2);

    CM = model.compositionMatrix;
    RM = model.reactionMatrix;
    GM =  model.gasReactionMatrix;
    SM = model.solidReactionMatrix;
    SPM = model.surfacePotentialMatrix;
    
    
    logSurfAct = logSurfaceAcitivityCoefficients;
    
    components = cellfun(@(x) exp(x), logComponents, 'UniformOutput', false);
    masterComponents = cellfun(@(x) exp(x), logMasterComponents,'UniformOutput', false);

    logK = model.logReactionConstants;

    eqs   = cell(1, model.nR + model.nMC);
    names = cell(1, model.nR + model.nMC);
    types = cell(1, model.nR + model.nMC);

    %% calculate ionic strength
    ionDum = 0;
    nP = model.nP;

    CV = model.chargeVector;
    eInd = strcmpi('e-', model.speciesNames);
    CV(1,eInd) = 0;
    
    for i = 1 : model.nC
        ionDum = ionDum + (CV(1,i).^2.*components{i}).*litre/mol;
    end
    ion = cell(1,model.nC);
    [ion{:}] = deal((1/2)*abs(ionDum));
    
    %% calculate acitivity coefficient by davies equation
    pg = cell(1,model.nC);
    for i = 1 : model.nC
        pg{i} = log(10).*-A.*CV(1,i).^2 .* (ion{i}.^(1/2)./(1 + ion{i}.^(1/2)) - 0.3.*ion{i});
        if CV(1,i) == 0
            pg{i} = ion{i}*0.1;
        end
    end
    
    %% mol fractions
    surfMat = repmat(model.surfMaster, 1, model.nC).*CM;
    surfTest = logical(sum(surfMat));
    
    moleFraction = components;
    
    for i = 1 : model.nC
        if surfTest(i)
            surfDen = 0;
            surfNum = 0;
            for j = 1 : model.nMC
                surfNum = surfNum + CM(j,i).*model.surfMaster(j);
                surfDen = surfDen + double(logical(CM(j,i).*model.surfMaster(j)))*masterComponents{j};
            end
            moleFraction{i} = (surfNum./surfDen).*components{i};
        end

    end
    logMoleFraction = cellfun(@(x) log(x), moleFraction, 'UniformOutput', false);
    
    
    %% reaction matrix
    for i = 1 : model.nR  
        
        eqs{i} = -logK{i}(:);

        % component contribution
        for k = 1 : model.nC
            eqs{i} = eqs{i} + RM(i, k).*(pg{k} + logMoleFraction{k});
        end
        
        % potential contribution
        for k = 1 : model.nP
            eqs{i} = eqs{i} + SPM(i, k).*logSurfAct{k};
        end
        
        % gas reactions
        for k = 1 : model.nG
            eqs{i} = eqs{i} + GM(i,k).*logPartialPressures{k};
        end
        
        % solid reactions
        for k = 1 : model.nS
            eqs{i} = eqs{i} + SM(i,k).*logSaturationIndicies{k};
        end
        
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
            masssum = masssum + CM(i,k).*components{k};
        end
     
        
        eqs{j} = log(masssum) - logMasterComponents{i};

        names{j} = ['Conservation of ', model.elementNames{i}] ;
    end
    
    
    %% surface potentials
    if ~isempty(model.surfInfo)
        call = 0;
        for i = 1 : numel(model.surfaces.groupNames)

            groupNames = model.surfaces.groupNames{i};

            if ismember(model.surfaces.scm{i},{'langmuir','ie'})
                call = call + 1;
                continue
            end

            funcNames = model.surfaces.masterNames{i};

            sig_0 = 0;
            sig_1 = 0;
            sig_2 = 0;

            for j = 1 : numel(funcNames)
                mInd = strcmpi(funcNames{j}, model.surfInfo.master);

                % grab the correct info
                S = model.surfInfo.s{mInd};
                a = model.surfInfo.a{mInd};
                C = model.surfaces.c{i-call};

                % number of species associated with surface
                nSp = numel(model.surfInfo.species{mInd});
                SpNames = model.surfInfo.species{mInd};
                charge = model.surfInfo.charge{mInd};

                switch model.surfaces.scm{i}
                    case 'tlm'
                        % calculate surface charges
                        for k = 1 : nSp
                            SpInd = strcmpi(SpNames{k}, model.speciesNames);
                            sig_0 = sig_0 + (F./(S.*a)).*charge{k}(1).*components{SpInd};
                            sig_1 = sig_1 + (F./(S.*a)).*charge{k}(2).*components{SpInd};
                            sig_2 = sig_2 + (F./(S.*a)).*charge{k}(3).*components{SpInd};
                        end

                    case 'ccm'
                        % calculate surface charge
                        for k = 1 : nSp
                            SpInd = strcmpi(SpNames{k}, model.speciesNames);
                            sig_0 = sig_0 + (F./(S.*a)).*charge{k}(1).*components{SpInd};
                        end
                end
            end

            switch model.surfaces.scm{i}
                case 'tlm'
                    mysinh = @(x) exp(x)./2 - exp(-x)./2;
                    
                    P2Ind = strcmpi([groupNames '_ePsi_2'], model.surfaceActivityCoefficientNames);
                    P1Ind = strcmpi([groupNames '_ePsi_1'], model.surfaceActivityCoefficientNames);
                    P0Ind = strcmpi([groupNames '_ePsi_0'], model.surfaceActivityCoefficientNames);

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
                    Pind = cellfun(@(x) ~isempty(x), regexpi(model.surfaceActivityCoefficientNames, groupNames));
                    eqs{end+1} = -sig_0 + (R*T/F).*logSurfAct{Pind}.*C(:,1);
                    names{end+1} = ['-s + Psi*C ,' groupNames];
                    types{end+1} = [];
            end

        end

    end

    [types{:}] = deal('cell');
    
    %% combination matrix
    for i = 1 : model.nLC
        combSum = 0;
        for k = 1 : model.nC
            combSum = combSum + model.combinationMatrix(i,k).*components{k};
        end
        if any(combSum<=0) || any(combinationComponents{i}<=0)
            eqs{end+1} = (combSum - combinationComponents{i});
        else
            eqs{end+1} = (log(combSum) - log(combinationComponents{i}));
        end
        
        names{end + 1} = [model.combinationNames{i}] ;
        types{end + 1} = 'cell';
    end

end
