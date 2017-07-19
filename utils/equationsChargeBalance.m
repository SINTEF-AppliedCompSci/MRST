function [eqs, names, types] = equationsChargeBalance(model, comps, masterComps)
    

    try 
        T = model.getProp(state, 'temperature');
    catch
        T = 298;
    end
    An  = 6.0221413*10^23;       	% avagadros number [#/mol]
    F   = 9.64853399e4;             % Faraday's Constant [C/mol]
    R   = 8.3144621;             	% Gas Constant [J/(K mol)]
    e_o = 8.854187817620e-12;       % permitivity of free space [C/Vm]
    e_w = 87.740 - 0.4008.*(T-273.15) + 9.398e-4.*(T-273.15).^2 - 1.410e-6*(T-273.15).^3;% Dielectric constant of water
    A   = 1.82e6*(e_w.*T).^(-3/2);
    
    RM = model.ReactionMatrix;
    CM = model.CompositionMatrix;

    CVCind = strcmpi(model.CVC, model.MasterCompNames)';
    
    CM = [CM , CVCind];
    RM = [RM , zeros(model.nR,1)];
    
    nC = model.nC + 1;
        
    ChargeVector = model.ChargeVector;
    ChargeVector = [ChargeVector, zeros(1, nC - size(ChargeVector,2))];
    
    logcomps = cellfun(@(x)log(x), comps(1:end-1),'UniformOutput', false);
    
    logK = model.LogReactionConstants;

    eqs   = cell(1, model.nR + model.nMC + 1);
    names = cell(1, model.nR + model.nMC + 1);
    types = cell(1, model.nR + model.nMC + 1);


    
    % calculate ionic strength
    ionDum = 0;
    for i = 1 : nC
        ionDum = ionDum + (ChargeVector(1,i).^2.*comps{i}).*litre/mol;
    end
    ion = cell(1,model.nC);
    [ion{:}] = deal((1/2)*abs(ionDum));
    
    % activity coefficients
    pg = cell(1,model.nC);
    
    for i = 1 : nC-1
        pg{i} = log(10).*-A.*ChargeVector(1,i)'.^2 .* (ion{i}.^(1/2)./(1 + ion{i}.^(1/2)) - 0.3.*ion{i});
    end

    % reaction matrix,
    for i = 1 : model.nR
        eqs{i} = - logK(i);
        for k = 1 : model.nC-1
            eqs{i} = eqs{i} + RM(i, k).*(pg{k} + logcomps{k});
        end
        names{i} = model.rxns{i};
    end
    
    % composition matrix
    for i = 1 : model.nMC
        j = model.nR + i;
        eqs{j} = - masterComps{i};
        for k = 1 : nC
            eqs{j} = eqs{j} + CM(i,k).*comps{k};
        end
        names{j} = ['Conservation of ', model.MasterCompNames{i}] ;
    end
    

    eqs{end} = 0;
    for k = 1 : nC
        eqs{end} = eqs{end} + ChargeVector(1,k).*comps{k};
    end
    names{end} = 'charge balance equation' ;
        
    if ~isempty(model.surfInfo)
        for i = 1 : numel(model.surfaces.groupNames)

            surfName = model.surfaces.groupNames{i};

            if strcmpi(model.surfaces.scm{i},'langmuir')
                call = call + 1;
                continue
            end

            mNames = model.surfaces.speciesNames(i,:);
            mNames = mNames(cellfun(@(x) ~isempty(x), mNames));

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

                    sig_2 = sig_2 + -(8*10^3*R*T.*ion{end}.*e_o*e_w).^(0.5).*mysinh(logcomps{P2Ind}./2);
                    
                    
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

