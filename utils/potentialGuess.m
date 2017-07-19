function [state] = potentialGuess(model, state)
    
%             comps = cellfun(@(x) x*litre/mol, comps,'UniformOutput', false);
%             masterComps = cellfun(@(x) x*litre/mol, masterComps,'UniformOutput', false);

    try 
        T = model.getProp(state, 'temperature');
    catch
        T = 298;
    end
    
    An  = 6.0221413*10^23;       	% avagadros number [#/mol]
    F   = 9.64853399e4;             % Faraday's Constant [C/mol]
    R   = 8.3144621;             	% Gas Constant [J/(K mol)]
    e_o = 8.854187817620e-12;       % permitivity of free space [C/Vm]
    e_w = 87.740 - 0.4008.*(T-273.15) + 9.398e-4.*(T-273.15).^2 - 1.410e-6.*(T-273.15).^3;% Dielectric constant of water
    A   = 1.82e6*(e_w.*T).^(-3/2);
    
    nC = numel(model.CompNames);
    comps = cell(1, nC);
    [comps{:}] = model.chemicalInputModel.getProps(state,model.CompNames{:}); 
    
    ionDum = 0;
    for i = 1 : nC
        ionDum = ionDum + (model.ChargeVector(1,i).^2.*comps{i}).*litre/mol;
    end
    ion = cell(1,model.nC);
    [ion{:}] = deal((1/2)*abs(ionDum));
    
    if ~isempty(model.surfInfo)
        
        for i = 1 : numel(model.surfaces.groupNames)
            
            surfName = model.surfaces.groupNames{i};
            
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
                    
                        % diffuse layer charge
                        mysinh = @(x) exp(x)./2 - exp(-x)./2;
                        myarcsinh = @(x) log(x + (x.^2 + 1).^(1/2));

                        sig_2d = -sig_1 - sig_0 - sig_2;

                        P2 = myarcsinh(-sig_2d./(8*10^3*R.*T.*ion{end}.*e_o.*e_w).^(0.5)).*2;
                        state = model.chemicalInputModel.setProp(state, ['log' surfName '_ePsi_2'], P2);

                        P1 = F./(R.*T).*(-(sig_2+ sig_2d)./C(:,2) + R.*T.*P2./F);
                        state = model.chemicalInputModel.setProp(state, ['log' surfName '_ePsi_1'], P1);

                        P0 = F./(R.*T).*(sig_0./C(:,1) + R.*T.*P1./F);
                        state = model.chemicalInputModel.setProp(state, ['log' surfName '_ePsi_0'], P0);
                        
                case 'ccm'
                        % explicitly calculate what the potential should be
                        Po = sig_0./C(:,1);
                        logPo = F.*Po./(R.*T);

                        state = model.chemicalInputModel.setProp(state, ['log' surfName '_ePsi'], logPo);
            end
            
        end
    end
        

    
end

