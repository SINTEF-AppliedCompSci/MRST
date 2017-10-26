function [state] = potentialGuess(model, state)
    
    T = model.getProp(state, 'temperature');

    An  = 6.0221413*10^23;       	% avagadros number [#/mol]
    F   = 9.64853399e4;             % Faraday's Constant [C/mol]
    R   = 8.3144621;             	% Gas Constant [J/(K mol)]
    e_o = 8.854187817620e-12;       % permitivity of free space [C/Vm]
    e_w = 87.740 - 0.4008.*(T-273.15) + 9.398e-4.*(T-273.15).^2 - 1.410e-6.*(T-273.15).^3;% Dielectric constant of water
    A   = 1.82e6*(e_w.*T).^(-3/2);
    
    nC = model.nC;
    components = cell(1, nC);
    [components{:}] = model.chemicalInputModel.getProps(state,model.speciesNames{:}); 
    
    CV = model.chargeVector;
    eInd = strcmpi('e-', model.speciesNames);
    CV(1,eInd) = 0;
    
    ionDum = 0;
    for i = 1 : nC
        ionDum = ionDum + (CV(1,i).^2.*components{i}).*litre/mol;
    end
    ion = cell(1,model.nC);
    [ion{:}] = deal((1/2)*abs(ionDum));
    
    %%
    if ~isempty(model.surfInfo)
        call = 0;
        for i = 1 : numel(model.surfaces.groupNames)
            
            groupName = model.surfaces.groupNames{i};
            
            if ismember(model.surfaces.scm{i},{'langmuir','ie'})
                call = call + 1;
                continue
            end
            
            funcNames = model.surfaces.masterNames{i};
            
            sig_0 = 0;
            sig_1 = 0;
            sig_2 = 0;
                        
            for j = 1 : numel(funcNames)
                mInd = ismember(funcNames{j}, model.surfInfo.master);
                                        
                % grab the correct info
                S = model.surfInfo.s{mInd};
                a = model.surfInfo.a{mInd};
                C = model.surfaces.c{i-call};

                % number of species associated with the surface
                nSp = numel(model.surfInfo.species{mInd});
                SpNames = model.surfInfo.species{mInd};
                charge = model.surfInfo.charge{mInd};

                switch model.surfaces.scm{i}
                    case 'tlm'

                        for k = 1 : nSp
                            SpInd = strcmpi(SpNames{k}, model.speciesNames);
                            sig_0 = sig_0 + (F./(S.*a)).*charge{k}(1).*components{SpInd};
                            sig_1 = sig_1 + (F./(S.*a)).*charge{k}(2).*components{SpInd};
                            sig_2 = sig_2 + (F./(S.*a)).*charge{k}(3).*components{SpInd};
                        end

                    case 'ccm'

                        for k = 1 : nSp
                            SpInd = strcmpi(SpNames{k}, model.speciesNames);
                            sig_0 = sig_0 + (F./(S.*a)).*charge{k}(1).*components{SpInd};
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
                        state = model.setProp(state, ['log' groupName '_ePsi_2'], P2);

                        P1 = F./(R.*T).*(-(sig_2+ sig_2d)./C(:,2) + R.*T.*P2./F);
                        state = model.setProp(state, ['log' groupName '_ePsi_1'], P1);

                        P0 = F./(R.*T).*(sig_0./C(:,1) + R.*T.*P1./F);
                        state = model.setProp(state, ['log' groupName '_ePsi_0'], P0);
                        
                case 'ccm'
                        % explicitly calculate what the potential should be
                        Po = sig_0./C(:,1);
                        logPo = F.*Po./(R.*T);

                        state = model.setProp(state, ['log' groupName '_ePsi'], logPo);
            end
            
        end
    end
        

    
end

