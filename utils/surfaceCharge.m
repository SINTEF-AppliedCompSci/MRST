function [state, model] = surfaceCharge(model, state)

    T = model.getProp(state, 'temperature');

    An  = 6.0221413*10^23;       	% avagadros number [#/mol]
    F   = 9.64853399e4;             % Faraday's Constant [C/mol]
    R   = 8.3144621;             	% Gas Constant [J/(K mol)]
    e_o = 8.854187817620e-12;       % permitivity of free space [C/Vm]
    e_w = 87.740 - 0.4008.*(T-273.15) + 9.398e-4.*(T-273.15).^2 - 1.410e-6.*(T-273.15).^3;% Dielectric constant of water
    A   = 1.82e6*(e_w.*T).^(-3/2);
    
    nC = numel(model.componentNames);
    components = cell(1, nC);
    [components{:}] = model.chemicalInputModel.getProps(state,model.componentNames{:}); 
    
    if ~isempty(model.surfInfo)
        
        call = 0;
        for i = 1 : numel(model.surfInfo.master)

            if strcmpi(model.surfInfo.scm{i},'langmuir')
                call = call + 1;
                continue
            end
            
            % grab the correct info
            S = model.surfInfo.s{i};
            a = model.surfInfo.a{i};
            C = model.surfInfo.c{i-call};

            % surface funcitonal group name
            surfName = model.surfInfo.master{i}; 

            % number of species associated with surface
            nSp = numel(model.surfInfo.species{i});
            SpNames = model.surfInfo.species{i};
            charge = model.surfInfo.charge{i};


            switch model.surfInfo.scm{i}
                case 'tlm'

                    % calculate surface and IHP charge
                    sig_0 = 0;
                    sig_1 = 0;
                    sig_2 = 0;
                    for j = 1 : nSp
                        SpInd = strcmpi(SpNames{j}, model.componentNames);
                        sig_0 = sig_0 + charge{j}(1).*components{SpInd};
                        sig_1 = sig_1 + charge{j}(2).*components{SpInd};
                        sig_2 = sig_2 + charge{j}(3).*components{SpInd};
                    end
                    sig_0 = (F./(S.*a)).*sig_0;
                    state = model.setProp(state, [surfName '_sig_0'], sig_0);
                    
                    sig_1 = (F./(S.*a)).*sig_1;
                    state = model.setProp(state, [surfName '_sig_1'], sig_1);

                                        
                    sig_2 = (F./(S.*a)).*sig_2;
                    
                    sig_2d = -sig_1 - sig_0 - sig_2;

                    state = model.setProp(state, [surfName '_sig_2'], sig_2 + sig_2d);
                                        


                case {'ccm','langmuir'}

                    % calculate surface charge
                    sig = 0;
                    for j = 1 : nSp
                        SpInd = strcmpi(model.surfInfo.species{i}{j}, model.componentNames);
                        sig = sig + model.surfInfo.charge{i}{j}.*components{SpInd};
                    end
                    sig = (F./(S.*a)).*sig;

                    state = model.setProp(state, [surfName '_sig'], sig);

                case 'ie'
                    warning('Ion exchange surfaces inherently have no charge, and therefore the surface charge is 0. Ignoring.')

            end
        end
    end
        


end
