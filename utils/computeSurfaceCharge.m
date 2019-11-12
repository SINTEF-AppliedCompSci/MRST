function [state, model] = computeSurfaceCharge(model, state)

    chemsys = model.chemicalSystem;

    T = model.getProp(state, 'temperature');

    An  = 6.0221413*10^23;       	% avagadros number [#/mol]
    F   = 9.64853399e4;             % Faraday's Constant [C/mol]
    R   = 8.3144621;             	% Gas Constant [J/(K mol)]
    e_o = 8.854187817620e-12;       % permitivity of free space [C/Vm]
    e_w = 87.740 - 0.4008.*(T-273.15) + 9.398e-4.*(T-273.15).^2 - 1.410e-6.*(T-273.15).^3;% Dielectric constant of water
    A   = 1.82e6*(e_w.*T).^(-3/2);
    
    nC = chemsys.nC;
    species = cell(1, nC);
    [species{:}] = model.getProps(state, chemsys.speciesNames{:}); 
    
    if ~isempty(chemsys.surfInfo)
        
        % Create surfaceCharge field for variable state, if not existing, and assign default values
        if ~isfield(state, 'surfaceCharges')
            elementsforsize = model.getProp(state, 'elements');
            ncells = size(elementsforsize, 1);
            clear elementsforsize
            ncomp = numel(chemsys.surfaceChargeNames);
            state.surfaceCharges = ones(ncells, ncomp);
        end
        
        for i = 1 : numel(chemsys.surfInfo.master)

            % grab the correct info
            S = chemsys.surfInfo.s{i};
            a = chemsys.surfInfo.a{i};

            % surface funcitonal group name
            surfName = chemsys.surfInfo.master{i}; 

            % number of species associated with surface
            nSp = numel(chemsys.surfInfo.species{i});
            SpNames = chemsys.surfInfo.species{i};
            charge = chemsys.surfInfo.charge{i};


            switch chemsys.surfInfo.scm{i}
                case 'tlm'

                    % calculate surface and IHP charge
                    sig_0 = 0;
                    sig_1 = 0;
                    sig_2 = 0;
                    for j = 1 : nSp
                        SpInd = strcmpi(SpNames{j}, chemsys.speciesNames);
                        sig_0 = sig_0 + charge{j}(1).*species{SpInd};
                        sig_1 = sig_1 + charge{j}(2).*species{SpInd};
                        sig_2 = sig_2 + charge{j}(3).*species{SpInd};
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
                        SpInd = strcmpi(chemsys.surfInfo.species{i}{j}, chemsys.speciesNames);
                        sig = sig + chemsys.surfInfo.charge{i}{j}.*species{SpInd};
                    end
                    sig = (F./(S.*a)).*sig;

                    state = model.setProp(state, [surfName '_sig'], sig);

            end
        end
    end
        


end
