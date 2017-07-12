function [state] = surfaceCharge(model, state)

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
    
    ChargeVector = model.ChargeVector;
    
    ChargeVector = [ChargeVector, zeros(1, nC-numel(ChargeVector))];
    
    ionDum = 0;
    for i = 1 : nC
        ionDum = ionDum + (ChargeVector(1,i).^2.*comps{i}).*litre/mol;
    end
    ion = cell(1,model.nC);
    [ion{:}] = deal((1/2)*ionDum);
    
    if ~isempty(model.surfInfo)
        
        call = 0;
        for i = 1 : numel(model.surfInfo.master)

            if strcmpi(model.surfInfo.scm{i},'langmuir')
                call = call + 1;
                continue
            end
            
            % grab the correct info
            S = model.surfInfo.s{i}*gram/(meter)^2;
            a = model.surfInfo.a{i}*litre/gram;
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
                        SpInd = strcmpi(SpNames{j}, model.CompNames);
                        sig_0 = sig_0 + charge{j}(1).*comps{SpInd}.*litre/mol;
                        sig_1 = sig_1 + charge{j}(2).*comps{SpInd}.*litre/mol;
                        sig_2 = sig_2 + charge{j}(3).*comps{SpInd}.*litre/mol;
                    end
                    sig_0 = (F./(S.*a)).*sig_0*mol/litre;
                    state = model.setProp(state, [surfName '_sig_0'], sig_0);
                    
                    sig_1 = (F./(S.*a)).*sig_1*mol/litre;
                    state = model.setProp(state, [surfName '_sig_1'], sig_1);

                                        
                    sig_2 = (F./(S.*a)).*sig_2*mol/litre;
                    
                    sig_2d = -sig_1 - sig_0 - sig_2;

                    state = model.setProp(state, [surfName '_sig_2'], sig_2 + sig_2d);
                                        


                case {'ccm','langmuir'}

                    % calculate surface charge
                    sig = 0;
                    for j = 1 : nSp
                        SpInd = strcmpi(model.surfInfo.species{i}{j}, model.CompNames);
                        sig = sig + model.surfInfo.charge{i}{j}.*comps{SpInd}.*litre/mol;
                    end
                    sig = (F./(S.*a)).*sig*mol/litre;

                    state = model.setProp(state, [surfName '_sig'], sig);

                                        
                    % explicitly calculate what the potential should be
                    Po = sig./C(:,1);
                    logPo = F.*Po./(R.*T);

                    state = model.chemicalInputModel.setProp(state, ['log' surfName '_Psi'], logPo);

                case 'ie'
                    error('Ion exchange surfaces inherently have no charge, and therefore the surface charge is 0.')

            end
        end
    end
        


end
