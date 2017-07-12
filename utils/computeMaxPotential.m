function [names, mins, maxs] = computeMaxPotential(model, state)
    
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
    
    nC = numel(model.CompNames);
    
    ChargeVector = model.ChargeVector;
    
    ChargeVector = [ChargeVector, zeros(1, nC-numel(ChargeVector))];
    
    comps = cell(1, nC);
    [comps{:}] = model.getProps(state,model.CompNames{:}); 
    
    ionDum = 0;
    for i = 1 : nC
        ionDum = ionDum + (ChargeVector(1,i).^2.*comps{i}).*litre/mol;
    end
    ion = cell(1,model.nC);
    [ion{:}] = deal((1/2)*ionDum);
    
    names ={};
    maxs = {};
    mins = {};
   
    myarcsinh = @(x) log(x + (x.^2 + 1).^(1/2));
%     myarcsinh = @(x) x;
    
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

            switch model.surfInfo.scm{i}
                case 'tlm'

                    
                    sig_2 = (F./(S.*a)).*model.getProp(state, surfName)*litre/mol;
                    sig_0 = -sig_2;
                    
                    P2 = R.*T./F.*myarcsinh(-sig_2./(8*10^3*R.*T.*ion{end}.*e_o.*e_w).^(0.5)).*2;

                    P1 = -(sig_2)./C(:,2) + P2;
                    
                    P0 = sig_0./C(:,1) + P1;
                                        
                    names{end+1} = [surfName '_ePsi_0'];
                    names{end+1} = [surfName '_ePsi_1'];
                    names{end+1} = [surfName '_ePsi_2'];
                    
                    maxs{end+1} = exp(F.*P0./(R.*T));
                    maxs{end+1} = exp(F.*P1./(R.*T));
                    maxs{end+1} = exp(F.*P2./(R.*T));
                    
                    mins{end+1} = exp(F.*-P0./(R.*T));
                    mins{end+1} = exp(F.*-P1./(R.*T));
                    mins{end+1} = exp(F.*-P2./(R.*T));
                    


                case 'ccm'

                    % calculate surface charge
                    max_sig = (F./(S.*a)).*model.getProp(state, surfName)*litre/mol;
                    min_sig = -max_sig;
                    
                    % explicitly calculate what the potential should be
                    max_P0 = max_sig./C(:,1);
                    min_P0 = min_sig./C(:,1);


                    names{end+1} = ['log' surfName '_ePsi'];
                    
                    maxs{end+1} = exp(F.*max_P0./(R*T));
                 
                    mins{end+1} = exp(F.*min_P0./(R*T));

                    
                case 'langmuir'

            end
        end
    end
        

    
end

