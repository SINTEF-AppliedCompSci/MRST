function [names, mins, maxs] = computeMaxPotential(model, state)
    
    chemsys = model.chemicalSystem;
    
    T = model.getProp(state, 'temperature');
    
    An  = 6.0221413*10^23;       	% avagadros number [#/mol]
    F   = 9.64853399e4;             % Faraday's Constant [C/mol]
    R   = 8.3144621;             	% Gas Constant [J/(K mol)]
    e_o = 8.854187817620e-12;       % permitivity of free space [C/Vm]
    e_w = 87.740 - 0.4008.*(T-273.15) + 9.398e-4.*(T-273.15).^2 - 1.410e-6.*(T-273.15).^3;% Dielectric constant of water
    
    nC = numel(chemsys.speciesNames);
    
    CV = chemsys.chargeVector;
    eInd = strcmpi('e-', chemsys.speciesNames);
    CV(1,eInd) = 0;
    
    comps = cell(1, nC);
    [comps{:}] = model.getProps(state,chemsys.speciesNames{:}); 
    
    ionDum = 0;
    for i = 1 : nC
        ionDum = ionDum + (CV(1,i).^2.*comps{i}).*litre/mol;
    end
    ion = cell(1,chemsys.nC);
    [ion{:}] = deal((1/2)*ionDum);
    
    names ={};
    maxs = {};
    mins = {};
   
    myarcsinh = @(x) log(x + (x.^2 + 1).^(1/2));
    
    if ~isempty(chemsys.surfInfo)
              
        for i = 1 : numel(chemsys.surfaces.groupNames)
            
            groupName = chemsys.surfaces.groupNames{i};
            
            SPNames = chemsys.surfaces.masterNames{i};
            
            sig_0 = 0;
            sig_1 = 0;
            sig_2 = 0;
                        
            for j = 1 : numel(SPNames)
                mInd = strcmpi(SPNames{j}, chemsys.surfInfo.master);
                                        
                % grab the correct info
                S = chemsys.surfInfo.s{mInd};
                a = chemsys.surfInfo.a{mInd};
                C = chemsys.surfaces.c{i};


                switch chemsys.surfaces.scm{i}
                    case 'tlm'

                        % calculate surface charges
                        sig_2 = sig_2 + (F./(S.*a)).*model.getProp(state, SPNames{j});

                        
                    case 'ccm'

                        % calculate surface charges
                        sig_0 = sig_0 + (F./(S.*a)).*model.getProp(state, SPNames{j});
                end
            end
            

            lim = log(realmax)-10;
            
            switch chemsys.surfaces.scm{i}
                case 'tlm'
                    
                    % diffuse layer charge
                    mysinh = @(x) exp(x)./2 - exp(-x)./2;
                    myarcsinh = @(x) log(x + (x.^2 + 1).^(1/2));
                    
                    sig_0 = -sig_2;
                    
                    P2 = R.*T./F.*myarcsinh(-sig_2./(8*10^3*R.*T.*ion{end}.*e_o.*e_w).^(0.5)).*2;

                    P1 = -(sig_2)./C(:,2) + P2;
                    
                    P0 = sig_0./C(:,1) + P1;
                    
                    upLim = exp(lim.*ones(numel(P0),1));
                    loLim = exp(-lim.*ones(numel(P0),1));
            
                    names{end+1} = [groupName '_ePsi_0'];
                    names{end+1} = [groupName '_ePsi_1'];
                    names{end+1} = [groupName '_ePsi_2'];
                    

                    
                    maxs{end+1} = max([exp(F.*-P0./(R.*T)), exp(F.*P0./(R.*T))], [], 2);
                    maxs{end} = min([maxs{end}, upLim],[],2 );
                    
                    mins{end+1} = min([exp(F.*-P0./(R.*T)), exp(F.*P0./(R.*T))], [], 2);
                    mins{end} = max([mins{end}, loLim],[], 2);
                    
                    maxs{end+1} = max([exp(F.*-P1./(R.*T)), exp(F.*P1./(R.*T))], [], 2);
                    maxs{end} = min([maxs{end}, upLim],[],2 );
                    
                    mins{end+1} = min([exp(F.*-P1./(R.*T)), exp(F.*P1./(R.*T))], [], 2);
                    mins{end} = max([mins{end}, loLim],[], 2);
                    
                    maxs{end+1} = max([exp(F.*-P2./(R.*T)), exp(F.*P2./(R.*T))], [], 2);
                    maxs{end} = min([maxs{end}, upLim],[],2 );
                    
                    mins{end+1} = min([exp(F.*-P2./(R.*T)), exp(F.*P2./(R.*T))], [], 2);
                    mins{end} = max([mins{end}, loLim],[], 2); 
                    
                case 'ccm'

                    % calculate surface charge
                    max_sig = sig_0;
                    min_sig = -max_sig;
                    
                    upLim = exp(lim.*ones(numel(sig_0),1));
                    loLim = exp(-lim.*ones(numel(sig_0),1));
                    
                    % explicitly calculate what the potential should be
                    max_P0 = max_sig./C(:,1);
                    min_P0 = min_sig./C(:,1);


                    names{end+1} = [groupName '_ePsi'];
                    
                    maxs{end+1} = max([exp(F.*max_P0./(R.*T)), exp(F.*min_P0./(R.*T))], [], 2);
                    maxs{end} = min([maxs{end}, upLim], [],2 );

                    mins{end+1} = min([exp(F.*max_P0./(R.*T)), exp(F.*min_P0./(R.*T))], [], 2);
                    mins{end} = max([mins{end}, loLim], [], 2); 

            end
            
        end
        
        

    end
        

    
end

