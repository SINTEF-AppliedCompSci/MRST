function [state, model] = surfacePotential(model, state)

    T = model.getProp(state, 'temperature');

    F   = 9.64853399e4;             % Faraday's Constant [C/mol]
    R   = 8.3144621;             	% Gas Constant [J/(K mol)]
    
    if ~isempty(model.surfInfo)
        
        for i = 1 : numel(model.surfInfo.master)
            
            % surface funcitonal group name
            surfName = model.surfInfo.master{i}; 
            
            if ismember(model.surfInfo.scm{i},{'langmuir','ie'})
                continue
            end

            switch model.surfInfo.scm{i}
                case 'tlm'

                    preP0 = model.getProp(state, ['log' surfName '_ePsi_0']); 
                    state = model.setProp(state, [surfName '_Psi_0'], R.*T./F.*preP0);
                    
                    preP1 = model.getProp(state, ['log' surfName '_ePsi_1']); 
                    state = model.setProp(state, [surfName '_Psi_1'], R.*T./F.*preP1);
                    
                    preP2 = model.getProp(state, ['log' surfName '_ePsi_2']); 
                    state = model.setProp(state, [surfName '_Psi_2'], R.*T./F.*preP2);
                                       
                case 'ccm'

                    preP0 = model.getProp(state, ['log' surfName '_ePsi']); 
                    state = model.setProp(state, [surfName '_Psi'], R.*T./F.*preP0);
                    
            end
        end
    end
        


end
