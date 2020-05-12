classdef AvgMPFA < NTPFA
% Average AvgMPFA
    
    properties
        
        % T is a cell structure containing two matrices which compute the fluxes. Each
        % matrix maps cell pressure values to internal face flux value. We have
        % one matrix for each side of the face (the first matrix for
        % G.faces.neighbors(: , 1) and the second for G.faces.neighbors(:, 2)).
        
        T 
        
    end

    methods

        function avgtpfa = AvgMPFA(model, varargin)
            
            avgtpfa = avgtpfa@NTPFA(model, varargin{:});
            G = model.G;
            OSflux = avgtpfa.OSflux;
            avgtpfa.T = AvgTransNTPFA(G, OSflux);
            
        end

        function v = gradient(ntpfa, pressure)
            
            grad = cell(1, ntpfa.nph);
            T = ntpfa.T;
            
            for i = 1 : ntpfa.nph
                
                v = -0.5*(T{1}*pressure - T{2}*pressure);
                v = v .* ntpfa.scale;
                grad{i} = v;
                
            end
        end
    end

end
