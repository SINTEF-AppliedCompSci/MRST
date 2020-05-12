classdef AvgNTPFA < NTPFA
% Average NTPFA 
    properties
        % T is a cell with two matrice from cell pressure values to internal face flux value.
        % One matrix for each side of the face.
        T 
    end

    methods

        function avgtpfa = AvgNTPFA(model, varargin)
            
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
