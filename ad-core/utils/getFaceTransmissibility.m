function T = getFaceTransmissibility(G, rock, deck, varargin)
    if nargin == 2 || isempty(deck)
        deck = struct('GRID', struct());
    end
    
    if mod(numel(varargin), 2) == 1
        varargin = {deck, varargin{:}};
        deck = struct('GRID', struct());
    end
    T = computeTrans(G, rock, varargin{:});
    
    % Multiply inn transmissibility multipliers
    m = computeTranMult(G, deck.GRID);
    if ~isempty(m)
        T = m.*T;
    end
    % Reduce half-face trans for face transmissibility
    cf = G.cells.faces(:,1);
    nf = G.faces.num;
    T  = 1 ./ accumarray(cf, 1./T, [nf, 1]);
    
    % Treat fault multipliers
    flt = processFaults(G, deck.GRID);
    if ~isempty(flt)
        for i = 1:numel(flt)
            ff = flt(i).faces;
            % Max of zero is to avoid negative values or NaNs
            T(ff) = T(ff).*max(flt(i).mult, 0);
        end
    end
end