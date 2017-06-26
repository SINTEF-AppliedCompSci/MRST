function kr = mergeHalfFaceRelperm(model, kr, varargin)
    opt = struct('type', 'cell');
    [opt, extra] = merge_options(opt, varargin{:});
    
    for ph = 1:numel(kr)
        switch lower(opt.type)
            case 'cell'
                n = model.G.cells.num;
                mapping = model.operators.N;
                kr_next = cell(n, 1);
            case 'face'
                n = size(model.operators.N, 1);
                mapping = repmat((1:n)', 1, 2);
                kr_next = cell(n, 1);
            otherwise
                error('Unknown merge strategy');
        end
        kr_old = kr{ph}.reservoir;
        mapping = mapping(:);
        for i = 1:n
            ix = (mapping == i);
            
            kr_next{i} = mergeTables(kr_old{ix});
        end        
        kr{ph}.reservoir = kr_next;
    end
end

function t = mergeTables(varargin)
    S = [];
    kr = [];
    for i = 1:nargin
        S = [S; varargin{i}.S];
        kr = [kr; varargin{i}.kr];
    end
    [S, ix] = sort(S);
    kr = kr(ix);
    
    t.S = S;
    t.kr = kr;
end