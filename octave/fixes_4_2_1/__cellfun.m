function d = cellfun(fn, data, varargin)
    if nargin > 2
        asCell = varargin{2};
    else
        asCell = false;
    end
    
    if asCell
        d = cell(size(data));
        parfor i = 1:numel(data)
            d{i} = fn(data{i});
        end
    else
        d = zeros(size(data));
        parfor i = 1:numel(data)
            d(i) = fn(data{i});
        end

    end
end