function W = setWellSign(W)
    for i = 1:numel(W)
        if isempty(W(i).sign)
            warning('mrst:badWellSign', ...
                ['Well ', W(i).name, ' has empty sign, guessing type']);
            W(i).sign = determineSign(W(i));
        end
    end
end

function s = determineSign(w)
    if strcmpi(w.type, 'pressure') || strcmpi(w.type, 'bhp')
        s = 0;
    else
        s = 1 - 2*(w.val < 0);
    end
end