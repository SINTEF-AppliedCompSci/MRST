function g = processAdjointGradients(grad, ws, varargin)
opt = struct('controlIx', [], 'controlNames', {{'bhp', 'rate', 'orat', 'wrat', 'grat', 'lrat'}});
opt = merge_options(opt, varargin{:});
% load data if ws is a ResultHandler
if isa(ws, 'ResultHandler')
    ws = ws(ws.getValidIds);
end

[nw, ns] = deal(numel(ws{1}), numel(ws));
ng = cellfun(@numel, grad);
if any(diff(ng(:))) %numel(grad{1}) ~= numel(grad{end}) 
    % some wells have been switced off, fill in zeros
    for k = 1:numel(grad)
        wellstats = [ws{k}.status];
        if ~all(wellstats)
            tmp = zeros(nw, 1);
            tmp(wellstats) = grad{k};
            grad{k} = tmp;
        end
    end
end
grad = horzcat(grad{:});

controls = cellfun(@(x){x.type}', ws, 'UniformOutput', false);
controls = horzcat(controls{:});

g = struct();
for k = 1:numel(opt.controlNames)
    nm = opt.controlNames{k};
    g.(nm) = strcmp(controls, nm).*grad;
    if ~isempty(opt.controlIx)
        g.(nm) = g.(nm)*sparse((1:ns)', opt.controlIx, 1, ns, opt.controlIx(end));
    end
end

end