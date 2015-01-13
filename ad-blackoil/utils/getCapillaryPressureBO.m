function [pcOW, pcOG] = getCapillaryPressureBO(f, sW, sG)
    [pcOW, pcOG] = deal(0);
    
    % Check for capillary pressure (p_cow)
    if isfield(f, 'pcOW') && ~isempty(sW)
        pcOW  = f.pcOW(sW);
    end
    
    % Check for capillary pressure (p_cog)
    if isfield(f, 'pcOG') && ~isempty(sG) && nargin > 1
        pcOG  = f.pcOG(sG);
    end
end