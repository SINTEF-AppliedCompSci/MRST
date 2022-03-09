function [isL, isV, isLV] = getPhaseFlagGeothermal(p, h, hL, hV, pcrit)
    isL  = value(h) < hL;% & ~isSC; % Liquid
    isV  = value(h) > hV;% & ~isSC; % Vapor
    isSC = value(p) > pcrit;      % Supercritical
    isV  = isV | (isSC & ~isL);            % Label supercritical as vapor
    isLV = ~isL & ~isV;           % Two-phase
end