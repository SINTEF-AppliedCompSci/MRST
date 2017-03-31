function fluid = addSolventProperties(fluid, varargin)

opt = struct('satMisc', @(sG, sS) sS./(sS + sG));

opt = merge_options(opt, varargin{:});

%% Small solvent saturation: Immiscible. Solvent as foruth component

krGT = @(sG, sS) fluid.krG(sG, sS);

krS = @(sG, sS) krGT(sG, sS).*(sS./(sS + sG));
krG = @(sG, sS) krGT(sG, sS).*(sG./(sS + sG));

%% Solvent displacing oil: Miscible

krO = @(sO, sG, sS) sO./(sO + sG + sS).*krOW(sO + sG + sS);
krGT = @(sO, sG, sS) (sS + sG)./(sO + sG + sS).*krOW(sO + sG + sS);

krS = @(sG, sS) krGT(sG, sS).*(sS./(sS + sG));
krG = @(sG, sS) krGT(sG, sS).*(sG./(sS + sG));