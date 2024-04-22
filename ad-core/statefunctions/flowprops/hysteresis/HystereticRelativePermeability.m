classdef HystereticRelativePermeability < BaseRelativePermeability
    properties
        % hysteresisArgs = {};
        sg_min = [];
        sg_max = [];
        sg_tmax = [];
        nreg_imbibition = 1;
        minSat = 0.05;
    end
    
    methods
        function gp = HystereticRelativePermeability(model, varargin)
            gp@BaseRelativePermeability(model, varargin{:});
            nc = model.G.cells.num;
            if isempty(gp.regions)
                reg_defaulted = true;
                if isfield(model.rock, 'regions')
                    r = model.rock.regions;
                    if isfield(r, 'saturation')
                        gp.regions = r.saturation;
                        reg_defaulted = false;
                    end
                end
                if reg_defaulted
                    gp.regions = ones(nc, 1);
                end
            end
            gp = gp.dependsOn({'s', 'sMax'}, 'state');
            
            regi = repmat(2, nc, 1);
            if isfield(model.rock, 'regions')
                r = model.rock.regions;
                if isfield(r, 'imbibition')
                    regi = r.imbibition;
                end
            end
            regi  = unique(regi);
            nregi = numel(regi);
            reg = regi-nregi;
            f = model.fluid;
            % Maximum and minimum gas saturation
            if isfield(f, 'sWcon')
                sgmax = 1 - f.sWcon;
                if numel(sgmax) ~= nregi  % temporary
                    sgmax = ones(1, nregi)*sgmax;
                end
                sgmin = zeros(1, nregi);
                sgtmax = sgmax;
                warning('Min gas saturation (primary drainage) set to 0')
            elseif isfield(f, 'krPts')
                sgtmax   = f.krPts.g(regi,2)';
                sgmax      = f.krPts.g(regi, 3)';
                sgmin      = f.krPts.g(reg, 1)';
            else
                sgmax = ones(1, nregi);
                sgtmax = sgmax;
                sgmin = zeros(1, nregi);
                warning('Max gas saturation (primary drainage) set to 1.')
                warning('Min gas saturation (primary drainage) set to 0')
            end
            assert(numel(sgmax) == nregi);

            % Add to fluid object
            % kriVars = cell(1, nregi);
            for n = 1:nregi
                % kriVars{n} = {sgmax(n), sgmin(n), f.krG{n}, f.krG{n+nregi}, sgtmax(n)};
            end
            gp.sg_max = sgmax;
            gp.sg_min = sgmin;
            gp.sg_tmax = sgtmax;
            gp.nreg_imbibition = nregi;
        end

        % --------------------------------------------------------------------%
        function krg = evaluateGasRelPerm(prop, model, sg, state)
            sgMax = model.getProps(state, 'sgMax');
            sgMax = max(sg, sgMax);
            krg = evaluatePhaseRelativePermeabilityWithHysteresis(prop, model, 'g', sg, sgMax);
        end

        function kr = evaluatePhaseRelativePermeabilityWithHysteresis(prop, model, phase, s, sMax) %LS
            f   = model.fluid;
            fn  = f.(['kr', upper(phase)]);
            if isfield(f, 'krHyst')
                nregi = prop.nreg_imbibition;
                regions = model.rock.regions.imbibition;
                [sample, isAD] = getSampleAD(s, sMax);
                kr = zeros(numel(regions), 1);
                if isAD
                    kr = prop.AutoDiffBackend.convertToAD(kr, sample);
                end
                for reg = 1:numel(unique(regions))
                    idCells = regions == reg + (min(regions)-1);
                    isRegi  = f.krHyst == reg + (min(regions)-1);
                    if any(isRegi)
                        kriVars =  {prop.sg_max(reg), prop.sg_min(reg), f.krG{reg}, f.krG{reg+nregi}, prop.sg_tmax(reg)};
                        fni = @(sg, sgi) scanRelPerm(sg, sgi, kriVars, prop.minSat, f.ehystr);
                        kr(idCells) = evaluateFunctionCellSubsetReg(prop, fni, regions, s(idCells), sMax(idCells));
                    else
                        kr(idCells) = evaluateFunctionCellSubsetReg(prop, fn{reg}, regions, s(idCells));
                    end
                end
                %______________________________________________________________
            else  % no hysteresis
                regions = model.rock.regions.saturation;
                kr      = evaluateFunctionCellSubsetReg(prop, fn, regions, s);
            end
        end
    end
end

function [kri, tol, minSat, sgt] = scanRelPerm(sg, sgi, kriVars, minSat, opts)
    % Compute scanning relative permeability curves.
    % sg  = current gas saturation
    % sgi = maximum gas saturation achieved in the current run.
    % __________________________________________________________
    nci = numelValue(sg);

    % Set tolerance for flow reversal and initialize
    assert(nci == numelValue(sgi))   % same n for both sg and sgi
    tol  = 1e-3;
    kri  = nan(nci, 1);
    AD = getSampleAD(sg, sgi);
    if isa(AD, 'GenericAD')
        kri = double2GenericAD(kri, AD);
    elseif isa(AD, 'ADI')
        kri = double2ADI(kri, AD);
    end

    % Inputs
    hystModel = opts{2};
    sgmx      = kriVars{1};
    sgmn      = kriVars{2};
    krG       = kriVars{3};
    krGib     = kriVars{4};
    sgtmx     = kriVars{5};

    % Index for cells where scanning curves are computed. We set a
    % minimum of 5% gas saturation to activate hysteresis.
    imb  = all([value(sg) + tol < value(sgi), ...
                value(sgi) > minSat+sgmn], 2);

    if any(imb)
        kri(~imb) = krG(sg(~imb));
        sgiv  = sgi(imb);
        sgv   = sg(imb);

        if hystModel == 2           % Killough (1976) nonwetting phase.
            a     = opts{4};        % modif. param. for improved converg.
            A     = 1 + a*(sgmx - sgiv);
            C     = 1/(sgtmx - sgmn) - 1/(sgmx-sgmn);
            sgt   = sgmn + (sgiv-sgmn)./(A + C*(sgiv-sgmn));
            snorm = sgtmx + (sgv-sgt)*(sgmx-sgtmx)./(sgiv-sgt);
            v = krGib(snorm).*krG(sgiv)./krG(sgmx);
            kri(imb) = v;

            isAbove = kri > krG(value(sg)) & imb;
            if any(isAbove)
                kri(isAbove) = krG(sg(isAbove));
                % kri(isAbove) = krG(value(sg(isAbove))) - 1e-5;
            end
        else
            error(['Indicated hysteresis model in item 2 of EHYSTR' ...
                ' keyword in the .DATA input file is not supported.'])
        end

        % Checks
        isBelow = all([sgv < sgt, value(kri(imb)) ~= 0], 2);
        assert(~any(isBelow));
        kriv = value(kri);
        assert(all(~isnan(kriv)));                  % no nan values
        assert(all(kriv >= 0));                     % 0 or positive values
        assert(all(kriv <= krG(value(sgmx))));      % max kr not exceeded
    else
        kri = krG(sg);                      % all cells on drainage curve
    end
end
%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
