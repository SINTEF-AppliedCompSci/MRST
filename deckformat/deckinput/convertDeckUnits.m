function deck = convertDeckUnits(deck, varargin)
%Convert ECLIPSE/FrontSim input deck units to MRST conventions
%
% SYNOPSIS:
%   deck = convertDeckUnits(deck)
%   deck = convertDeckUnits(deck, 'pn1', pv1, ...)
%
% PARAMETERS:
%   deck    - An ECLIPSE/FrontSim input deck as defined by function
%             'readEclipseDeck'.
%
%   'pn'/pv - A list of 'key'/value pairs defining optional parameters.
%             The supported options are:
%               - verbose -- Whether or not to emit informational messages
%                            if a section is skipped.
%                            Logical.  Default value: verbose = mrstVerbose
%               - outputUnit -- Optional output unit. Valid options are 
%                               'METRIC', 'FIELD', 'LAB', 'PVT_M', 'SI'
%                               Default value: 'SI'
%
% RETURNS:
%   deck - An ECLIPSE/FrontSim deck where all quantities affected by
%          varying units of measurements have been converted to MRST's
%          strict, SI only, conventions.
%
% SEE ALSO:
%   `readEclipseDeck`, `mrstVerbose`.

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

   opt = struct('verbose', mrstVerbose, 'outputUnit', 'SI');
   opt = merge_options(opt, varargin{:});

   assert (isstruct(deck) && isfield(deck, 'RUNSPEC'), ...
           'Input does not appear to be a valid deck');

   inputUnit = deckUnitName(deck.RUNSPEC);
   u = unitConversionFactors(inputUnit, opt.outputUnit);
   if ~isfield(deck, 'PCUNIT')
       % On the first conversion, we store the pressure unit. This may be
       % used to interpret PC curves as dimensionless J-functions instead.
       deck.PCUNIT = u.press;
   end
   for sect = reshape(fieldnames(deck), 1, [])
      s = sect{1};

      switch s
         case {'RUNSPEC', 'GRID', 'PROPS', 'SOLUTION', 'SCHEDULE'}
            cvrt     = str2func(['convert', s]);
            deck.(s) = cvrt(deck.(s), u);

         case 'UnhandledKeywords'
            continue   % MRST specific

         otherwise
            dispif(opt.verbose, ...
                   'No converter needed in section ''%s''.\n', s);
      end
   end
end

%--------------------------------------------------------------------------
% Private helpers follow
%--------------------------------------------------------------------------

function unitName = deckUnitName(rspec)
validUnitNames = {'METRIC', 'FIELD', 'LAB', 'PVT_M', 'SI'};
occurance = isfield(rspec, validUnitNames);
if sum(occurance) ~= 1
    error(id('USys:Unknown'), ...
        ['Input unit system must be either METRIC, ', ...
        'FIELD, LAB, PVT_M, or SI.']);
else
    unitName = validUnitNames{occurance};
end
end


%--------------------------------------------------------------------------

function rspec = convertRUNSPEC(rspec, u)
   usys = { 'METRIC', 'FIELD', 'LAB', 'PVT_M' 'SI'};

   rspec = rmfield(rspec, usys(isfield(rspec, usys)));

   rspec.(upper(u.unit_out)) = true;
end

%--------------------------------------------------------------------------

function grid = convertGRID(grid, u)
   if isempty(grid), return; end

   for kw = reshape(fieldnames(grid), 1, [])
      key = kw{1};

      switch key
         case {'PERMX' , 'PERMXY', 'PERMXZ', ...
               'PERMYX', 'PERMY' , 'PERMYZ', ...
               'PERMZX', 'PERMZY', 'PERMZ' , ...
               'PERMXX', 'PERMYY', 'PERMZZ'}
            grid.(key) = convertFrom(grid.(key), u.perm);

         case {'DXV'   , 'DYV'   , 'DZV'   , 'DEPTHZ', ...
               'DX'    , 'DY'    , 'DZ'    , 'TOPS'  , ...
               'COORDX', 'COORDY', 'COORDZ',           ...
               'COORD' , 'ZCORN', 'DEPTH'                     }
            grid.(key) = convertFrom(grid.(key), u.length);

         case 'MAPAXES'
            unt = u.length;

            if isfield(grid, 'MAPUNITS')
               mapstr = 'METRES';
               switch grid.MAPUNITS
                  case 'METRES'
                      unt = meter;
                  case 'FEET'
                      unt = ft;
                  otherwise
                     warning('MapUnits:Unknown', ...
                             'Unknown map units ''%s''. Using METRES', ...
                             grid.MAPUNITS);
               end
               grid.MAPUNITS = mapstr;
            end

            grid.(key) = convertFrom(grid.(key), unt);

         case 'MAPUNITS'
            continue;  % Handled in 'MAPAXES'.

         case {'MINPV', 'MINPVV', 'PORV'}
            grid.(key) = convertFrom(grid.(key), u.liqvol_r);

         case 'NNC'
            grid.(key)(:,7) = convertFrom(grid.(key)(:,7), u.trans);

         case {'TRANX', 'TRANY', 'TRANZ'}
            grid.(key) = convertFrom(grid.(key), u.trans);

         case {'THCONR'}
            grid.(key) = convertFrom(grid.(key), u.rockcond);

         case {'PINCH', 'PINCHREG'}
            i    = [1, 3];
            data = convertFrom([ grid.(key){:,i} ], u.length);

            grid.(key)(:, i) = reshape(num2cell(data), [], numel(i));

            clear i data

         case {'SIGMAV', 'SIGMA'}
              grid.(key) = convertFrom(grid.(key), u.invarea);

         case {'DZMTRXV', 'DZMTRX'}
              grid.(key) = convertFrom(grid.(key), u.length);
              
         case 'JFUNC'
             for i = 2:3
                grid.(key){i} = convertFrom(grid.(key){i}, u.jsurftens);
             end

         case {'cartDims',                          ...  % MRST specific
               'ACTNUM'  , 'NTG', 'PORO', 'MULTPV', ...
               'FLUXNUM' ,                          ...
               'FAULTS'  , 'MULTFLT' ,              ...
               'GDORIENT',                          ...  % Strings only
               'MULTNUM' ,                          ...
               'MULTX'   , 'MULTX_'  , ...
               'MULTY'   , 'MULTY_'  , ...
               'MULTZ'   , 'MULTZ_'  , ...
               'OPERNUM' , 'PINCHNUM','INIT'}
            continue;  % Pure scalars (no unit of measure)

         otherwise
            error(id('GRID:NoConverter'), ...
                  'No known unit conversion in GRID for ''%s''.', key);
      end
   end
end

%--------------------------------------------------------------------------

function props = convertPROPS(props, u)
   for kw = reshape(fieldnames(props), 1, [])
      key = kw{1};

      switch key
         case {'DENSITY', 'SDENSITY'}
            props.(key) = convertFrom(props.(key), u.density);

         case 'RKTRMDIR'
            continue; %dimensionless

         case 'ROCK'
            unt         = [u.press, u.compr, u.compr, u.compr, 1, 1];
            unt         = unt(1 : size(props.(key), 2));
            props.(key) = convertFrom(props.(key), unt);

         case 'ROCKTAB'
            unt         = [u.press, ones(1, size(props.(key){1}, 2)-1)];
            for t = 1:numel(props.(key))
                props.(key){t} = convertFrom(props.(key){t}, unt);
            end

         case 'SURFST'
            unt         = [u.concentr, u.surf_tension];
            for t = 1 : numel(props.(key))
                props.(key){t} = convertFrom(props.(key){t}, unt);
            end

         case 'SURFADS'
            unt         = [u.concentr, 1];
            for t = 1 : numel(props.(key))
                props.(key){t} = convertFrom(props.(key){t}, unt);
            end

         case 'SURFROCK'
            unt         = [1, u.concentr];
            props.(key) = convertFrom(props.(key), unt);

         case 'SURFVISC'
            unt         = [u.concentr, u.viscosity];
            for t = 1 : numel(props.(key))
                props.(key){t} = convertFrom(props.(key){t}, unt);
            end

         case {'SPECHEAT'}
            unt         = [u.temp, repmat(u.massheatcapacity,1,3)];
            for t = 1 : numel(props.(key))
                props.(key){t} = convertFrom(props.(key){t}, unt);
            end

         case {'VISCREF'}
            unt         = [u.press, u.gasvol_s/u.liqvol_s];
            for t = 1 : numel(props.(key))
                props.(key){t} = convertFrom(props.(key){t}, unt);
            end

         case {'OILVISCT', 'WATERVICT'}
            unt         = [u.temp, u.viscosity];
            for t = 1 : numel(props.(key))
                props.(key){t} = convertFrom(props.(key){t}, unt);
            end

        case {'SPECROCK'} 
            unt         = [u.temp, u.volumheatcapacity];
            for t = 1 : numel(props.(key))
                props.(key){t} = convertFrom(props.(key){t}, unt);
            end

         case 'MW'
            unt         = u.mass / u.mol;
            props.(key) = convertFrom(props.(key), unt);

        case 'PLYMAX'
            unt         = [u.concentr, u.concentr];
            props.(key) = convertFrom(props.(key), unt);

         case 'PARACHOR'
            % Always the same units
            unt         = (((centi*meter)^3)/mol)*(dyne/(centi*meter)).^(-1/4);
            props.(key) = convertFrom(props.(key), unt);

         case 'PLYROCK'
            unt         = [1, 1, u.concentr, 1, 1];
            props.(key) = convertFrom(props.(key), unt);

         case {'PLYVISC', 'PLYADS'}
            unt         = [u.concentr, 1];
            for t = 1 : numel(props.(key))
                props.(key){t} = convertFrom(props.(key){t}, unt);
            end

         case 'PLYSHEAR'
            unt = [u.length/u.time, 1];
            for t = 1 : numel(props.(key))
               props.(key){t} = convertFrom(props.(key){t}, unt);
            end

         case 'PLYSHLOG'
             unt = [u.concentr, u.concentr, 1];
             % the last one should be TEMPERATURE, while we DONOT support
             % unit conversion for temperature at the moment.
             props.(key).refcondition = convertFrom(props.(key).refcondition, unt);
             if isfield(props, 'SHRATE')
                 unt = [1/u.time, 1];
             else
                 unt = [u.length/u.time, 1];
             end
             for t = 1 : numel(props.(key).data)
                props.(key).data{t} = convertFrom(props.(key).data{t}, unt);
             end

         case { 'PCRIT', 'PCW', 'IPCW', 'PCG', 'IPCG' }
            unt         = u.press;
            props.(key) = convertFrom(props.(key), unt);

         case 'PVCDO'
            unt         = [u.press, 1, u.compr, u.viscosity, u.compr];
            props.(key) = convertFrom(props.(key), unt);

         case {'PVDG', 'PVDS'}
            unt = [u.press, u.gasvol_r/u.gasvol_s, u.viscosity];
            for t = 1 : numel(props.(key))
               props.(key){t} = convertFrom(props.(key){t}, unt);
            end

         case 'PVDO'
            assert (iscell(props.(key)));

            unt = [u.press, u.liqvol_r/u.liqvol_s, u.viscosity];
            for t = 1 : numel(props.(key))
               props.(key){t} = convertFrom(props.(key){t}, unt);
            end

         case 'PVTG'
            assert (iscell(props.(key)));

            uk = u.press;
            ud = [u.liqvol_s/u.gasvol_s, ...
                  u.gasvol_r/u.gasvol_s, u.viscosity];

            for t = 1 : numel(props.(key))
               props.(key){t}.key  = convertFrom(props.(key){t}.key , uk);
               props.(key){t}.data = convertFrom(props.(key){t}.data, ud);
            end

         case 'PVTO'
            uk = u.gasvol_s / u.liqvol_s;
            ud = [u.press, u.liqvol_r/u.liqvol_s, u.viscosity];

            for t = 1 : numel(props.(key))
               props.(key){t}.key  = convertFrom(props.(key){t}.key , uk);
               props.(key){t}.data = convertFrom(props.(key){t}.data, ud);
            end

         case 'PVTW'
            unt         = [u.press, u.liqvol_r/u.liqvol_s, ...
                           u.compr, u.viscosity, u.compr];
            props.(key) = convertFrom(props.(key), unt);

         case 'RSCONSTT'
            unt         = [u.gasvol_s/u.liqvol_s, u.press];
            props.(key) = convertFrom(props.(key), unt);

          case {'RTEMP', 'RTEMPA'}
            d =  convertFrom(props.(key) + u.tempoffset, u.temp);
            if isnan(d)
                d = 288.706*Kelvin();
            end
            props.(key) = d;

         case 'STCOND'
            d = props.(key);
            d(1) = convertFrom(d(1) + u.tempoffset, u.temp);
            d(2) = convertFrom(d(2), u.press);
            def = [288.7100, 1.01325*barsa];
            d(isnan(d)) = def(isnan(d));
            props.(key) = d;

         case {'SGFN', 'SWFN'}
            unt = [1, 1, u.press];

            for t = 1 : numel(props.(key))
               props.(key){t} = convertFrom(props.(key){t}, unt);
            end

         case {'SGOF', 'SWOF', 'SGWFN', 'SLGOF'}
            assert (iscell(props.(key)));

            unt = [1, 1, 1, u.press];

            for t = 1 : numel(props.(key))
               props.(key){t} = convertFrom(props.(key){t}, unt);
            end
            
         case 'PMISC'
            unt = [u.press, 1];
            for t = 1 : numel(props.(key))
               props.(key){t} = convertFrom(props.(key){t}, unt);
            end

         case 'TCRIT'
            unt         = u.temp;
            props.(key) = convertFrom(props.(key), unt);
            
         case {'TEMPVD', 'ZMFVD'}
            for i = 1:numel(props.(key))
                d = props.(key){i};
                if strcmp(key, 'TEMPVD')
                    d(:, 2) = convertFrom(d(:, 2) + u.tempoffset, u.temp);
                end
                d(:, 1) = convertFrom(d(:, 1), u.length);
                props.(key){i} = d;
            end

         case 'VCRIT'
            unt         = u.density / (u.mass*u.mol);
            props.(key) = convertFrom(props.(key), unt);

         case 'ZCRIT'
            unt         = 1;
            props.(key) = convertFrom(props.(key), unt);

         case {'SOF2', 'SOF3', 'STONE', 'STONE1', 'STONE2', ...
               'SIMPLE', 'TLMIXPAR', 'PLMIXPAR', 'SHRATE', ...
               ...
               'KRW'   , 'KRO'    , 'KRG'  ,           ...
               'SWL'   ,            'ISWL' ,           ...
               'SWLX'  , 'SWLX-'  , 'ISWLX', 'ISWLX-', ...
               'SWLY'  , 'SWLY-'  , 'ISWLY', 'ISWLY-', ...
               'SWLZ'  , 'SWLZ-'  , 'ISWLZ', 'ISWLZ-', ...
               'SWLPC' , 'EHYSTR' ,          'ISWLPC', ...
               ...
               'SWCR'  ,            'SWU'  ,           ...
               'SWCRX' , 'SWCRX-' , 'SWUX' , 'SWUX-' , ...
               'SWCRY' , 'SWCRY-' , 'SWUY' , 'SWUY-' , ...
               'SWCRZ' , 'SWCRZ-' , 'SWUZ' , 'SWUZ-' , ...
               'ISWCR' ,            'ISWU' ,           ...
               'ISWCRX', 'ISWCRX-', 'ISWUX', 'ISWUX-', ...
               'ISWCRY', 'ISWCRY-', 'ISWUY', 'ISWUY-', ...
               'ISWCRZ', 'ISWCRZ-', 'ISWUZ', 'ISWUZ-', ...
               ...
               'SGL'   ,            'ISGL' ,           ...
               'SGLX'  , 'SGLX-'  , 'ISGLX', 'ISGLX-', ...
               'SGLY'  , 'SGLY-'  , 'ISGLY', 'ISGLY-', ...
               'SGLZ'  , 'SGLZ-'  , 'ISGLZ', 'ISGLZ-', ...
               'SGLPC' ,                     'ISGLPC', ...
               ...
               'SGCR'  ,            'SGU'  ,           ...
               'SGCRX' , 'SGCRX-' , 'SGUX' , 'SGUX-' , ...
               'SGCRY' , 'SGCRY-' , 'SGUY' , 'SGUY-' , ...
               'SGCRZ' , 'SGCRZ-' , 'SGUZ' , 'SGUZ-' , ...
               'ISGCR' ,            'ISGU' ,           ...
               'ISGCRX', 'ISGCRX-', 'ISGUX', 'ISGUX-', ...
               'ISGCRY', 'ISGCRY-', 'ISGUY', 'ISGUY-', ...
               'ISGCRZ', 'ISGCRZ-', 'ISGUZ', 'ISGUZ-', ...
               ...
               'SOL'   ,            'ISOL' ,           ...
               'SOLX'  , 'SOLX-'  , 'ISOLX', 'ISOLX-', ...
               'SOLY'  , 'SOLY-'  , 'ISOLY', 'ISOLY-', ...
               'SOLZ'  , 'SOLZ-'  , 'ISOLZ', 'ISOLZ-', ...
               ...
               'SOCR'  ,            'SOU'  ,           ...
               'SOCRX' , 'SOCRX-' , 'SOUX' , 'SOUX-' , ...
               'SOCRY' , 'SOCRY-' , 'SOUY' , 'SOUY-' , ...
               'SOCRZ' , 'SOCRZ-' , 'SOUZ' , 'SOUZ-' , ...
               'ISOCR' ,            'ISOU' ,           ...
               'ISOCRX', 'ISOCRX-', 'ISOUX', 'ISOUX-', ...
               'ISOCRY', 'ISOCRY-', 'ISOUY', 'ISOUY-', ...
               'ISOCRZ', 'ISOCRZ-', 'ISOUZ', 'ISOUZ-', ...
               ...
               'SWATINIT', ...
               ...
               'SOWCR' , 'SOGCR'  , ...
               'ISOWCR', 'ISOGCR' , ...
               ...
               'KRW'   ,            'IKRW'   ,             ...
               'KRWX'  , 'KRWX-'  , 'IKRWX'  , 'IKRWX-'  , ...
               'KRWY'  , 'KRWY-'  , 'IKRWY'  , 'IKRWY-'  , ...
               'KRWZ'  , 'KRWZ-'  , 'IKRWZ'  , 'IKRWZ-'  , ...
               ...
               'KRWR'  ,            'IKRWR'  ,             ...
               'KRWXR' , 'KRWRX-' , 'IKRWXR' , 'IKRWRX-' , ...
               'KRWRY' , 'KRWRY-' , 'IKRWRY' , 'IKRWRY-' , ...
               'KRWRZ' , 'KRWRZ-' , 'IKRWRZ' , 'IKRWRZ-' , ...
               ...
               'KRG'   ,            'IKRG'   ,             ...
               'KRGX'  , 'KRGX-'  , 'IKRGX'  , 'IKRGX-'  , ...
               'KRGY'  , 'KRGY-'  , 'IKRGY'  , 'IKRGY-'  , ...
               'KRGZ'  , 'KRGZ-'  , 'IKRGZ'  , 'IKRGZ-'  , ...
               ...
               'KRGR'  ,            'IKRGR'  ,             ...
               'KRGXR' , 'KRGRX-' , 'IKRGXR' , 'IKRGRX-' , ...
               'KRGRY' , 'KRGRY-' , 'IKRGRY' , 'IKRGRY-' , ...
               'KRGRZ' , 'KRGRZ-' , 'IKRGRZ' , 'IKRGRZ-' , ...
               ...
               'KRO'   ,            'IKRO'   ,             ...
               'KROX'  , 'KROX-'  , 'IKROX'  , 'IKROX-'  , ...
               'KROY'  , 'KROY-'  , 'IKROY'  , 'IKROY-'  , ...
               'KROZ'  , 'KROZ-'  , 'IKROZ'  , 'IKROZ-'  , ...
               ...
               'KRORW' ,            'IKRORW' ,             ...
               'KRORWX', 'KRORWX-', 'IKRORWX', 'IKRORWX-', ...
               'KRORWY', 'KRORWY-', 'IKRORWY', 'IKRORWY-', ...
               'KRORWZ', 'KRORWZ-', 'IKRORWZ', 'IKRORWZ-', ...
               ...
               'KRORG' ,            'IKRORG' ,             ...
               'KRORGX', 'KRORGX-', 'IKRORGX', 'IKRORGX-', ...
               'KRORGY', 'KRORGY-', 'IKRORGY', 'IKRORGY-', ...
               'KRORGZ', 'KRORGZ-', 'IKRORGZ', 'IKRORGZ-', ...
               ...
               'SURFCAPD', ...
               ...
               'ACF', 'BIC', 'CNAMES', 'ROCKOPTS', 'EOS', 'PRCORR', ...
               'SSHIFT',...
               ...
               'MISC', 'MSFN', 'SSFN', 'SCALECRS', ...
               ...
               'GRAVITY', ...
               }
            continue;  % Dimensionless

         otherwise
            error(id('PROPS:NoConverter'), ...
                  'No known unit conversion in PROPS for ''%s''.', key);
      end
   end
end

%--------------------------------------------------------------------------

function soln = convertSOLUTION(soln, u)
   for kw = reshape(fieldnames(soln), 1, [])
      key = kw{1};

      switch key
         case 'AQUCT'                                %  7  8  9 10 11
            unt = [1, u.length, u.press u.perm, 1,    ...
                   1./u.press, u.length, u.length, 1, ...
                   1, 1, 0.0, NaN];
            soln.(key) = convertFrom(soln.(key), unt);

         case 'AQUANCON'
            soln.(key)(:,9) = ...
               cellfun(@(c) convertFrom(c, u.length ^ 2), ...
                       soln.(key)(:,9), 'UniformOutput', false);

         case 'AQUFETP'
            unt = [1, u.length, u.press, u.liqvol_s, ...
                  1 / u.press, u.liqvol_s / u.time / u.press, ...
                  1, u.mass / (u.liqvol_s), u.temp];
            
            soln.(key) = convertFrom(soln.(key), unt);

         case 'EQUIL'                                %  7  8  9 10 11
            unt = [repmat([u.length, u.press], [1, 3]), 1, 1, 1, 1, 1];

            unt        = unt(1 : size(soln.(key), 2));
            soln.(key) = convertFrom(soln.(key), unt);

         case 'FIELDSEP'
            unt = [1, u.temp, u.press, 1, 1, 1, 1, 1, u.temp, u.press];
            off = [0, u.tempoffset, 0, 0, 0, 0, 0, 0, u.tempoffset, 0];
            def = [nan, 288.71, 1*atm];
            ncol = size(soln.(key), 2);
            act = 1:ncol;
            unt = unt(act);
            off = off(act);
            soln.(key) = convertFrom(soln.(key) + off, unt);
            for i = 1:min(ncol, 3)
                defaulted = isnan(soln.(key)(:, i));
                soln.(key)(defaulted, i) = def(i);
            end

         case 'DATUM'
            soln.(key) = convertFrom(soln.(key), u.length);

         case {'PBVD', 'PDVD'}
            unt = [u.length, u.press];

            for reg = 1 : numel(soln.(key))
               soln.(key){reg} = convertFrom(soln.(key){reg}, unt);
            end

         case {'PBUB', 'PRESSURE'}
               soln.(key) = convertFrom(soln.(key), u.press);

         case 'RSVD'
            unt = [u.length, u.gasvol_s/u.liqvol_s];

            for reg = 1 : numel(soln.(key))
               soln.(key){reg} = convertFrom(soln.(key){reg}, unt);
            end

         case 'RVVD'
            unt = [u.length, u.liqvol_s/u.gasvol_s];

            for reg = 1 : numel(soln.(key))
               soln.(key){reg} = convertFrom(soln.(key){reg}, unt);
            end

         case 'RS'
            soln.(key) = convertFrom(soln.(key), u.gasvol_s / u.liqvol_s);

         case 'RV'
            soln.(key) = convertFrom(soln.(key), u.liqvol_s / u.gasvol_s);

         case {'SGAS', 'SOIL', 'SWAT', 'XMF', 'YMF', 'ZMF'}
            continue;  % Dimensionless

         case {'RPTRST', 'RPTSOL', 'OUTSOL'}
            continue;  % Reporting controls.  No units.

         case 'TEMPI'
            soln.(key) = convertFrom(soln.(key), u.temp);

         case 'THPRES'
            unt = [1, 1, u.press];
            soln.(key) = convertFrom(soln.(key), unt);

         otherwise
            error(id('SOLUTION:NoConverter'), ...
                  'No known unit conversion in SOLUTION for ''%s''.', key);
      end
   end
end

%--------------------------------------------------------------------------

function schd = convertSCHEDULE(schd, u)
   schd.step.val = convertFrom(schd.step.val, u.time);

   if isfield(schd, 'control')
      for c = 1 : numel(schd.control)
         schd.control(c) = convertControl(schd.control(c), u);
      end
   end
end

%--------------------------------------------------------------------------

function ctrl = convertControl(ctrl, u)
   for kw = reshape(fieldnames(ctrl), 1, [])
      key = kw{1};

      switch key
         case 'COMPDAT'
            if ~isempty(ctrl.(key))
               ctrl.(key) = convertCompDat(ctrl.(key), u);
            end

         case 'COMPSEGS'
            if ~isempty(ctrl.(key))
               ctrl.(key) = convertCompSegs(ctrl.(key), u);
            end

         case 'WCONHIST'
            if ~isempty(ctrl.(key))
               ctrl.(key) = convertWconHist(ctrl.(key), u);
            end

         case 'WCONINJ'
            if ~isempty(ctrl.(key))
               ctrl.(key) = convertWConInj(ctrl.(key), u);
            end

         case 'WCONINJE'
            if ~isempty(ctrl.(key))
               ctrl.(key) = convertWConInje(ctrl.(key), u);
            end

         case 'WCONINJH'
            if ~isempty(ctrl.(key))
               ctrl.(key) = convertWConInjh(ctrl.(key), u);
            end

         case 'WCONPROD'
            if ~isempty(ctrl.(key))
               ctrl.(key) = convertWConProd(ctrl.(key), u);
            end

         case 'WPOLYMER'
            if ~isempty(ctrl.(key))
               ctrl.(key) = convertWPolymer(ctrl.(key), u);
            end

         case 'WSURFACT'
            if ~isempty(ctrl.(key))
               ctrl.(key) = convertWSurfact(ctrl.(key), u);
            end
            
         case 'WSOLVENT'
            if ~isempty(ctrl.(key))
               ctrl.(key) = convertWSolvent(ctrl.(key), u);
            end

         case 'GCONINJE'
            if ~isempty(ctrl.(key))
               ctrl.(key) = convertGconInje(ctrl.(key), u);
            end

         case 'GCONPROD'
            if ~isempty(ctrl.(key))
               ctrl.(key) = convertGConProd(ctrl.(key), u);
            end

         case 'GECON'
            if ~isempty(ctrl.(key))
               ctrl.(key) = convertGEcon(ctrl.(key), u);
            end

         case 'WELSPECS'
            if ~isempty(ctrl.(key))
               ctrl.(key) = convertWellSpecs(ctrl.(key), u);
            end

         case 'WELSEGS'
            if ~isempty(ctrl.(key))
               ctrl.(key) = convertWellSegs(ctrl.(key), u);
            end

         case {'GRUPTREE', 'WGRUPCON', 'VAPPARS', ...
               'RPTSCHED', 'RPTRST', 'OUTSOL'}
            continue; % No conversion necessary

         case 'GRUPNET'
            if ~isempty(ctrl.(key))
               ctrl.(key) = convertGrupNet(ctrl.(key), u);
            end

         case 'WTEMP'
            if ~isempty(ctrl.(key))
               ctrl.(key) = convertWTemp(ctrl.(key), u);
            end

         case 'DRSDT'
            unt           = (u.gasvol_s / u.liqvol_s) / u.time;
            ctrl.(key){1} = convertFrom(ctrl.(key){1}, unt);

         case { 'VFPINJ', 'VFPPROD' }
            ctrl.(key) = feval(['convert', key], ctrl.(key), u);

         otherwise
            error(id('SCHEDULE:NoConverter'), ...
                  'No known unit conversion in SCHEDULE for ''%s''.', key);
      end
   end
end

%--------------------------------------------------------------------------

function compdat = convertCompDat(compdat, u)
   c = [8      , 9       , 10             , 14      ];
   u = [u.trans, u.length, u.perm*u.length, u.length];

   for n = 1 : numel(c)
      compdat(:, c(n)) = cellfun(@(x) convertFrom(x, u(n)), ...
                                 compdat(:, c(n)), 'UniformOutput', false);
   end
end

%--------------------------------------------------------------------------

function compsegs = convertCompSegs(compsegs, u)
   %     Start     End       Depth     Thermal Length
   c = [ 5       , 6       , 9       , 10       ];
   u = [ u.length, u.length, u.length, u.length ];

   for w = 1 : size(compsegs, 1)
      for n = 1 : numel(c)
         compsegs{w, 2}(:, c(n)) = ...
            cellfun(@(x) convertFrom(x, u(n)), ...
                    compsegs{w, 2}(:, c(n)), ...
                    'UniformOutput', false);
      end
   end
end

%--------------------------------------------------------------------------

function wch = convertWconHist(wch, u)
   lrat = u.liqvol_s / u.time;
   grat = u.gasvol_s / u.time;

   c = [4   , 5   , 6   , 9      , 10     , 11  ];
   u = [lrat, lrat, grat, u.press, u.press, grat];

   c = c(c <= size(wch, 2));

   for n = 1 : numel(c)
      wch(:, c(n)) = cellfun(@(x) convertFrom(x, u(n)), ...
                             wch(:, c(n)), 'UniformOutput', false);
   end
end

%--------------------------------------------------------------------------

function wconinj = convertWConInj(wconinj, u)
   % Surface flow rate
   is_gas = strcmpi(wconinj(:, 2), 'gas');
   is_oil = strcmpi(wconinj(:, 2), 'oil');

   r_unit         = repmat(u.liqvol_s / u.time, [numel(is_gas), 1]);
   r_unit(is_gas) = u.gasvol_s / u.time;

   wconinj(:, 5)  = arrayfun(@convertFrom,                   ...
                             vertcat(wconinj{:, 5}), r_unit, ...
                             'UniformOutput', false);

   % Vapourised gas concentration
   if any(is_gas | is_oil)
      R_unit = repmat(u.liqvol_s / u.gasvol_s, [numel(is_gas), 1]);
      R_unit(is_oil) = u.gasvol_s / u.liqvol_s;

      wconinj(:, 12) = arrayfun(@convertFrom,                    ...
                                vertcat(wconinj{:, 12}), R_unit, ...
                                'UniformOutput', false);
   end

   % Handle default pressure
   i = vertcat(wconinj{:,9}) == -1;
   wconinj(i,9) = { 100e3*psia };
   wconinj(~i, 9) = cellfun(@(x) convertFrom(x, u.press), ...
                            wconinj(~i,9), 'UniformOutput', false);

   % General case
   c = [6                , 7, 10     ];
   u = [u.liqvol_r/u.time, 1, u.press];

   for n = 1 : numel(c)
      wconinj(:, c(n)) = cellfun(@(x) convertFrom(x, u(n)), ...
                                 wconinj(:, c(n)), 'UniformOutput', false);
   end
end

%--------------------------------------------------------------------------

function wconinj = convertWConInje(wconinj, u)
   assert (size(wconinj, 2) == 14, ...
           'WCONINJE data must consist of 14 items per record.');

   % Surface flow rate
   is_gas = strcmpi(wconinj(:, 2), 'gas');
   is_oil = strcmpi(wconinj(:, 2), 'oil');

   r_unit         = repmat(u.liqvol_s / u.time, [numel(is_gas), 1]);
   r_unit(is_gas) = u.gasvol_s / u.time;

   wconinj(:, 5)  = arrayfun(@convertFrom,                   ...
                             vertcat(wconinj{:, 5}), r_unit, ...
                             'UniformOutput', false);

   % Vapourised gas concentration
   if any(is_gas | is_oil)
      R_unit = repmat(u.liqvol_s / u.gasvol_s, [numel(is_gas), 1]);
      R_unit(is_oil) = u.gasvol_s / u.liqvol_s;

      wconinj(:, 10) = arrayfun(@convertFrom,                    ...
                                vertcat(wconinj{:, 10}), R_unit, ...
                                'UniformOutput', false);
   end

   % Handle default pressure
   i = isnan(vertcat(wconinj{:,7}));
   wconinj( i,7) = { 100e3*psia };
   wconinj(~i,7) = cellfun(@(x) convertFrom(x, u.press), ...
                           wconinj(~i,7), 'UniformOutput', false);

   % General case.
   % Items 12:14 are dimensionless numbers.  Just skip them.
   %
   c = [6                , 8      , 11                   ];
   u = [u.liqvol_r/u.time, u.press, u.gasvol_s/u.liqvol_s];

   for n = 1 : numel(c)
      wconinj(:, c(n)) = cellfun(@(x) convertFrom(x, u(n)), ...
                                 wconinj(:, c(n)), 'UniformOutput', false);
   end
end

%--------------------------------------------------------------------------

function wconinj = convertWConInjh(wconinj, u)
   % Surface flow rate
   is_gas = strcmpi(wconinj(:, 2), 'gas');
   is_oil = strcmpi(wconinj(:, 2), 'oil');

   r_unit         = repmat(u.liqvol_s / u.time, [numel(is_gas), 1]);
   r_unit(is_gas) = u.gasvol_s / u.time;

   wconinj(:, 4)  = arrayfun(@convertFrom,                   ...
                             vertcat(wconinj{:, 4}), r_unit, ...
                             'UniformOutput', false);

   % Vapourised gas concentration
   if any(is_gas | is_oil)
      R_unit = repmat(u.liqvol_s / u.gasvol_s, [numel(is_gas), 1]);
      R_unit(is_oil) = u.gasvol_s / u.liqvol_s;

      wconinj(:, 8) = arrayfun(@convertFrom,                    ...
                                vertcat(wconinj{:, 8}), R_unit, ...
                                'UniformOutput', false);
   end

   % General case
   c = [5      , 6      ];
   u = [u.press, u.press];

   for n = 1 : numel(c)
      wconinj(:, c(n)) = cellfun(@(x) convertFrom(x, u(n)), ...
                                 wconinj(:, c(n)), 'UniformOutput', false);
   end
end

%--------------------------------------------------------------------------

function wcp = convertWConProd(wcp, u)
   resv = u.liqvol_r / u.time;
   lrat = u.liqvol_s / u.time;
   grat = u.gasvol_s / u.time;

   c = [4   , 5   , 6   , 7   , 8   , 9      , 10     ];
   u = [lrat, lrat, grat, lrat, resv, u.press, u.press];

   for n = 1 : numel(c)
      wcp(:, c(n)) = cellfun(@(x) convertFrom(x, u(n)), ...
                             wcp(:, c(n)), 'UniformOutput', false);
   end
end

%--------------------------------------------------------------------------

function wcp = convertWPolymer(wcp, u)

   c = [2        , 3        ];
   u = [u.density, u.density];

   for n = 1 : numel(c)
      wcp(:, c(n)) = cellfun(@(x) convertFrom(x, u(n)), ...
                             wcp(:, c(n)), 'UniformOutput', false);
   end
end

%--------------------------------------------------------------------------

function wcp = convertWSurfact(wcp, u)

   wcp(:, 2) = cellfun(@(c) convertFrom(c, u.concentr), ...
                       wcp(:, 2), 'UniformOutput', false);
end

%--------------------------------------------------------------------------

function wcp = convertWSolvent(wcp, varargin)

   wcp(:, 2) = cellfun(@(c) convertFrom(c, 1), ...
                       wcp(:, 2), 'UniformOutput', false);
end

%--------------------------------------------------------------------------

function wcp = convertWTemp(wcp, u)

   wcp(:, 2) = cellfun(@(c) convertFrom(c, u.temp), ...
                       wcp(:, 2), 'UniformOutput', false);
end

%--------------------------------------------------------------------------

function wcp = convertGConProd(wcp, u)
   resv = u.liqvol_r / u.time;
   lrat = u.liqvol_s / u.time;
   grat = u.gasvol_s / u.time;

   % Skip "calorific rate target or upper limit" (item 17).
   % We're isothermal...
   c = [3   , 4   , 5   , 6   , 14  , 15  , 16  , 17 , 18  , 19   ];
   u = [lrat, lrat, grat, lrat, resv, resv, grat, inf, grat, lrat ];

   for n = 1 : numel(c)
      wcp(:, c(n)) = cellfun(@(x) convertFrom(x, u(n)), ...
                             wcp(:, c(n)), 'UniformOutput', false);
   end
end

%--------------------------------------------------------------------------

function gconinj = convertGconInje(gconinj, u)
   % Surface flow rate
   is_gas = strcmpi(gconinj(:, 2), 'gas');

   r_unit         = repmat(u.liqvol_s / u.time, [numel(is_gas), 1]);
   r_unit(is_gas) = u.gasvol_s / u.time;

   gconinj(:, 4)  = arrayfun(@convertFrom,                   ...
                             vertcat(gconinj{:, 4}), r_unit, ...
                             'UniformOutput', false);

   % General case
   c = [5         , 13        ];
   u = [u.liqvol_r, u.gasvol_s] ./ u.time;

   for n = 1 : numel(c)
      gconinj(:, c(n)) = cellfun(@(x) convertFrom(x, u(n)), ...
                                 gconinj(:, c(n)), 'UniformOutput', false);
   end
end

%--------------------------------------------------------------------------

function gecon = convertGEcon(gecon, u)
   orat = u.liqvol_s / u.time;
   grat = u.gasvol_s / u.time;

   c = [ 2    , 3    ];
   u = [ orat , grat ];

   for n = 1 : numel(c)
      gecon(:, c(n)) = cellfun(@(x) convertFrom(x, u(n)), ...
                               gecon(:, c(n)), 'UniformOutput', false);
   end
end

%--------------------------------------------------------------------------

function grpnet = convertGrupNet(grpnet, u)
   % Only item 2 needs unit conversion.  The remaining items are strings or
   % pure numbers.
   grpnet(:, 2) = cellfun(@(x) convertFrom(x, u.press), grpnet(:, 2), ...
                          'UniformOutput', false);
end

%--------------------------------------------------------------------------

function wspecs = convertWellSpecs(wspecs, u)
   c = [5       , 7       ];
   u = [u.length, u.length];

   for n = 1 : numel(c)
      wspecs(:, c(n)) = cellfun(@(x) convertFrom(x, u(n)), ...
                                wspecs(:, c(n)), 'UniformOutput', false);
   end
end

%--------------------------------------------------------------------------

function wsegs = convertWellSegs(wsegs, u)
   for w = 1 : size(wsegs, 1)
      wsegs{w,2}.header   = convertWellSegsHeader (wsegs{w,2}.header, u);
      wsegs{w,2}.segments = convertWellSegsRecords(wsegs{w,2}.segments, u);
   end
end

%--------------------------------------------------------------------------

function header = convertWellSegsHeader(header, u)
   % Note: Subtract 1 from columns for well name stored outside of header

   %     Depth     Length    Eff. Vol  X coord   Y coord
   c = [ 2       , 3       , 4       , 8       , 9        ] - 1;
   u = [ u.length, u.length, u.volume, u.length, u.length ];

   for n = 1 : numel(c)
      header(c(n)) = cellfun(@(x) convertFrom(x, u(n)), ...
                             header(c(n)), 'UniformOutput', false);
   end
end

%--------------------------------------------------------------------------

function segments = convertWellSegsRecords(segments, u)
   %     Length    Depth     Int. Diam Roughness
   c  = [ 5       , 6       , 7       , 8       ];
   uu = [ u.length, u.length, u.length, u.length];

   %         Cross   Volume    X coord   Y coord   Wall Area
   c  = [c , 9     , 10      , 11      , 12      , 13    ];
   uu = [uu, u.area, u.volume, u.length, u.length, u.area];

   for n = 1 : numel(c)
      segments(:, c(n)) = ...
         cellfun(@(x) convertFrom(x, uu(n)), ...
                 segments(:, c(n)), 'UniformOutput', false);
   end
end

%--------------------------------------------------------------------------

function vfpinj = convertVFPINJ(vfpinj, u)
   for i = 1 : numel(vfpinj)
      if isempty(vfpinj{i}), continue, end

      [vfpinj{i}, usys] = convertVFPCommon(vfpinj{i}, u);

      vfpinj{i}.BHP = convertFrom(vfpinj{i}.BHP, usys.press);
   end
end

%--------------------------------------------------------------------------

function vfpprod = convertVFPPROD(vfpprod, u)
   for i = 1 : numel(vfpprod)
      if isempty(vfpprod{i}), continue, end

      if ~strcmp(vfpprod{i}.QID, 'BHP')
         error('VFPProd:Unsupp:Quant', ...
              ['Physical quantity ''%s'' not supported in VFPPROD ', ...
               'table %d/%d'], vfpprod{i}.QID, i, numel(vfpprod));
      end

      [vfpprod{i}, usys] = convertVFPCommon(vfpprod{i}, u);

      % Note: GFR, WFR are dimensionless ratios (with the exception of
      % molar weights, which are not supported here presently).

      if ~ ((numel(vfpprod{i}.ALQ) == 1) && (vfpprod{i}.ALQ == 0))
         switch vfpprod{i}.ALQID
            case {'', ' '},        alq_unit = 1;
            case 'GRAT',           alq_unit = usys.gasvol_s / usys.time;
            case 'IGLR',           alq_unit = 1;
            case 'TGLR',           alq_unit = 1;
            case {'DENG', 'DENO'}, alq_unit = usys.density;
            case 'BEAN',           alq_unit = usys.length;
            case 'PUMP',           error('Pump rating not known');
            case 'COMP',           error('Compressor unit not known');
            otherwise,             alq_unit = 1;
         end

         vfpprod{i}.ALQ = convertFrom(vfpprod{i}.ALQ, alq_unit);
      end

      vfpprod{i}.Q = convertFrom(vfpprod{i}.Q, usys.press);
   end
end

%--------------------------------------------------------------------------

function [vfp, usys] = convertVFPCommon(vfp, u)
   if strcmp(vfp.USYS, 'USYS')
      usys = u;
   else
      usys = unitConversionFactors(vfp.USYS);
   end

   if any(strcmp(vfp.FLOID, {'OIL', 'WAT', 'LIQ'}))
      funit = usys.liqvol_s;
   elseif strcmp(vfp.FLOID, 'GAS')
      funit = usys.gasvol_s;
   else
      error('Unsupported FLOID: %s', vfp.FLOID);
   end

   vfp.depth = convertFrom(vfp.depth, usys.length);
   vfp.THP   = convertFrom(vfp.THP  , usys.press);
   vfp.FLO   = convertFrom(vfp.FLO  , funit);
   vfp.USYS  = 'SI';
end

%--------------------------------------------------------------------------

function s = id(s)
   s = ['convertDeckUnits:', s];
end
