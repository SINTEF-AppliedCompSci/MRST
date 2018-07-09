function deck = readPROPS(fid, dirname, deck)
% deck = readPROPS(fid, dirname, deck)
%
% Otherwise intentionally undocumented.


%{
Copyright 2009-2018 SINTEF ICT, Applied Mathematics.

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

   [dims, ntsfun, ntpvt, ntmisc, ntrocc, ncomp] = get_dimensions(deck);

   [prp, miss_kw] = get_state(deck);

   kw = getEclipseKeyword(fid);
   in_section = ischar(kw);
   while in_section
      if isfield(prp, kw)
         error('Keyword ''%s'' is already defined.', kw);
      end

      switch kw
         case 'ACF'
            prp.(kw) = readVector(fid, kw, ncomp);

         case 'BOX'
            boxKeyword(fid);

         case 'BIC'
            prp.(kw) = readVector(fid, kw, ncomp*(ncomp-1)/2);

         case 'CNAMES'
            tmpl = arrayfun(@(x) sprintf('COMP_%d', x), 1:ncomp, ...
                                                'UniformOutput', false);
            prp.(kw) = readDefaultedKW(fid, tmpl, 'NRec', 1); clear tmpl

         case 'ENDBOX'
            endboxKeyword;
            
          case 'EOS'
              prp.(kw) = readDefaultedRecord(fid, {'PR'});

          case 'PRCORR'
              prp.PRCORR = true;
              
         case 'DENSITY'
            tmpl(1:3) = { '0.0' };
            data      = readDefaultedKW(fid, tmpl, 'NRec', ntpvt);
            prp.(kw)  = to_double(data);  clear tmpl

         case 'SDENSITY'
            tmpl(1) = { '0.0' };
            data      = readDefaultedKW(fid, tmpl, 'NRec', ntpvt);
            prp.(kw)  = to_double(data);  clear tmpl

         case 'SURFST'
            prp.(kw) = readImmisciblePVTTable(fid, ntpvt, 2);

         case 'SURFCAPD'
            prp.(kw) = readRelPermTable(fid, kw, ntsfun, 2);

         case 'SURFVISC'
            prp.(kw) = readImmisciblePVTTable(fid, ntpvt, 2);

         case 'SURFADS'
           prp.(kw) = readRelPermTable(fid, kw, ntsfun, 2);

         case 'SURFROCK'
            tmpl(1:2) = { 'NaN' };
            data      = readDefaultedKW(fid, tmpl, 'NRec', ntsfun);
            prp.(kw)  = to_double(data);  clear tmpl

         case {'MISC', 'PMISC'}
            prp.(kw) = readRelPermTable(fid, kw, ntsfun, 2);

         case 'MSFN'
            prp.(kw) = readRelPermTable(fid, kw, ntsfun, 3);

         case {'MW', 'PCRIT', 'TCRIT', 'VCRIT', 'ZCRIT'}
            prp.(kw) = readVector(fid, kw, ncomp);

         case 'PLYADS'
            prp.(kw) = readRelPermTable(fid, kw, ntsfun, 2);

         case 'PLYMAX'
            tmpl(1:2) = { '0.0' };
            data      = readDefaultedKW(fid, tmpl, 'NRec', ntmisc);
            prp.(kw)  = to_double(data);  clear tmpl

         case 'PLYROCK'
            tmpl(1:5) = { 'NaN' };
            data      = readDefaultedKW(fid, tmpl, 'NRec', ntsfun);
            prp.(kw)  = to_double(data);  clear tmpl

         case 'PLYVISC'
            prp.(kw) = readImmisciblePVTTable(fid, ntpvt, 2);

         case 'PLYSHEAR'
            prp.(kw) = readImmisciblePVTTable(fid, ntpvt, 2);

         case 'PLYSHLOG'
             tmpl(1:3) = { 'NaN' };
             refcondition = readDefaultedKW(fid, tmpl, 'NRec', ntpvt);
             prp.(kw).refcondition = to_double(refcondition);
             clear tmpl refcondition;
             data = readImmisciblePVTTable(fid, ntpvt, 2);
             prp.(kw).data = data;
             clear data;

         case 'SHRATE'
             tmpl = { 'NaN' };
             data = readDefaultedKW(fid, tmpl, 'NRec', ntpvt);
             if ~isempty(data)
                 data = to_double(data);
             else
                 data = repmat(4.8, [ ntpvt, 1 ]);
             end
             data(isnan(data)) = 4.8;
             prp.(kw) = data; clear tmpl; clear data;

         case 'PVCDO'
            tmpl(1:5) = { '0.0' };
            data      = readDefaultedKW(fid, tmpl, 'NRec', ntpvt);
            prp.(kw)  = to_double(data);  clear tmpl

         case {'PVDG', 'PVDO', 'PVDS'}
            prp.(kw) = readImmisciblePVTTable(fid, ntpvt, 3);

         case 'PVCO'
            tbl = readImmisciblePVTTable(fid, ntpvt, 6);
            assert(~isfield(prp, 'PVTO'));
            prp.PVTO = expandPVCOintoPVTO(tbl, ntpvt);

         case {'PVTG', 'PVTO'}
            prp.(kw) = readMisciblePVTTable(fid, ntpvt, 3, kw);
%            prp.(kw) = feval(['complete', kw], prp.(kw));

         case 'PVTW'
            tmpl(1:5) = { '0.0' };
            data      = readDefaultedKW(fid, tmpl, 'NRec', ntpvt);
            prp.(kw)  = to_double(data);  clear tmpl

         case 'RKTRMDIR'
            prp.(kw) = true;

         case 'ROCKOPTS'
            tmpl = { 'PRESSURE', 'NOSTORE', 'PVTNUM', 'DEFLATION' };
            prp.(kw) = readDefaultedKW(fid, tmpl, 'NRec', 1);    clear tmpl

         case 'ROCK'
            nrec = ntpvt;

            if isfield(prp, 'ROCKOPTS')
               nrec = ntrocc;
            end

            tmpl(1:6) = { 'NaN' };
            data      = readDefaultedKW(fid, tmpl, 'NRec', nrec);
            prp.(kw)  = to_double(data);  clear tmpl nrec

         case 'SPECHEAT'
            prp.(kw) = readRelPermTable(fid, kw, ntsfun, 4);

         case 'SPECROCK'
            prp.(kw) = readRelPermTable(fid, kw, ntsfun, 2);

         case 'ROCKTAB'
            % Number of columns is 5 if RKTRMDIR is specified, 3 otherwise.
            ncol     = 3 + 2*isfield(prp, 'RKTRMDIR');
            prp.(kw) = readRelPermTable(fid, kw, ntrocc, ncol);  clear ncol

         case 'RSCONSTT'
            tmpl     = { '-1.0', '-1.0' };
            data     = readDefaultedKW(fid, tmpl, 'NRec', ntpvt);
            prp.(kw) = to_double(data);   clear tmpl

         case 'SOF2'
            prp.(kw) = readRelPermTable(fid, kw, ntsfun, 2);

         case {'SGFN', 'SWFN', 'SOF3', 'SSFN'}
            prp.(kw) = readRelPermTable(fid, kw, ntsfun, 3);

         case {'SGOF', 'SGWFN', 'SLGOF', 'SWOF'}
            prp.(kw) = readRelPermTable(fid, kw, ntsfun, 4);

         case {'STONE' , 'STONE1', 'STONE2', 'SIMPLE'}
            prp.(kw) = true;

         case {'SWL'   ,            'ISWL' ,           ...
               'SWLX'  , 'SWLX-'  , 'ISWLX', 'ISWLX-', ...
               'SWLY'  , 'SWLY-'  , 'ISWLY', 'ISWLY-', ...
               'SWLZ'  , 'SWLZ-'  , 'ISWLZ', 'ISWLZ-', ...
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
               'SOCR'  ,            'SOU'  ,           ...
               'SOCRX' , 'SOCRX-' , 'SOUX' , 'SOUX-' , ...
               'SOCRY' , 'SOCRY-' , 'SOUY' , 'SOUY-' , ...
               'SOCRZ' , 'SOCRZ-' , 'SOUZ' , 'SOUZ-' , ...
               'ISOCR' ,            'ISOU' ,           ...
               'ISOCRX', 'ISOCRX-', 'ISOUX', 'ISOUX-', ...
               'ISOCRY', 'ISOCRY-', 'ISOUY', 'ISOUY-', ...
               'ISOCRZ', 'ISOCRZ-', 'ISOUZ', 'ISOUZ-', ...
               ...
               'SOWCR' , 'SOGCR'  , ...
               'ISOWCR', 'ISOGCR' , ...
               ...
               'SWATINIT'
               }
            prp = readGridBoxArray(prp, fid, kw, prod(dims));

         case {'TLMIXPAR', 'PLMIXPAR'}
            % Note: The following warning is given in ECLIPSE 2013.2 when
            % running a polymer model with TLMIXPAR: THE TLMIXPAR KEYWORD
            % IS NO LONGER COMPATIBLE WITH THE POLYMER FLOOD MODEL. USE THE
            % PLMIXPAR KEYWORD INSTEAD.
            tmpl(1:2) = { 'NaN' };
            data      = readDefaultedKW(fid, tmpl, 'NRec', ntmisc);
            prp.(kw)  = to_double(data);  clear tmpl

         case {'ADD', 'COPY', 'EQUALS', 'MAXVALUE', ...
               'MINVALUE', 'MULTIPLY'}
            prp = applyOperator(prp, fid, kw);

         case {'ECHO', 'NOECHO'}
            kw = getEclipseKeyword(fid);
            continue;  % Ignore.  Not handled in MRST

         %-----------------------------------------------------------------
         % Sectioning keywords below.  Modifies flow of control.
         % Don't change unless absolutely required...
         %
         case 'END'
            % Logical end of input deck.
            %
            in_section = false;
            deck.PROPS = prp;

            % Restore default input box at end of section
            gridBox(defaultBox);

         case 'REGIONS'
            % Read next section (i.e., 'PROPS' -> 'REGIONS' -> 'SOLUTION')
            in_section = false;

            deck = set_state(deck, prp, miss_kw);

            % Restore default input box at end of section
            gridBox(defaultBox);

            deck = readREGIONS(fid, dirname, deck);

         case 'SOLUTION'
            % Read next section (i.e., 'PROPS' -> 'SOLUTION', no 'REGIONS')
            in_section = false;

            deck = set_state(deck, prp, miss_kw);

            % Restore default input box at end of section
            gridBox(defaultBox);

            deck = readSOLUTION(fid, dirname, deck);

         case 'INCLUDE'
            % Handle 'INCLUDE' (recursion).
            deck = set_state(deck, prp, miss_kw);

            deck = readEclipseIncludeFile(@readPROPS, fid, dirname, deck);

            % Prepare for additional reading.
            [prp, miss_kw] = get_state(deck);

          otherwise
            if ischar(kw)
               miss_kw = [ miss_kw, { kw } ];  %#ok
            end
      end

      % Get next keyword.
      kw = getEclipseKeyword(fid);
      in_section = in_section && ischar(kw);
   end

   deck = set_state(deck, prp, miss_kw);
end

%--------------------------------------------------------------------------

function [dims, ntsfun, ntpvt, ntmisc, ntrocc, ncomp] = get_dimensions(deck)
   assert (isstruct(deck) && isfield(deck, 'RUNSPEC') && ...
           isstruct(deck.RUNSPEC));

   dims = reshape(deck.RUNSPEC.cartDims, 1, []);

   [ntsfun, ntpvt, ntmisc, ncomp] = deal(1);
   ntrocc = -1;

   if isfield(deck.RUNSPEC, 'TABDIMS')
      ntsfun = deck.RUNSPEC.TABDIMS(1);  assert (ntsfun >= 1);
      ntpvt  = deck.RUNSPEC.TABDIMS(2);  assert (ntpvt  >= 1);
      ntrocc = deck.RUNSPEC.TABDIMS(13);
   end

   if isfield(deck.RUNSPEC, 'ROCKCOMP')
      ntrocc = max(deck.RUNSPEC.ROCKCOMP{2}, ntrocc);
      assert (ntrocc >= 1);
   end

   if isfield(deck.RUNSPEC, 'MISCIBLE')
      ntmisc = deck.RUNSPEC.MISCIBLE{1};  assert (ntmisc >= 1);
   end

   if isfield(deck.RUNSPEC, 'COMPS')
      ncomp = deck.RUNSPEC.COMPS;
   end
end

%--------------------------------------------------------------------------

function v = to_double(v)
   convert = @(s) sscanf(regexprep(s, '[dD]', 'e'), '%f');

   if ischar(v)
      v = convert(v);
   else
      v = cellfun(convert, v);
   end
end

%--------------------------------------------------------------------------

function [prp, miss_kw] = get_state(deck)
   prp     = deck.PROPS;
   miss_kw = deck.UnhandledKeywords.PROPS;
end

%--------------------------------------------------------------------------

function deck = set_state(deck, prp, miss_kw)
   deck.PROPS                   = prp;
   deck.UnhandledKeywords.PROPS = unique(miss_kw);
end
