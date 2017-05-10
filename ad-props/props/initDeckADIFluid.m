function fluid = initDeckADIFluid(deck, varargin)


 % this only work for full deck or first region.
 opt.G = [];
 opt.set_method='from_deck';
 opt.regs=struct('PVTNUM',  1, ...
                 'SATNUM',  1,...
                 'SURFNUM', 1,...
                 'IMBNUM',  1,...
                 'ROCKNUM', 1);
 opt = merge_options(opt,varargin{:});

 reg.PVTNUM = [];
 reg.PVTINX = ':';
 reg.SURFNUM = [];
 reg.SURFINX = ':';
 reg.SATINX=':';
 reg.SATNUM = [];
 reg.SATINX = ':';
 reg.IMBNUM = [];
 reg.IMBINX = ':';
 reg.ROCKNUM = [];
 reg.ROCKINX = ':';
 switch opt.set_method
     case 'from_deck'
         reg = handleRegions(deck,opt.G);%, varargin{:});
     case 'first_region'
         %
     case 'given_regions_ongrid'
         %
         nc=G.cells.num;
         oo=ones(nc,1);
         fname=fieldnames(reg);
         for i=1:numel(fname)
             if(isfield(opt.regs,fname{i}))
                 if(numel(opts.regs.(fname{i}))==nc)
                     reg.(fname{i})= opt.regs.(fname{i})*oo;
                 else
                     reg.(fname{i})= opt.regs.(fname{i});
                 end
             end
         end
     case 'one_region'
         %
         fname=fieldnames(reg);
         for i=1:numel(fname)
             if(isfield(opt.regs,fname{i}))
                 assert(numel(opt.regs.(fname{i}))==1);
                 reg.(fname{i})= opt.regs.(fname{i});
             end
         end
         
     otherwise
         error('No surch region method')
 end
%reg = handleRegions(deck,opt.G);%, varargin{:});
fluid = [];
%props
props = deck.PROPS;
fns = fieldnames(props);
for k = 1:numel(fns)
    fn   = fns{k};
    if doAssign(fn)
        asgn = str2func(['assign',fn]);
        try
            fluid = asgn(fluid, props.(fn), reg);
        catch  %#ok
            warning(msgid('Assign:Failed'), ...
                'Could not assign property ''%s''.', fn)
        end
    end
end
fluid = assignRelPerm(fluid);
end

function flag = doAssign(propNm)
% Properties not resulting in individual functions
excpt = {'SWL'   ,'SWCR'   ,'SWU' , ...
         'SGL'   ,'SGCR'   ,'SGU' , ...
         'SOWCR' ,'SOGCR'  , ...
         'CNAMES', 'BIC'   , 'ACF', ...
         'PCRIT' , 'TCRIT' , 'VCRIT',...
         'MW', 'ZCRIT', ...
         'ISWL'  ,'ISWCR'  ,'ISWU', ...
         'ISGL'  ,'ISGCR'  ,'ISGU', ...
         'ISOWCR','ISOGCR'};
flag = ~any( strcmp(propNm , excpt) );
end

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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
