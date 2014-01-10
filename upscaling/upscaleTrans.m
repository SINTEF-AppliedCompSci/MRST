function [HT_cg,T_cg,cgwells,report]=upscaleTrans(cg,T_fine,varargin)
%% Calculate upscaled transmissibilities for a coarse model
%
% SYNOPSIS:
%   [HT_cg,T_cg,cgwells,upscaled]=upscaleTrans(cg,T_fine,varargin)
%
% PARAMETERS:
%   cg       - coarse grid using the new-coarsegrid module structure
%
%   T_fine   - transmissibility on fine grid
%
%   'pn'/pv - List of 'key'/value pairs for supplying optional parameters.
%             The supported options are
%               - Verbose -- Whether or not to display informational
%                            messages throughout the computational process.
%                            Logical.  Default value: Verbose = false
%                            (don't display any informational messages).
%               - bc_method  choise of global boundary condition method.
%                            Valid choises is
%                            'bc_simple' (default),
%                            'bc' need  option 'bc'
%                            'wells_simple', need option 'wells' with one
%                            well configuration and give an upscaled well
%                            cgwells as output
%                            'wells' need 'wells as' cellarray and do not
%                            upscale wells
%             - match_method how to match to the set of global boundary
%                            conditions options are
%                            'max_flux'
%                             lsq_flux'
%             - fix_trans     'true/false' set negative and zero
%                             transmissibility to lowest positive found
%
%
% RETURNS:
%   HT_cg - halfface transmissibilities
%    T_cg - face transmissibilities
%    cgwells - coarse grid well with upscaled well indices
%    upscaled - detailed infomation from the upscaling procedure used for
%               testing
% COMMENTS:
%
%
% SEE ALSO:
%       coarsegrid module, grid_structure
%

%{
Copyright 2009-2014 SINTEF ICT, Applied Mathematics.

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


opt     = struct('verbose',              mrstVerbose, ...
                 'bc_method',        'bc_simple',       ...
                 'bc',[],...
                 'wells',[],...
                 'match_method',       'max_flux',...
                 'fix_trans',           false);
opt     = merge_options(opt, varargin{:});
require agglom;
%if(~isempty(opt.bc))
%   if(~iscell(opt.bc))
%      opt.bc={opt.bc};
%   end
%end
% choose global boundary conditions to use
switch opt.bc_method
   case 'bc_simple'
      % simple pressure in 3 directions
      bc{1} = pside([], cg.parent, 'West', 100*barsa);
      bc{1} = pside(bc{1}, cg.parent, 'East', 200*barsa);
      bc{2} = pside([], cg.parent, 'North', 100*barsa);
      bc{2} = pside(bc{2}, cg.parent, 'South', 200*barsa);
      if(cg.griddim>2)
        bc{3} = pside([], cg.parent, 'Top', 100*barsa);
        bc{3} = pside(bc{3}, cg.parent, 'Bottom', 200*barsa);
      end
      well_cases = cell(numel(bc),1);
      cgwells_cases=[];
   case 'bc'
      % given bc conditions
      if(~isempty(opt.bc))
         bc=opt.bc;
      else
         error('bc_method==bc need non empty bc  option')
      end
      well_cases = cell(numel(bc),1);
      cgwells_cases=[];
   case 'wells_simple'
      % use all linearly independent well configurations for incompressible
      % flow
      if(~isempty(opt.wells))
         num_w=numel(opt.wells);
         A=ones(num_w);
         rates=null(A);
         rates=rates/darcy;
         w_cases=size(rates,2);
         well_cases=cell(w_cases,1);
         bc=cell(w_cases,1);
         for i=1:w_cases
            well_cases{i}=opt.wells;
            for j=1:numel(opt.wells)
               well_cases{i}(j).type='rate';
               well_cases{i}(j).val=rates(j,i);
            end
            % make solution unique
            well_cases{i}(end).type='bhp';
            well_cases{i}(end).val=300*barsa;

         end
      else
         error('bc_method==wells need non empty wells option')
      end
      for i=1:numel(well_cases)
        cgwells_cases{i}=makeCGWells(cg,well_cases{i});
      end
   case 'wells'
      % given wells conditions, need to be valid for incompressible flow
      if(~isempty(opt.wells))
         if(iscell(opt.wells))
            well_cases=opt.wells;
         else
            error('wells option need to be cell array for bc_method=wells')
         end
      else
         error('bc_method==bc need non empty bc  option')
      end
      bc = cell(numel(opt.wells),1);
      %% this is nessesary for now
      %assert(numel(opt.wells)==1)
      %cgwells=makeCGWells(cg,opt.wells);
      for i=1:numel(opt.wells)
        cgwells_cases{i}=makeCGWells(cg,opt.wells{i});
      end
   otherwise
      error('This bc_method is not implemented')
end
% number of cases
nr_global=numel(bc);
% use imple single phase fluid
fluid=initSingleFluid('mu' ,    1*centi*poise     , ...
                      'rho', 0*kilogram/meter^3);
% init structures
state=initState(cg.parent, [], 100*barsa);
upscaled.pressure=nan(cg.cells.num,nr_global);
upscaled.flux=nan(cg.faces.num,nr_global);
upscaled.trans=nan(cg.faces.num,nr_global);
if(isfield(state,'facePressure'))
   upscaled.facePressure=nan(cg.faces.num,nr_global);
   upscaled.htrans=nan(numel(cg.cells.faces(:,1)),nr_global);
end
upscaled=cell(nr_global,1);
flux=nan(cg.faces.num,nr_global);
dpf=nan(cg.faces.num,nr_global);
trans=nan(cg.faces.num,nr_global);
htrans=nan(numel(cg.cells.faces(:,1)),nr_global);
dphf=nan(numel(cg.cells.faces(:,1)),nr_global);
% make coarse grid well structure with out well indices which are to be
% upscaled

% loop over global boundary conditions
for i=1:nr_global
   % finescale solve
   state=incompTPFA(state, cg.parent, T_fine,fluid,'bc',bc{i},'wells',well_cases{i});
   % upscale solution and wellSol
   if(~isempty(cgwells_cases))
      cgwells=cgwells_cases{i};
   else
      cgwells=[];
   end
   upscaled{i} = calculateUpscaledSolution(cg,state,cgwells);%,well_cases{i});
   % calculate effectie trans for the case
   internal = ~any((cg.faces.neighbors==0),2);
   upscaled{i}.trans=nan(cg.faces.num,1);
   upscaled{i}.trans(internal)=...
      upscaled{i}.flux(internal)./(upscaled{i}.pressure(cg.faces.neighbors(internal,2))-...
      upscaled{i}.pressure(cg.faces.neighbors(internal,1)));
   upscaled{i}.trans=-upscaled{i}.trans*fluid.properties();

   if(isfield(state,'facePressure'))
      % calculate half face trans if possible
      cg_cellno=rldecode((1:cg.cells.num)',diff(cg.cells.facePos));
      upscaled{i}.htrans=upscaled{i}.flux(cg.cells.faces(:,1))./...
         ((upscaled{i}.facePressure(cg.cells.faces(:,1))-upscaled{i}.pressure(cg_cellno))...
         .*(2*(cg_cellno==cg.faces.neighbors(cg.cells.faces(:,1),1))-1));
      upscaled{i}.htrans=-fluid.properties()*upscaled{i}.htrans;
      htrans(:,i)=upscaled{i}.htrans;
      dphf(:,i)=(upscaled{i}.facePressure(cg.cells.faces(:,1))-upscaled{i}.pressure(cg_cellno)).*...
         (2*(cg_cellno==cg.faces.neighbors(cg.cells.faces(:,1),1))-1);
   end
   if(~isempty(cgwells))
      % calculate upscaled well indices
     for j=1:numel(cgwells);
        press=repmat(upscaled{i}.wellSol(j).pressure,numel(cgwells(j).cells),1);
        upscaled{i}.WI{j}=(upscaled{i}.wellSol(j).flux)./(press-upscaled{i}.pressure(cgwells(j).cells));
        upscaled{i}.WI{j}=upscaled{i}.WI{j}*fluid.properties();
        fluxwells{j}(:,i)=upscaled{i}.wellSol(j).flux;
        dpwells{j}(:,i)=(press-upscaled{i}.pressure(cgwells(j).cells));
     end
   end
   flux(:,i)=upscaled{i}.flux;
   dpf(internal,i)=upscaled{i}.pressure(cg.faces.neighbors(internal,2))-...
       upscaled{i}.pressure(cg.faces.neighbors(internal,1));
   trans(:,i)=upscaled{i}.trans;
   %% test if caurse solve give same solution
   %assert(isempty(bc{i}));
   if(isempty(bc{i}))
       cgbc=[];
   else
      % Extract subfaces on coarse interfaces
      [nsub, sub]     = subFaces(cg.parent, cg);
      [sgn, coarse_f] = signOfFineFacesOnCoarseFaces(cg.parent, cg, nsub, sub);

      % Coarse boundary conditions
      cgbc = convertBC2Coarse(bc{i},cg.parent, cg, nsub, sub, coarse_f, sgn);


   end
   cgwells_tmp=cgwells;
   for j=1:numel(cgwells_tmp);
       cgwells_tmp(j).WI=upscaled{i}.WI{j};
       cgwells_tmp(j).type=well_cases{i}(j).type;
       cgwells_tmp(j).val=well_cases{i}(j).val;
   end
   cgtrans=upscaled{i}.trans;
   cgtrans(~isfinite(cgtrans))=0;
   if(isfield(state,'facePressure'))
       htrans_tmp=upscaled{i}.htrans;
       %htrans_tmp(abs(upscaled{i}.flux)<1e-6*max(abs(upscaled{i}.flux)))=0
       trans_tmp=accumarray(cg.cells.faces(:,1), 1./upscaled{i}.htrans);
       trans_tmp=1./trans_tmp;
       cgtrans(~internal)=trans_tmp(~internal);
   end
   cgtrans(~isfinite(cgtrans))=min(T_fine)*1e-6;
   %% do not manage to to the check for all cases.
   % use this for debuging
   if(false)
       %cgtrans(abs(upscaled{i}.flux)<1e-6*max(abs(upscaled{i}.flux)))=0;
       cgstate_tmp=struct('flux',upscaled{i}.flux,...
           'pressure',upscaled{i}.pressure,...
           's',zeros(cg.cells.num,1));
       cgstate=incompTPFA(cgstate_tmp, cg, cgtrans, fluid,'bc',cgbc,'wells',cgwells_tmp,'use_trans',true);
       if(max(cgstate.pressure)>min(cgstate.pressure))
           %assert(all(abs((cgstate.pressure-upscaled{i}.pressure)/(max(cgstate.pressure) - min(cgstate.pressure)))<1e-6))
       end
       if(any(upscaled{i}.flux(internal)))
           assert(all(abs(cgstate.flux(internal)-upscaled{i}.flux(internal))./max(abs(cgstate.flux))< 1e-6))
       end
   end
   assert(all(isfinite(cgtrans)))
   upscaled{i}.trans=cgtrans;
   trans(:,i)=cgtrans;

end
% make a choise of which transmisiblities to use on each face. For now use
% the one generated with largest flux
switch opt.match_method
    case 'max_flux'
        % use the transmisibility of the case with largest flux over the face
        [mm,choise]=max(abs(flux),[],2);
        ch_ind=sub2ind(size(htrans),[1:numel(cg.cells.faces(:,1))]',choise(cg.cells.faces(:,1)));
        HT_cg=htrans(ch_ind);
        ch_ind=sub2ind(size(trans),[1:cg.faces.num]',choise);
        T_cg=trans(ch_ind);
        if(~isempty(cgwells))
            for j=1:numel(cgwells);
                bb=cellfun(@(x) x.WI{j},upscaled,'UniformOut',false);
                %wind=horzcat(upscaled{:}.WI{j});
                wind=[bb{:}];
                bb=cellfun(@(x) x.wellSol(j).flux,upscaled,'UniformOut',false);
                %wflux=horzcat(upscaled{:}.wellSol(j).flux);
                wflux=[bb{:}];
                [mm,choise]=max(abs(wflux),[],2);
                ch_ind=sub2ind(size(wind),[1:size(wind,1)]',choise);
                cgwells(j).WI=wind(ch_ind);
            end
        end


    case 'lsq_flux'
        % minimize least square flux error over the set of cases
        HT_cg=-fluid.properties()*sum(flux(cg.cells.faces(:,1),:),2)./sum(dphf,2);
        T_cg =-fluid.properties()*sum(flux,2)./sum(dpf,2);
        if(~isempty(cgwells))
            for j=1:numel(cgwells);
                cgwells(j).WI= sum(fluxwells{j},2)./sum(dpwells{j},2);
                cgwells(j).WI=cgwells(j).WI*fluid.properties();
            end
        end
        T_tmp=1./accumarray(cg.cells.faces(:,1),1./HT_cg,[cg.faces.num,1]);
        T_cg(~internal)=T_tmp(~internal);
        if(any(~isfinite(HT_cg)) || any(~isfinite(T_cg )))
            warning('Some transmisibilities undesided set to zero')
            HT_cg(~isfinite(HT_cg))=0;
            T_cg(~isfinite(T_cg))=0;
        end
    otherwise
        error('This match method is not implemented')
end

%% fix some transmissibilities
if(~strcmp(opt.bc_method,'bc'))
    if(isfield(state,'facePressure'))
        hf_ext=~internal(cg.cells.faces(:,1));
        T_cg(cg.cells.faces(hf_ext))=HT_cg(hf_ext);
    else
        warning('facePressure not present set transmisibility of boundary faces to zero')
        T_cg(~internal)=0;
    end

    if(any(strcmp(opt.bc_method,{'wells','wells_simple'})))
        warning('Set boundary face trans to zero')
        hf_ext=~internal(cg.cells.faces(:,1));
        HT_cg(hf_ext)=0;
        T_cg(~internal)=0;
    end
    if(opt.fix_trans==true)
        % Set negative and zero transmissibility to the minimum positive
        % transmissibility.
        HT_min = min(HT_cg(HT_cg>0));
        HT_cg=max(HT_cg,HT_min);
        T_cg =max(HT_cg,HT_min);
        if(~isempty(cgwells))
            for j=1:numel(cgwells);
                cgwells(j).WI=max(cgwells(j).WI,0);
            end
        end
    end
    if(~isempty(opt.wells))
        if(~iscell(opt.wells))
            for j=1:numel(opt.wells);
                cgwells(j).type=opt.wells(j).type;
                cgwells(j).val=opt.wells(j).val;
            end
        else
            %cgwells=cgwells_str;
        end
    end
end

assert(all(isfinite(T_cg)))

report=struct('cg_state',upscaled,'state',state);
end
function cgwells  = makeCGWells(cg,wells)
    cgwells=[];
    for i=1:numel(wells);
       fcells= wells(i).cells;
       cgcells=cg.partition(fcells);
       %nperf=cellfun('prodofsize',[wells.cells]);
       %wno=rldecode(1:numel(wells),nperf);
       tab=sortrows([cgcells,fcells,[1:numel(fcells)]']);
       [cells,n]=rlencode(tab(:,1));
       fcellspos =cumsum([1;n]);
       pno=rldecode([1:numel(cells)]',diff(fcellspos));
       if(cg.griddim>2)
        hpos=accumarray(pno,cg.parent.cells.centroids(tab(:,2),3)...
              .*cg.parent.cells.volumes(tab(:,2)))./...
            accumarray(pno,cg.parent.cells.volumes(tab(:,2)));
       else
           hpos=0;
       end
       cgwells=[cgwells,struct('cells'   , cells,   ...
             'type'    , wells(i).type,             ...
             'val'     , wells(i).val,              ...
             'r'       , wells(i).r,                ...
             'dir'     , wells(i).dir,              ...
             'WI'      , nan(numel(cells),1),       ...
             'dZ'      , hpos-wells(i).refDepth,    ...
             'name'    , wells(i).name,             ...
             'compi'   , wells(i).compi,            ...
             'refDepth', wells(i).refDepth,         ...
             'sign'    , wells(i).sign,...
             'fcells'  , tab(:,2),...%not needed???
             'fcellspos',fcellspos,...
             'fperf',tab(:,3))];
   end
end
function upscaled = calculateUpscaledSolution(cg,state,cgwells)
    cg_volumes= accumarray(cg.partition,cg.parent.cells.volumes);
    %
    cfacesno  = rldecode([1:cg.faces.num]',diff(cg.faces.connPos));
    % calculate coarse face areas
    cg_areas   = accumarray(cfacesno,cg.parent.faces.areas(cg.faces.fconn));
    pressure   = accumarray(cg.partition,state.pressure.*cg.parent.cells.volumes)./cg_volumes;
    % calculate total flux over coarse face
    cfsign     = fineToCoarseSign(cg);
    flux = accumarray(cfacesno,state.flux(cg.faces.fconn).*cfsign);
    if(isfield(state,'facePressure'))
          facePressure = accumarray(cfacesno,state.facePressure(cg.faces.fconn).*cg.parent.faces.areas(cg.faces.fconn))./cg_areas;
    end
    % make upscaled well and well quantities
    if(~isempty(cgwells))
       wellSol=state.wellSol;
       for i=1:numel(cgwells);
           %% need more refinement
          pno=rldecode([1:numel(cgwells(i).cells)]',diff(cgwells(i).fcellspos));
          %wellSol=rmfield(wellSol,'flux');
          wellSol(i).flux=accumarray(pno,state.wellSol(i).flux(cgwells(i).fperf));
       end
    else
       wellSol=[];
    end

    upscaled=struct('pressure',pressure,...
       'flux',flux,...
       'facePressure',facePressure,...
       'wellSol',wellSol);
end


