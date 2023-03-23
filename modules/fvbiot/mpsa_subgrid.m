function out = mpsa_subgrid(G,stiffness,mode,varargin)
% Discretize linear elasticity by a multipoint stress approximation.
%
% This wrapper allows for discretizaiton of parts of the domain, as opposed
% to the full discretization computed by mpsa.m 
% So far, the only version implemented handles the case where the
% computation in memory bound (see comment in mpsa.m). A similar approach
% could be used if the permeability or geometry has changed in part of the
% domain, however, at the moment this has not been implemented.
%
% Parameters:
% G - mrst grid structure, see http://www.sintef.no/projectweb/mrst/
%     To get an overview of the data format, type 'help grid_structure'
% stiffness - stiffness matrix, as created by by shear_normal_stress.
% mode - basically a parameter describing the reason for calling this
%    function and not mpfa.m. At the moment, this parameter defaults to 1,
%    which is the memory bound case
%
% Optional parametrs, defined as keyword - value pairs:
% 'eta' - Location of continuity point on the half edges, defined
%   according to Aavatsmark 2002. Between 0 and 1. Default value is 0, for
%   simplices, the value should be 1/3 (will give symmetric method, see
%   Klausen, Radu Eigestad 2008).
% bc - boundary conditions, as defined by the MRST function addBC. If none
%    are provided, homogeneous Neumann conditions will be assigned.
% invertBlocks - method for inverting block diagonal systems. Should be
%    either 'matlab' or 'mex'. The former is pure matlab, which will be
%    slow for large problems. Mex is substantially faster, but requires
%    that matlab has access to a mex / c compiler.
% maxMemory - maximum memory assigned to storing the inverse gradient
%    operator in mpfa. Note that although this gives a substantial part of
%    the total memory consumption, in particular for large grids, maxMemory
%    does NOT control total memory usage.
%
%{
Copyright 2015-2016, University of Bergen.

This file is part of FVBiot.

FVBiot is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

FVBiot is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this file.  If not, see <http://www.gnu.org/licenses/>.
%} 

opt = struct('MaxMemory',100e6,...
             'invertBlocks','mex',...
             'eta',0, ...
             'bc', []);

opt = merge_options(opt,varargin{:});

switch mode
    case 1
        % Memory bound
        Nd = G.griddim;
        
        [cno,nno, ~, fno, subfno] = createSubcellMapping(G);
        [~, ind] = unique(subfno);
        hfi = bsxfun(@minus,repmat(subfno(ind) * Nd,1,Nd),fliplr(1:Nd)-1);
        fi = bsxfun(@minus,repmat(fno(ind) * Nd,1,Nd),fliplr(1:Nd)-1);
        hf2f = sparse(fi,hfi,1);

        blocksz = accumarray(nno(ind),1);
        neq = 2 * G.griddim * blocksz;
        
        bnds = cumsum(neq.^2);
        lower = 1;
        mem = opt.MaxMemory;
        threshold = mem;
        iter = 0;
        while lower < G.nodes.num
            upper = find(bnds < threshold,1,'last');
            if upper == lower
                error('You need to spend sufficient memory for at least one block')
            end
            
            nodes = lower : upper;
            d = mpsa(G,stiffness,nodes,'returnHalfEdgeValues',1,'invertBlocks',opt.invertBlocks,...
                     'eta', opt.eta, 'bc', opt.bc);
            if lower == 1
                out.stress = hf2f * d.stress;
                out.boundStress = hf2f * d.boundStress;
                out.gradP = d.gradP;   
                out.divD = d.divD;
                out.stabDelta = d.stabDelta;
            else
                out.stress = out.stress + hf2f * d.stress;
                out.boundStress = out.boundStress + hf2f * d.boundStress;
                out.gradP = out.gradP + d.gradP;    
                out.divD = out.divD + d.divD;
                out.stabDelta = out.stabDelta + d.stabDelta;
            end
            lower = upper + 1;
            threshold = threshold + mem;
            iter = iter + 1;
        end
        disp(['Finished mpsa reduced memory in ' num2str(iter) ' iterations'])
        out.div = vectorDivergence(G);
        out.A = out.div * out.stress;
    case 2
        error('Partial update of discretization not yet implemented')
        
        
        
    otherwise
        error(['Unknown mode ' mode])
end

