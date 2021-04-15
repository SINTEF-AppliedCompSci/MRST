classdef IndexArray
% Class used to store multiple-indices.
%
%
% SYNOPSIS:
%   IndexArray(structtbl, varargin)
%
% DESCRIPTION: Class that is used to store and manipulate arrays of indices
% (multi-indices). Such array describes a subspace of a tensor space
%
% PARAMETERS:
%   structtbl - Structure from which the index array is constructed
%   varargin  - 
%
% RETURNS:
%   class instance
%
% EXAMPLE:
%   `computeVagTrans`
%
% SEE ALSO:
%   `TensorProd`, `TensorMap`, `SparseTensor`.
    
    properties
        
        % NxD array of multiple-indices. Each column corresponds to a given indexing.
        inds
        % Names of the indices given in a cell array of dimension D
        fdnames 
        
        tblname % not really implemented for the moment
        parents % not really implemented for the moment
        isvirtual % The multiple indices are not set up. We keep the fdnames. This
                  % option is used in an optimization setup when the local
                  % indexing is taking care of by the user.
        vnum      % store size of array if it is virtual
    end
   
    methods
        
        function tbl = IndexArray(structtbl, varargin)

        % The standard input for the constructor takes is a MATLAB struct object
        % structtbl. The field names of structtbl are used to setup the index
        % names tbl.fdnames. The values of each field provides an index vector
        % given in the indexing corresponding to the field name. These vectors
        % should have same length and they are used to set up tbl.inds.
        %
        % To set the IndexArray manually, set structtbl=[]
            opt = struct('fdnames'  , []   , ...
                         'isvirtual', false, ...
                         'tblname'  , []   , ...
                         'inds'     , [], ...
                         'num', []);
            opt = merge_options(opt, varargin{:}); 
            
            if ~isempty(structtbl)
                fdnames = fieldnames(structtbl)';
                nfdnames = numel(fdnames);
                num = size(structtbl.(fdnames{1}), 1);
                inds = zeros(num, nfdnames);
                for i = 1 : nfdnames
                    inds(:, i) = structtbl.(fdnames{i});
                end
                
                tbl.fdnames = fdnames;
                tbl.inds    = inds;
            else
                fdnames = opt.fdnames;
                inds    = opt.inds;
                tbl     = tbl.setup(fdnames, inds);
                tbl.parents = [];
            end

            if opt.isvirtual
                tbl.isvirtual = true;
                vnum = opt.num;
                assert(~isempty(vnum), 'In case of virtual table, the size of the index array should be given');
                tbl.vnum = vnum;
            end
            tbl.tblname = opt.tblname;    
        end

        function n = num(tbl)
            
        % Gives the number of multiple-indices (size of columns in inds) or the value stored in the table is virtual
            if tbl.isvirtual
                n = tbl.vnum;
            else
                n = size(tbl.inds, 1);
            end
        end
        
        function tbl = setup(tbl, fdnames, inds)
        % Direct setup of IndexArray with given fnames and inds.
            assert(numel(unique(fdnames)) == numel(fdnames), ['repeted field names']);

            tbl.fdnames = fdnames;
            tbl.inds = inds;
            
        end
        
        function tbl = addInd(tbl, fdname, ind)
        % Add index in IndexArray. The new index is decribed by its name fdname and the
        % index vector ind.
            
            fdnames = tbl.fdnames;
            inds = tbl.inds;
            tblind = strcmp(fdname, fdnames);
            assert(~any(tblind), ['index array contains already field with that ' ...
                                'name']);
            assert(size(inds, 1) == size(ind, 1), ['input size index does not ' ...
                                'match in size']);
            
            inds    = [inds, ind];
            fdnames = {fdnames{:}, fdname};
            
            tbl.fdnames = fdnames;
            tbl.inds = inds;            
             
        end
        
        function tbl = duplicateInd(tbl, fdcell)
        %   tbl    - IndexArray
        %   fdcell - Names of the index to duplicate with names of the duplicated
        %   indices. The syntax is {'name', {'dupname1', 'dupname2'}}
            
            oldfd = fdcell{1};
            fd1 = fdcell{2}{1};
            fd2 = fdcell{2}{2};

            fdnames = tbl.fdnames;
            inds = tbl.inds;
            
            [isok, fdind] = ismember(oldfd, fdnames);
            assert(isok, 'field name not recognized');
            
            nfds = numel(fdnames);
            fdinds = [(1 : fdind), fdind, ((fdind + 1) : nfds)];
            
            fdnames = fdnames(fdinds);
            fdnames{fdind}     = fd1;
            fdnames{fdind + 1} = fd2;
            inds = inds(:, fdinds);
            
            tbl.fdnames = fdnames;
            tbl.inds = inds;
            
        end
        
        function inds = get(tbl, fdname)
        % get index vector for the indexing given by fdname
            fdnames = tbl.fdnames;
            tblinds = strcmp(fdname, fdnames);
            assert(any(tblinds), 'index array field name not recognized');
            inds = tbl.inds(:, tblinds);
            
        end    
        
        function inds = gets(tbl, getfdnames)
        % get the index vectors for the indexings given by fdnames
            fdnames = tbl.fdnames;
            inds = tbl.inds;
            [isok, getinds] = ismember(getfdnames, fdnames);
            assert(all(isok), 'index array field name not recognized');
            inds = tbl.inds(:, getinds);
            
        end    
        
        function tbl = removeInd(tbl, rmfdnames)
        % remove the indices given by the names rmfdnames
            fdnames = tbl.fdnames;
            inds = tbl.inds;
            for i = 1 : numel(rmfdnames)
                fdname = rmfdnames{i};
                [isok, ind] = ismember(fdname, fdnames);
                assert(isok, 'field does not exist');
                fdnames(ind) = [];
                inds(:, ind) = [];
            end
            tbl.fdnames = fdnames;
            tbl.inds = inds;
        end    

        function tbl = addLocInd(tbl, locindname)
        % Add local indexing (given by 1 : N ) and gives it the name locindname
            fdnames = tbl.fdnames;
            tbl.fdnames = {fdnames{:}, locindname};
            
            locind = (1 : tbl.num)';
            tbl.inds = [tbl.inds, locind];
            
        end    

        function print(tbl, varargin)
        % Display IndexArray in terminal
            fprintf('%s ', tbl.fdnames{:});
            fprintf('\n');
            
            switch nargin 
              case 1
                rg = (1 : tbl.num)';
                tblinds = (1 : numel(tbl.fdnames))';
              case 2
                rg = varargin{1};
                tblinds = (1 : numel(tbl.fdnames))';
              case 3
                rg = varargin{1};
                tblinds = strcmp(varargin{2}, tbl.fdnames);
            end
            
            display(tbl.inds(rg, tblinds));
        end
        
    end
    
        
   
end

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
