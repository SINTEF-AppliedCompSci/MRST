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
                if isempty(opt.inds)
                    inds = ones(0, numel(fdnames));
                else
                    inds = opt.inds;
                end
                tbl     = tbl.setup(fdnames, inds);
                tbl.parents = [];
            end

            if opt.isvirtual
                tbl.isvirtual = true;
                vnum = opt.num;
                assert(~isempty(vnum), 'In case of virtual table, the size of the index array should be given');
                tbl.vnum = vnum;
            else
                tbl.isvirtual = false;
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
            if size(ind, 1) == 1
                ind = ind*ones(size(inds, 1), 1);
            end
            
            assert(size(inds, 1) == size(ind, 1), ['input size index does not ' ...
                                'match in size']);
            
            inds    = [inds, ind];
            fdnames = {fdnames{:}, fdname};
            
            tbl.fdnames = fdnames;
            tbl.inds = inds;            
             
        end
        
        function tbl = duplicateInd(tbl, fdcell)
        %   tbl    - IndexArray
        %   fdcell - Names of the index to duplicate with names of the duplicated indices. The syntax is {'name', {'dupname1', 'dupname2'}}
            
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

        function tbl = set(tbl, fdname, inds)
        % replace the row for a given variable name, return error if fdname is not found
            assert(tbl.isvirtual == false, 'This function cannot be used for virtual tables');
            fdnames = tbl.fdnames;
            tblinds = strcmp(fdname, fdnames);
            assert(any(tblinds), 'index array field name not recognized');
            tbl.inds(:, tblinds) = inds;
        end
        
        function inds = get(tbl, fdname)
        % get index vector for the indexing given by fdname
            assert(tbl.isvirtual == false, 'This function cannot be used for virtual tables');
            fdnames = tbl.fdnames;
            tblinds = strcmp(fdname, fdnames);
            assert(any(tblinds), 'index array field name not recognized');
            inds = tbl.inds(:, tblinds);
            
        end    
        
        function inds = gets(tbl, getfdnames)
            
        % get the index vectors for the indexings given by fdnames
            assert(tbl.isvirtual == false, 'This function cannot be used for virtual tables');
            fdnames = tbl.fdnames;
            inds = tbl.inds;
            [isok, getinds] = ismember(getfdnames, fdnames);
            assert(all(isok), 'index array field name not recognized');
            inds = tbl.inds(:, getinds);
            
        end    

        function tbl = filterfds(tbl, fdnames)
        % create table with only the indices given by fdnames, by removing columns. Note that the function does not
        % check for uniqueness of the rows and should therefore be used with care.

            newnameinds = [];
            newnames = {};
            for ifd = 1 : numel(fdnames)
                fdname = fdnames{ifd};
                if iscell(fdname)
                    newnameinds(end + 1) = ifd;
                    newnames{end + 1} = fdname{2};
                    fdnames{ifd} =  fdname{1};
                end
            end
            [isok, fdinds] = ismember(fdnames, tbl.fdnames);
            assert(all(isok), 'field does not exist');

            if ~isempty(newnameinds)
                fdnames(newnameinds) = newnames;
            end
            
            tbl.fdnames = fdnames;
            tbl.inds = tbl.inds(:, fdinds);
            
        end
        
        function tbl = removefds(tbl, fdnames)
        % create table where the column given by fdnames are removed. Note that the function does not
        % check for uniqueness of the rows and should therefore be used with care.

            [isok, fdinds] = ismember(fdnames, tbl.fdnames);
            assert(all(isok), 'field does not exist');
            
            tbl.fdnames(fdinds) = [];
            tbl.inds(:, fdinds) = [];
            
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
            opt = struct('range', (1 : tbl.num)', ...
                         'fdnames', {tbl.fdnames});

            opt = merge_options(opt, varargin{:});

            [found, indcol] = ismember(opt.fdnames, tbl.fdnames);
            assert(all(found), 'some field names were not found');
            
            indrow = opt.range;
            
            try

                t = array2table(tbl.inds(indrow, indcol));
                t.Properties.VariableNames = tbl.fdnames(indcol);
                disp(t)
                
            catch
                
                fprintf('%s ', tbl.fdnames{indcol});
                fprintf('\n');
                display(tbl.inds(indrow, indcol));
                
            end
        end
        
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
