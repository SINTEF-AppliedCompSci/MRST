classdef TensorConvert
% Convert a sparse matrix from `fromTbl` to `toTbl` into a vector with sparsity given by pivottbl.
    properties
        
        fromTbl  % Table for column in the matrix
        toTbl    % Table for row in the matrix
        pivottbl % Table for the sparsity of the tensor

        replacefdsToTbl   % Possibility to change field names of first table
                          % (before setting up conversion)
        replacefdsFromTbl % Possibility to change field names of second table
                          % (before setting up conversion)

        % Dispatching indexes
        indTo   % dispatching indexes for toTbl from pivottbl
        indFrom % dispatching indexes for fromTbl from pivottbl
        inds    % dispatching indexes for the conversion, obtained from indTo and indFrom
        
        issetup % Flag is set to true is product has been set up.
        
        chunksize = 100000; % Chunk size for the computation of the product
        verbose = false;        
    end
    
    methods
        
        function tconv = TensorConvert(varargin)
            opts = struct('toTbl'            , [], ...
                          'fromTbl'          , [], ...
                          'replacefdsToTbl'  , [], ...
                          'replacefdsFromTbl', []);
            
            tconv = merge_options(tconv, varargin{:}); 

            fds = {'replacefdsToTbl', 'replacefdsFromTbl'};

            for ifd = 1 : numel(fds)
                fd = fds{ifd};
                if isempty(tconv.(fd))
                    tconv.(fd) = {};
                end                
            end

            tconv.issetup = false;
            
        end

        function tconv = setup(tconv)

            toTbl    = tconv.toTbl;
            fromTbl  = tconv.fromTbl;
            pivottbl = tconv.pivottbl;
            
            if ~isempty(tconv.replacefdsToTbl)
                toTbl = replacefield(toTbl, tconv.replacefdsToTbl);
            end
            
            if ~isempty(tconv.replacefdsFromTbl)
                fromTbl = replacefield(fromTbl, tconv.replacefdsFromTbl);
            end
            
            fdsTo   = toTbl.fdnames;
            fdsFrom = fromTbl.fdnames;
            pfds    = pivottbl.fdnames;

            assert(all(ismember(fdsTo, pfds)), ['There exist fields in toTbl that do ' ...
                                               'not belong to fields of pivottbl']);
            assert(all(ismember(fdsFrom, pfds)), ['There exist fields in fromTbl that do ' ...
                                               'not belong to fields of pivottbl']);
            map = TensorMap();
            map.fromTbl  = toTbl;
            map.toTbl    = pivottbl;
            map.mergefds = fdsTo;

            tconv.indTo = map.getDispatchInd();

            map = TensorMap();
            map.fromTbl  = fromTbl;
            map.toTbl    = pivottbl;
            map.mergefds = fdsFrom;

            tconv.indFrom = map.getDispatchInd();

            tconv = tconv.setupFromInds();
            
        end

        function tconv = setupFromInds(tconv)

            inds = sub2ind([tconv.toTbl.num, tconv.fromTbl.num], tconv.indTo, tconv.indFrom);
            
            tconv.inds    = inds;
            tconv.issetup = true;
            
        end
        
        function m = convert(tconv, M)

            assert(tconv.issetup, 'please setup TensorConvert first');
            
            inds      = tconv.inds;
            chunksize = tconv.chunksize;
            
            if ~isempty(chunksize)

                pivottbl = tconv.pivottbl;
                
                vsize   = numel(inds);
                nchunks = ceil(vsize/chunksize);

                m = sparse(vsize, 1);

                if tconv.verbose
                    fprintf('number of chunks %d\n', nchunks);
                end
                
                for ichunk = 1 : nchunks

                    if ichunk < nchunks
                        ind = (1 + (ichunk - 1)*chunksize) : ichunk*chunksize;
                    else
                        ind = (1 + (ichunk - 1)*chunksize) : vsize;
                    end
                    
                    m = m + sparse(ind, ones(numel(ind), 1), M(inds(ind)), pivottbl.num, 1);
                    
                end
                
                
            else
                
               m = M(inds);
                
            end
            
        end
        
    end

    
end
