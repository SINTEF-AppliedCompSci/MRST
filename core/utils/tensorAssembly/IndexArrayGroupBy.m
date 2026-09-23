classdef IndexArrayGroupBy < handle

    properties (SetAccess = immutable)

        tbl      % input Index Array
        gfdnames % grouping field names
        ofdnames % other field names

        bytbl    % Index Array for the groups (field names gfdnames)
        num      % number of groups (equal to bytbl.num)

        startinds % list of index starts for the group in the indexing of tbl
        endinds   % list of index ends for the group in the indexing of tbl

        %% helpers
        ind_ofdnames % index of ofdnames
    end

    properties 

        i % current index in the group list

    end

    methods

        function iagb = IndexArrayGroupBy(tbl, gfdnames)
            
            iagb.gfdnames = gfdnames;
            iagb.tbl      = tbl.sort(gfdnames, 'keepAllFields', true);
            iagb.bytbl    = iagb.tbl.proj(gfdnames);
            iagb.num      = iagb.bytbl.num;
            
            ind = ismember(tbl.fdnames, gfdnames);
            iagb.ind_ofdnames = find(~ind);
            iagb.ofdnames = tbl.fdnames(~ind);

            inds = iagb.tbl.gets(iagb.gfdnames);
            [~, n] = rlencode(inds, 1);

            iagb.startinds = [1; cumsum(n(1 : end - 1)) + 1];
            iagb.endinds   = cumsum(n);

            iagb.i = 0;
            
        end

        function reset(iagb)
            iagb.i = 0;
        end
        
        function ind = startInd(iagb)
        % index start for the group in the indexing of tbl
            ind = iagb.startinds(iagb.i);
        end
        
        function ind = endInd(iagb)
        % index end for the group in the indexing of tbl
            ind = iagb.endinds(iagb.i);
        end

        function isok = next(iagb)

            if (iagb.i + 1) > iagb.num
                isok = false;
                return
            else
                iagb.i = iagb.i + 1;
                isok = true;
            end
            
        end

        function ctbl = currenttbl(iagb)
        % current group
            ctbl = iagb.bytbl;
            ctbl.inds = ctbl.inds(iagb.i, :);
            
        end
        
        function gtbl = grouptbl(iagb)
            
            gtbl = IndexArray([], ...
                              'fdnames', iagb.tbl.fdnames(iagb.ind_ofdnames), ...
                              'inds', iagb.tbl.inds(iagb.startInd() : iagb.endInd(), iagb.ind_ofdnames));
            
        end

        function subtbl = subtbl(iagb)
            
            gtbl = iagb.grouptbl();
            ctbl = iagb.currenttbl();
            
            subtbl = crossIndexArray(ctbl, gtbl, {}, 'optpureproduct', true);
            
        end

        function n = subgroupnums(iagb)

            map = TensorMap();
            map.fromTbl  = iagb.tbl;
            map.toTbl    = iagb.bytbl;
            map.mergefds = iagb.gfdnames;
            map = map.setup();

            n = map.eval(ones(iagb.tbl.num, 1));
            
        end
        
    end
    
end

