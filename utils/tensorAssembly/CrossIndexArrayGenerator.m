classdef CrossIndexArrayGenerator 
    
    properties
        tbl1        % IndexArray for the first table
        tbl2        % IndexArray for the second table
        mergefds    % Field names for merging
         
        replacefds1 % Possibility to change field names of first table
                    % (before setting up map)
        replacefds2 % Possibility to change field names of second table
                    % (before setting up map)
    end
   
    methods
        
        function gen = CrossIndexArrayGenerator(varargin)
            
            map.tbl1 = [];
            map.tbl2 = [];
            map.replacefds1 = [];
            map.replacefds2 = [];
            map.mergefds = [];
            
            map = merge_options(map, varargin{:}); 
            
            if isempty(map.mergefds)
                map.mergefds = {};
            end            
            
        end
        
        function tbl = eval(gen)
            
            mergefds  = gen.mergefds;
            tbl1 = gen.tbl1;
            tbl2 = gen.tbl2;
            
            if ~isempty(gen.replacefds1)
                 tbl1 = replacefield(tbl1, gen.replacefds1);
            end
            
            if ~isempty(gen.replacefds2)
                 tbl2 = replacefield(tbl2, gen.replacefds2);
            end
            
            tbl = crossIndexArray(tbl1, tbl2, mergefds);

        end
        
    end
   
end
