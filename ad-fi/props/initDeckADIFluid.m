function fluid = initDeckADIFluid(deck, varargin)

reg = handleRegions(deck, varargin{:});
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
         'ISWL'  ,'ISWCR'  ,'ISWU', ...
         'ISGL'  ,'ISGCR'  ,'ISGU', ...
         'ISOWCR','ISOGCR'};
flag = ~any( strcmp(propNm , excpt) );
end