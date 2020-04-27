function v = propertynames(x)

  % disable warning in call to struct(x) below
  warn_id = 'Octave:classdef-to-struct';
  initial_warn_id_state = warning('query', warn_id).state;

  warning('off', warn_id);
  
  all_props = fieldnames(struct(x)); % includes public, protected and
                                     % private properties

  warning(initial_warn_id_state, warn_id);
  
  % keep only public properties
  v = {};
  for i = 1:numel(all_props)
    name = all_props{i};
    try
      tmp = x.(name);
      
      % this method was accessible, so keep it
      v = [v; {name}];
    catch
      
      % this method was inaccessible, i.e. not public.  Discard it
    end_try_catch
    
  end
  
end


