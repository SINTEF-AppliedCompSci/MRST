classdef StochasticConfigurations
    % Template class for holding the stochastic configurations relevant for
    % all ensemble members in an MRSTEnsemble
    
    properties
        num = inf % inf means that we can sample new ensemble members on the fly
        data
        
    end
    
    methods
        
        function stochasticConfigurations = StochasticConfigurations()
            %error('Template class not meant for direct use!');
        end
        
        function specificProblem = getProblem(stochasticConfigurations, baseProblem, seed)
            % getProblem based on the baseProblem, create a problem with 
            % configuration given by the seed, which is the seed for the
            % random generator or an index of data.
            error('Template class not meant for direct use!');
        end
    end
end
    
