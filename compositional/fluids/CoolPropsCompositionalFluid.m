classdef CoolPropsCompositionalFluid < CompositionalFluid
    % Create MRST compositional fluid by using the open CoolProp library
    % (can be downloaded from http://www.coolprop.org/)
    properties
        
    end
    
    methods
        function fluid = CoolPropsCompositionalFluid(names)
            if ischar(names)
                names = {names};
            end
            
            ncomp = numel(names);
            [Tcrit, Pcrit, rhocrit, acc, molarMass, T_ref, P_ref] = deal(zeros(1, ncomp));
            
            validChoices = CoolPropsCompositionalFluid.getFluidList();
            ok = ismember(lower(names), lower(validChoices));
            
            if ~all(ok)
                s ='Unable to create fluid. The following names were not known to coolprops: ';
                msg = [s, sprintf('%s ', names{~ok})];
                error(msg);
            end
            
            for i = 1:numel(names)
                isF = strcmpi(names{i}, validChoices);
                n = validChoices{isF};
                Tcrit(i) = CoolProp.Props1SI(n, 'Tcrit');
                Pcrit(i) = CoolProp.Props1SI(n, 'Pcrit');
                rhocrit(i) = CoolProp.Props1SI(n, 'rhocrit');
                acc(i) = CoolProp.Props1SI(n, 'acentric');
                molarMass(i) = CoolProp.Props1SI(n, 'molarmass');
                
                % assuming IUPAC STP ??
                T_ref(i) = 273.15;
                P_ref(i) = 100000*Pascal;
            end
            assert(all(isfinite(Tcrit)));
            assert(all(isfinite(Pcrit)));
            assert(all(isfinite(rhocrit)));
            assert(all(isfinite(acc)));
            assert(all(isfinite(molarMass)));
            Vcrit = molarMass./rhocrit;
            fluid = fluid@CompositionalFluid(names, Tcrit, Pcrit, Vcrit, acc, molarMass, T_ref, P_ref);
        end
    end
    methods (Static)
        function varargout = getFluidList()
            names = num2str(CoolProp.get_global_param_string('FluidsList'));
            names = strsplit(names, ',');
            if nargout > 0
                varargout{1} = names;
            else
                disp('Possible fluid choices are:');
                fprintf('%s\n', names{:});
            end
        end 
    end
end