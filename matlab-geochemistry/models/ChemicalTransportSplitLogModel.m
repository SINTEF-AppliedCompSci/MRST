classdef ChemicalTransportSplitLogModel < ChemicalTransportLogModel


    methods

        function model = ChemicalTransportSplitLogModel(G, rock, fluid, chemicalLogModel, ...
                                                        varargin)

            model = model@ChemicalTransportLogModel(G, rock, fluid, chemicalLogModel, ...
                                                    varargin{:});

        end

        function [problem, state]=  getEquations(model, state0, state, dt, ...
                                                        drivingForces, ...
                                                        varargin)

            opt = struct('Verbose', mrstVerbose, ...
                         'reverseMode', false,...
                         'resOnly', false,...
                         'iteration', -1);  % Compatibility only
            opt = merge_options(opt, varargin{:});
            itnum = opt.iteration;


            % First we solve the chemistry
            if ((itnum >= 2) & (itnum <= 3))
                chemModel = model.chemicalModel;
                [newstate, failure, report] = ...
                    chemModel.solveChemicalState(state);
                if report.Converged
                    state = newstate;
                end
            end

            % Then we assemble the full system of equations
            [problem, state] = getEquations@ChemicalTransportLogModel(model, ...
                                                              state0, state, ...
                                                              dt, drivingForces, ...
                                                              varargin{:});

        end

    end
end
