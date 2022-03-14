classdef myEclipseTwoPhaseTransferFunction < TransferFunction
    %% A corrected version of the Eclipse two-phase transfer function.
    properties
        shape_factor_object
        shape_factor_name
    end

    methods

        function transferfunction = myEclipseTwoPhaseTransferFunction(shape_factor_name,block_dimension)

            transferfunction = transferfunction@TransferFunction();
            transferfunction.nphases = 2;

            if (nargin<1)
                %% No information about shape factor is provided so we set it to 1
                transferfunction.shape_factor_name = 'ConstantShapeFactor';
                block_dimension = [1,1,1];
                shape_factor_value = 1;
                shape_factor_handle = str2func(transferfunction.shape_factor_name);
                transferfunction.shape_factor_object = shape_factor_handle(block_dimension,shape_factor_value);
            else
                if (nargin<2)
                    msg = ['ERROR: If you provide the shape factor name, make sure to',...
                          'provide also the fracture spacing, [lx,ly,lz]'];
                    error(msg);
                else
                    transferfunction.shape_factor_name = shape_factor_name;
                    shape_factor_handle = str2func(transferfunction.shape_factor_name);
                    transferfunction.shape_factor_object = shape_factor_handle(block_dimension);
                end
            end


        end


        function [Talpha] = calculate_transfer(ktf,model,fracture_fields,matrix_fields)

            %% All calculate_transfer method should have this call. This is a "sanity check" that
            % ensures that the correct structures are being sent to calculate the transfer
            ktf.validate_fracture_matrix_structures(fracture_fields,matrix_fields);

            %% The varibles
            pom = matrix_fields.pom;
            swm = matrix_fields.swm;
            pO = fracture_fields.pof;
            sW = fracture_fields.swf;

            %% Oil Saturations
            %Fracture Oil saturation
            sO  = 1 - sW;

            %Matrix Oil saturation
            som = 1-swm;


            %% Residual/connate saturations
            try
                swrm = model.fluid_matrix.swr;
                swrf = model.fluid.swr;

                %residual oil in presence of water
                snrm = model.fluid_matrix.snr;
                snrf = model.fluid.snr;
            catch
                swrm = 0;
                swrf = 0;

                %residual oil in presence of water
                snrm = 0;
                snrf = 0;
            end

            %% Pressures
            pcOW = 0;
            pcOWm = 0;

            if isfield(model.fluid, 'pcOW') && ~isempty(sW)
                pcOW  = model.fluid.pcOW(sW);
            end

            if isfield(model.fluid_matrix, 'pcOW') && ~isempty(swm)
                pcOWm  = model.fluid_matrix.pcOW(swm);
            end

            pwm = pom - pcOWm;
            pW = pO - pcOW;

             %% Evaluate Rel Perms
            %Rel perms for the transfer
            [krW, krO] = model.evaluateRelPerm({sW, sO},'medium','fracture');
            [krWm, krOm] = model.evaluateRelPerm({swm, som},'medium','matrix');

            %% Additional Properties
            km = model.rock_matrix.perm(:,1);
            muwm = model.fluid_matrix.muW(pwm);
            muom = model.fluid_matrix.muO(pom);
            bWm = model.fluid.bW(pwm);
            bOm = model.fluid.bO(pom);

            %matrix fluid densities
            rhow = model.fluid_matrix.rhoWS;
            rhon= model.fluid_matrix.rhoOS;

            %Gravity
            g = norm(gravity);

            %% Shape Factor
            %lx,ly,lz are the inverse of the fracture spacing given
            % by the user i.e. the matrix block dimensions
            [sigma]=ktf.shape_factor_object.calculate_shape_factor(model);

            %% Fluid heights

            %vertical matrix block height
            lz=ktf.shape_factor_object.block_dimension(:,3);

            sWe = (sW-swrf) ./(1-snrf-swrf);
            swme = (swm-swrm)./(1-snrm-swrm);
%             sWe = min(max(sWe, 0), 1);
%             swme = min(max(swme, 0), 1);
            hwf = sWe .* lz;
            hwm = swme .* lz;
         
            hnf = lz-hwf;
            hnm = lz-hwm;

            %here the matrix densities have been used
            delta_rho = (rhow-rhon);

            % Introduce pseudo-pressures to account for the pressure
            % gradient due to gravity, see Eclipse Technical Description  
            psi_of = pO + delta_rho*g*hnf;
            psi_wf = pW + delta_rho*g*hwf;
            psi_om = pom + delta_rho*g*hnm;
            psi_wm = pwm + delta_rho*g*hwm;

            dpsiw = (value(psi_wm-psi_wf)<=0);
            dpsio = (value(psi_om-psi_of)<=0);

            % Upstream weighing of mobilities
            krwt = krW.*dpsiw + krWm.*(~dpsiw);            
            krot = krO.*dpsio + krOm.*(~dpsio);

            %% Compute Transfer
            %(units 1/T)
            tw=(sigma.*bWm.*krwt./muwm).*( pW-pwm + g.*delta_rho.*(hwf-hwm));
            to=(sigma.*bOm.*krot./muom).*( pO-pom + g.*delta_rho.*(hnf-hnm));

            %% Note that we return a 2x1 Transfer since our model is 2ph
            Talpha{1} = tw;
            Talpha{2} = to;
            
        end

        function [] = validate_fracture_matrix_structures(ktf,fracture_fields,matrix_fields)
            %% We use the superclass to validate the structures of matrix/fracture variables
            validate_fracture_matrix_structures@TransferFunction(ktf,fracture_fields,matrix_fields);
        end

    end


end

%{
Copyright 2022 Geological Survey of Denmark and Greenland (GEUS).

Author: Nikolai Andrianov, nia@geus.dk.

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
