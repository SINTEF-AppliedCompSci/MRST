function constit = shear_normal_stress(Nc, Nd, mu, lambda, phi)

%if Nd == 1
%    Voigt = 1; 
%elseif Nd == 2
%    Voigt = [1 0 0 0; 0 0 0 1; 0 .5 .5 0];
%elseif Nd ==3
%    Voigt = [1 0 0 0 0 0 0 0 0; 0 0 0 0 1 0 0 0 0; 0 0 0 0 0 0 0 0 1; 0 0 0 0 0 .5 0 .5 0; 0 0 .5 0 0 0 .5 0 0; 0 .5 0 .5 0 0 0 0 0];
%end
%VoigtI = Voigt>0;


%{ 
Copyright 2014 Jan Nordbotten jan.nordbotten@math.uib.no
%} 
alpha = 0; %1e-10; 
constit = cell(Nc,1);


for iter1 = 1:Nc;
    for iter3 = 1:Nd %iter3 indicates the normal vector of the stress surface
        for iter4 = 1:Nd %iter4 indicates the unit directions of stress
            %The constitutive matrix gives the directional stress due to a given
            %displacement gradient
            if (iter3==iter4)
                temp = lambda(iter1)*eye(Nd,Nd) + phi(iter1)*(ones(Nd)-eye(Nd));
            else
                temp = phi(iter1)*eye(Nd);
            end
            if 1 % Nonsymmetric formulation
                temp(iter3,iter4)=temp(iter3,iter4)+(1+alpha)*mu(iter1);
                temp(iter4,iter3)=temp(iter4,iter3)+(1-alpha)*mu(iter1);
            else % Imposes pointwise symmetry - makes MPFA local system underdetermined:
                temp(iter3,iter4)=temp(iter3,iter4)+mu(iter1);
                temp(iter4,iter3)=temp(iter4,iter3)+mu(iter1);
            end
            
            
            constit{iter1}(sub2ind_loc([Nd,Nd],iter3, iter4),:) = reshape(temp,1,Nd^2);
        end
    end
end


end


