function x = getQuadraturePoints(degree, dim)

    a = sqrt(3/5);
    b = sqrt(1/3);

    switch dim
        case 1
            
            if degree <= 1
            
                x = 0;
                
            elseif degree == 2

                x = [-b; b];

            elseif degree == 3

                x = [-a; 0; a];

            end
            
        case 2
            
            if degree <= 1
                
                x = [0,0];
            
            elseif degree == 2

                x = [-a -a;
                      0  a;
                      a  -a];

            elseif degree == 3

                x = [ 0  0;
                      0  b;
                     -a -a;
                      a -a;
                     -a  a;
                      a  a];

            end
            
    end
    
end