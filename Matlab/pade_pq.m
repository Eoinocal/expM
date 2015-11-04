
function eA = pade_pq(A, p, q)

    eA = N(A)/D(A);

    function NpqA = N(A)
        
        coef_N = eye(size(A));
        NpqA = eye(size(A));
        
        for k = 1:p
%           NpqA = NpqA + (factorial(p+q-k)*factorial(p)) / ...
%                (factorial(p+q)*factorial(k)*factorial(p-k)) * A^k;

            coef_N = (coef(p-k+1)) / (coef(p+q-k+1) * coef(k)) * (coef_N*A);
            NpqA = NpqA + coef_N;
        end
    end

    function DpqA = D(A)
              
        coef_D = eye(size(A));
        DpqA = eye(size(A));  
        
        for k = 1:q      
%           DpqA = DpqA + (factorial(p+q-k)*factorial(q)) / ...
%                (factorial(p+q)*factorial(k)*factorial(q-k)) * (-A)^k;
         
            coef_D = (coef_D * coef(q-k+1)) / (coef(p+q-k+1) * coef(k)) * -A;
            DpqA = DpqA + coef_D;
        end
    end
    
    function x = coef(xx)
        if xx == 0
            x = 0;
        else
            x = xx;
        end
    end
end
