function [q, j] = optimal_qj(norm_A, tol)

    q = -1;
    j = -1;
    min = 100;
    
    for Q = 0:100
        for J = 0:100            
            val = 8*(norm_A / 2^J)^(2*Q) ...
                * factorial(Q)^2 / (factorial(2*Q) * factorial(2*Q+1));
            
            if (val < tol && Q+J < min)
                
                q = Q;
                j = J;
                min = q+j;
                
                fprintf('New pair (%i, %i)\n', q, j);
            end            
        end
    end
    q
    j
end

