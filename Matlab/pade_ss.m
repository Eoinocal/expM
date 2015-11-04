
function eA = pade_ss(A, q, j)

    A = A / 2^j;
    
    eA = pade_pq(A, q, q);
    
    for i = 1:j
        eA = eA*eA;
    end
    
end