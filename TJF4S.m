function [x, incr, I, iter] = TJF4S(U, k, h, E, tol, maxiter)
    U_j = U(:, 1);  
    U_jj = U(:, 2); 
    x0 = U_j(2:end-1); 
    
    r = 1;      
    lambda = -5;
    psi = 0;    
    iter = 0;
    incr = tol + 1;
    I = [];
    Func = @(u) funcionpr([U_j, [U_jj(1); u; U_jj(end)]], k, h, E);

    while incr > tol && iter < maxiter
        Fx = Func(x0);

        D = divdiff(Func, x0 + r * Fx, x0 - r * Fx);
        
        y = x0 - D \ Fx;
        Fy = Func(y);
        
        vk = (Fy' * Fy) / (Fx' * Fx);
        Kk = 1 / (1 + lambda * vk);
        pk = Kk * (1 + psi * vk);
        qk = 2 * Kk * vk;
        
        x = y - D \ (pk * Fy + qk * Fx);

        incr = norm(x - x0, inf)+norm(Func(x),inf);
        I = [I, incr];
        
        x0 = x;
        iter = iter + 1;
    end
end