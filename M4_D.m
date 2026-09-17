function [x, incr, I, iter] = M4_D(U, k, h, E, tol, maxiter)
    
    % U en el instante anterior.    
    U_j = U(:, 1);
    % U en el  mismo instante. Solo se conoce los valores extremos (cc).
    U_jj = U(:, 2);   
    
    x = U_j(2:end-1); 
    % Para la primera vez se necesita una x que no se conoce. Se añade una
    % pequeña perturbación sobre la primera que tenemos
    %x0 = x .*(rand(1,1)*0.05+0.95);
    x0 = x .* 1.5 + 1.0;
    nx = length(x);
    Id = eye(nx);
    iter = 0; incr = tol + 1; I = [];

    Func = @(u) funcionpr([U_j, [U_jj(1); u; U_jj(end)]], k, h, E);
    
    while incr > tol && iter < maxiter
        Fx = Func(x);
        
        % Gamma = -[x(k), x(k-1); F]^-1. Se calcula solo -[x(k), x(k-1); F] y luego
        % la inversa se calcula con la contrabarra.
        T_gamma = - divdiff(Func, x, x0);
        
        % w = x - Gamma * F(x)
        w = x + T_gamma \ Fx;
        
        % [w, x; F]
        T_wx = divdiff(Func, w, x);
        
        % y = x - [w, x; F]^-1 * F(x)
        y = x - T_wx \ Fx;
        
        % mu = I - [w, x; F]^-1 * [y, w; F]
        T_yw = divdiff(Func, y, w);
        mu = Id - T_wx \ T_yw;
        
        % H(mu) = I + mu + mu^2
        H_val = Id + mu + mu^2;
        
        % z = [y, x; F]^-1 * F(y)
        T_yx = divdiff(Func, y, x);
        Fy = Func(y);
        z = T_yx \ Fy;
                
        I = [I, incr];
        
        x0 = x;
        % x = y - H(mu) * z
        x = y - H_val * z;

        Fx=Func(x);
        incr = norm(x - x0, inf) + norm(Fx, inf);
        iter = iter + 1;
    end
end