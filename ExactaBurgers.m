function u = ExactaBurgers(X, T, Re, a, b, ci)
    nu = 1 / Re;
    L = b - a;
    n_y = 15000; 
    y = linspace(a - L, b + L, n_y);
    if ischar(ci) || isstring(ci)
        ci = str2func(ci);
    end
    u_initial = ci(y);
      
    if iscolumn(u_initial)
        u_initial = u_initial'; 
    end
    
    integral_u = cumtrapz(y, u_initial);
    log_Theta = -(1 / (2 * nu)) * integral_u; 
    
    u = zeros(size(X));
    
    for j = 1:size(T, 2)
        for i = 1:size(X, 1)
            x = X(i, j);
            t = T(i, j);
            
            if t == 0
                u(i, j) = ci(x); %Condición Inicial

            else
                log_G = -((y - x).^2) ./ (4 * nu * t);
                
                log_total = log_Theta + log_G;
                
                % Estabilización: Restamos el máximo para evitar Overflow/Underflow en la exponencial.
                % El factor de escala (exp(-max)) se cancela en la división final.
                max_log = max(log_total);
                product_Theta_G = exp(log_total - max_log);
                
                numerator_integral = trapz(y, (y - x) .* product_Theta_G);
                denominator_integral = trapz(y, product_Theta_G);
                
                if denominator_integral == 0
                    u(i, j) = 0;
                else
                    u(i, j) = - numerator_integral / (t * denominator_integral);
                end
            end
        end
    end
end