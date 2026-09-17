function H = divdiff_analitica(u, v, A, B, ua, ub)
    N = length(u);
    u = u(:); v = v(:);
    
    % Vectores desplazados
    v_next = [v(2:end); ub];
    v_prev = [ua; v(1:end-1)];

    diag_sup = -A + B * u(1:end-1);      % Posición j+1
    diag_inf = -A - B * u(2:end);        % Posición j-1 (El menos es por u_{i-1})
    diag_main = (1 + 2*A) + B * (v_next - v_prev);

    % Ensamblaje
H = spdiags([ [diag_inf; 0], diag_main, [0; diag_sup] ], -1:1, N, N); 
end