function TDD = divdiff(fun, X, Y)
    m = length(X); 
    TDD = zeros(m, m);
    tol = eps * 100;

    for j = 1 : m
        diff_val = Y(j) - X(j);
        
        if abs(diff_val) < tol
            TDD(:, j) = 0;
            TDD(j, j) = 1;
        else
            v1 = X; v1(1:j) = Y(1:j);
            v2 = X; v2(1:j-1) = Y(1:j-1);
            
            [Fv1, ~] = feval(fun, v1);
            [Fv2, ~] = feval(fun, v2);
            
            TDD(:, j) = (Fv1(:) - Fv2(:)) / diff_val;
        end
    end
end