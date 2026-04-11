function P = bernstein_path(alpha, Par)
    % alpha: [n x 1] vector of optimization variables (deviations)
    n = length(alpha);
    A = Par.A;      
    B = Par.B;     
    dc = Par.dc;       
    t = (0:dc:1)';   
    
    % Geometry setup
    dir_vec = B - A;
    L = norm(dir_vec);
    n_vec = [-dir_vec(2); dir_vec(1)] / L; 
    
    % --- Generalized Bernstein Basis ---
    % Degree of the polynomial k = n + 1 (since alpha1...alphan are internal nodes)
    % The start (t=0) and end (t=1) points are fixed to A and B.
    k = n + 1; 
    deviation = zeros(size(t));
    
    for i = 1:n
        % Bernstein Basis function B_{i,k}(t) = comb(k,i) * t^i * (1-t)^(k-i)
        % We use i=1 to n to exclude the fixed endpoints at i=0 and i=k+1
        coeff = nchoosek(k, i);
        basis = coeff .* (t.^i) .* ((1-t).^(k-i));
        deviation = deviation + alpha(i) * basis;
    end
    
    baseline = A' + t* dir_vec';
    P = baseline + deviation * n_vec';
end
