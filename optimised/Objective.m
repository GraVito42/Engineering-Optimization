function [cost] = Objective(alpha, Par)
    alpha1 = alpha(1);
    alpha2 = alpha(2);
    P = bernstein_path(alpha1,alpha2,Par);
    dt = Par.dt;
    [total_cost, details] = bernstein_cost(P, dt);
    cost = details.length;
end
