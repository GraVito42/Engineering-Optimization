function [cost] = Objective(alpha, Par)
    P = bernstein_path(alpha,Par);
    dt = Par.dt;
    [total_cost, details] = bernstein_cost(P, dt);
    cost = details.length;
end
