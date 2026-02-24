function [g ceq] = Constraint(alpha, Par)
    centers = Par.obs(:,1:2);
    radius = Par.obs(:,3);
    P = bernstein_path(alpha,Par);
    g = obstacle_distance(centers,radius,P);
    ceq = [];
end