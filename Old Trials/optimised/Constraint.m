function [g ceq] = Constraint(alpha, Par)
    alpha1 = alpha(1);
    alpha2 = alpha(2);
    centers = Par.obs(:,1:2);
    radius = Par.obs(:,3);
    P = bernstein_path(alpha1,alpha2,Par);
    g = obstacle_distance(centers,radius,P);
    ceq = [];
end