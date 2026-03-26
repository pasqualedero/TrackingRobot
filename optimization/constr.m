function [c,ceq] = constr(u, A, B, z_tilde_k, ur_k, r1)
    epsilon = 1e-2;
    z_tilde_k_1 = A * z_tilde_k + B * u - B * ur_k;
    c = norm(z_tilde_k_1) - (r1 - epsilon);  
    ceq = [];
end