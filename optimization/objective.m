function [f] = objective(u, A, B, z_tilde_k, ur_k)
    z_tilde_k_1 = A * z_tilde_k + B * u - B * ur_k;
    f = norm(z_tilde_k_1)^2 + 0.5 * norm(u)^2;
end