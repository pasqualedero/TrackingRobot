function [f] = objective(u, A, B, z_tilde_k, ur_k, u_pre, R_weight)
    % Calculate predicted state
    z_tilde_k_1 = A * z_tilde_k + B * u - B * ur_k;  
    % Quadratic tracking error + Control rate penalty
    f = (z_tilde_k_1' * z_tilde_k_1) + (u - u_pre)' * R_weight * (u - u_pre);
end