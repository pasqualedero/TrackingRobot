function [c,ceq] = constr(u, A, B, z_tilde_k, ur_k, r1)
    % Calculate predicted state
    z_tilde_k_1 = A * z_tilde_k + B * u - B * ur_k;
    
    % Purely quadratic constraint (z^T * z <= r^2). No epsilon needed!
    c = (z_tilde_k_1' * z_tilde_k_1) - (r1^2);  
    ceq = [];
end