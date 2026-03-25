clear; clc; close all;

%% Robot & Sim. Data
R = 0.021;      % m
D = 0.1047;     % m
omegaBar = 10;  % rad/s

Ts = 0.15;      % s
b = 0.15;       % Look-ahead distance of point B 
kf = 44;        % sim. Time

% Time Vector 
t = 0:Ts:kf;    % Array of discrete time steps
N = length(t);

% Pre-allocate variables for speed and ref.
d_norms = zeros(1, N);
xr = zeros(1,N);
yr = zeros(1,N);
theta_r = zeros(1,N);

%% Compute reference trajectory

for k = 1:N
    tk = t(k);
    
    % Reference Trajectory 
    xr(k) = 0.6*sin(tk/3.5); 
    yr(k) = 0.6*sin(tk/7);
    
    % 1st derivatives 
    x_dot = (0.6 / 3.5) * cos(tk / 3.5);
    y_dot = (0.6 / 7.0) * cos(tk / 7.0);
    
    % 2nd derivatives 
    x_ddot = -(0.6 / 3.5^2) * sin(tk / 3.5);
    y_ddot = -(0.6 / 7.0^2) * sin(tk / 7.0);
    
    % Unicycle Reference Velocities (Eq 11) 
    vr = sqrt(x_dot^2 + y_dot^2);
    
    % Prevent division by zero 
    if (x_dot^2 + y_dot^2) == 0
        wr = 0;
    else
        wr = (y_ddot * x_dot - x_ddot * y_dot) / (x_dot^2 + y_dot^2);
    end
    
    theta_r(k) = atan2(y_dot, x_dot);
    
    % Virtual Reference Input ur(k) (Eq 17) 
    T_FL_inv = [cos(theta_r(k)), -b * sin(theta_r(k)); 
                sin(theta_r(k)),  b * cos(theta_r(k))];
            
    ur = T_FL_inv * [vr; wr];
    
    % Disturbance d(k) (Eq 16) 
    B = eye(2) * Ts;
    d = -B * ur;
    
    % Euclidean norm (length) of the disturbance vector
    d_norms(k) = norm(d);
end

% maximum disturbance radius for Set D
radii.rD = max(d_norms);
fprintf('radius of the disturbance set D (rd): %.4f\n', radii.rD);

%% Plots
figure;
hold on;
title('Reference Trajectory')
plot(xr, yr, 'b-.', 'LineWidth', 1.5); 
xlabel('X Position (m)');
ylabel('Y Position (m)');
grid on;
xlim([-1 1]);
ylim([-1 1]);
axis square
hold off;

%% Offline Phase

% Calculate radii
radii.rU = (2 * omegaBar * R * b) / (sqrt(4 * b^2 + D^2));
fprintf('radius of set circ. approx. U [see (20)]: %.4f\n', radii.rU);

radii.rBU = radii.rU * Ts;
fprintf('radius of BU [see (20)]: %.4f\n', radii.rBU);

if radii.rD < radii.rBU
    disp('Assumption 1 has been validated!')
else
    disp(['Assumption 1 cannot be validated: controller may not have sufficient' ...
        'authority to overcome disturbance'])
end

% Initial.
q0 = [0.6; 0; pi];
z0 = q0(1:2) + [b * cos(q0(3)); b * sin(q0(3))];
zr0 = [xr(1) + b * cos(theta_r(1)); yr(1) + b * sin(theta_r(1))];
zTilde0 = z0 - zr0;

% Compute ROSC sets Ti
radii.T = zeros(1,500);
radii.T(1) = radii.rD;
i = 2;

while ~isCovered(zTilde0, radii.T(i-1))
    radii.T(i) = radii.T(i-1) - radii.rD + Ts * radii.rU;
    i = i + 1;
end

function boolean = isCovered(zTilde0, r)
    dist = norm(zTilde0,2);
    boolean = dist <= r;
end

radii.T = radii.T(1:i-1);

%% Online
