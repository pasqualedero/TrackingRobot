clear; clc; close all;

%% Robot & Sim. Data
R = 0.021;      % m
D = 0.1047;     % m
u1 = 10;  % rad/s
u2 = 7;
u3 = 5;

Ts = 0.15;      % s
b = 0.15;       % Look-ahead distance of point B 
kf = 44;        % sim. Time

% Time Vector 
t = 0:Ts:kf;    % Array of discrete time steps
N = length(t);

% Pre-allocate variables 
d_norms = zeros(1, N);

xr = zeros(1,N);
yr = zeros(1,N);
theta_r = zeros(1,N);

vr = zeros(1,N);
wr = zeros(1,N);

ur = zeros(2,N);

%% Compute reference trajectory
alpha = 0.06;
beta = 9;
dv = zeros(2,N);
for k = 1:N
    tk = t(k);
    
    % Reference Trajectory 
    xr(k) = alpha*sin(tk/(beta/2)); 
    yr(k) = alpha*sin(tk/beta);
    
    % 1st derivatives 
    x_dot = (alpha / (beta/2)) * cos(tk / (beta/2));
    y_dot = (alpha / beta) * cos(tk / beta);
    
    % 2nd derivatives 
    x_ddot = -(alpha / (beta/2)^2) * sin(tk / (beta/2));
    y_ddot = -(alpha / beta^2) * sin(tk / beta);
    
    % Unicycle Reference Velocities (Eq 11) 
    vr(k) = sqrt(x_dot^2 + y_dot^2);
    
    % Prevent division by zero 
    if (x_dot^2 + y_dot^2) == 0
        wr(k) = 0;
    else
        wr(k) = (y_ddot * x_dot - x_ddot * y_dot) / (x_dot^2 + y_dot^2);
    end
    
    theta_r(k) = atan2(y_dot, x_dot);
    
    % Virtual Reference Input ur(k) (Eq 17) 
    T_FL_inv = [cos(theta_r(k)), -b * sin(theta_r(k)); 
                sin(theta_r(k)),  b * cos(theta_r(k))];
            
    ur(:,k) = T_FL_inv * [vr(k); wr(k)];
    
    % Disturbance d(k) (Eq 16) 
    B = eye(2) * Ts;
    d = -B * ur(:,k);
    dv(:,k) = d;
    % Euclidean norm (length) of the disturbance vector
    d_norms(k) = norm(d);
end

% maximum disturbance radius for Set D
radii.rD = max(d_norms);
fprintf('radius of the disturbance set D (rd): %.4f\n', radii.rD);

figure;
hold on;
disturbances = plot(dv(1,:),dv(2,:));
Dset = ellipsoid(radii.rD^2 * eye(2));
plot(Dset);
grid on;
hold off


%% Offline Phase

% Calculate radii u1
radii.rU1 = (2 * u1 * R * b) / (sqrt(4 * b^2 + D^2));
fprintf('radius of set circ. approx. U [see (20)]: %.4f\n', radii.rU1);

radii.rBU1 = radii.rU1 * Ts;
fprintf('radius of BU [see (20)]: %.4f\n', radii.rBU1);

if radii.rD < radii.rBU1
    disp('Assumption 1 has been validated!')
else
    disp(['Assumption 1 cannot be validated: controller may not have sufficient' ...
        'authority to overcome disturbance'])
end

radii.rD2 = radii.rD;
% Calculate radii u2
radii.rU2 = (2 * u2 * R * b) / (sqrt(4 * b^2 + D^2));
fprintf('radius of set circ. approx. U [see (20)]: %.4f\n', radii.rU2);

radii.rBU2 = radii.rU2 * Ts;
fprintf('radius of BU [see (20)]: %.4f\n', radii.rBU2);

if radii.rD2 < radii.rBU2
    disp('Assumption 1 has been validated!')
else
    disp(['Assumption 1 cannot be validated: controller may not have sufficient' ...
        'authority to overcome disturbance'])
end

radii.rD3 = radii.rD2;
% Calculate radii u3
radii.rU3 = (2 * u3 * R * b) / (sqrt(4 * b^2 + D^2));
fprintf('radius of set circ. approx. U [see (20)]: %.4f\n', radii.rU3);

radii.rBU3 = radii.rU3 * Ts;
fprintf('radius of BU [see (20)]: %.4f\n', radii.rBU3);

if radii.rD3 < radii.rBU3
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

% Compute ROSC sets Ti for u1
radii.T1 = zeros(1,500);
radii.T1(1) = radii.rD;
i = 2;

while ~isCovered(zTilde0, radii.T1(i-1))
    radii.T1(i) = radii.T1(i-1) - radii.rD + Ts * radii.rU1;
    i = i + 1;
end

radii.T1 = radii.T1(1:i-1);
Nradii1 = length(radii.T1);

% Compute ROSC sets Ti for u2
radii.T2 = zeros(1,500);
radii.T2(1) = radii.rD2;
i = 2;

while ~isCovered(zTilde0, radii.T2(i-1))
    radii.T2(i) = radii.T2(i-1) - radii.rD2 + Ts * radii.rU2;
    i = i + 1;
end

radii.T2 = radii.T2(1:i-1);
Nradii2 = length(radii.T2);

% Compute ROSC sets Ti for u3
radii.T3 = zeros(1,500);
radii.T3(1) = radii.rD3;
i = 2;

while ~isCovered(zTilde0, radii.T3(i-1))
    radii.T3(i) = radii.T3(i-1) - radii.rD3 + Ts * radii.rU3;
    i = i + 1;
end

radii.T3 = radii.T3(1:i-1);
Nradii3 = length(radii.T3);


%% Online
X = zeros(3,N);
X(:,1) = q0;

Z = zeros(2,N);
Zr = zeros(2,N);
ZTilde = zeros(2,N);

U_opt = zeros(2,N);

Hd = [-1/u1   0        1/u1      0;
      0         -1/u1      0         1/u1]';
T = [R/2 R/2; R/D -R/D];
A = eye(2);

WR = zeros(1,N);
WL = zeros(1,N);

Imin = zeros(1,N);

% --- Add these before the loop ---
u_pre = [0; 0];          % Track the previous control action
R_weight = 0.5 * eye(2); % Weight for the control rate penalty

for k = 1:N
    % Compute zTilde, store
    theta_k = X(3,k);
    z_k = X(1:2,k) + [b * cos(theta_k); b * sin(theta_k)];
    zr_k = [xr(k) + b * cos(theta_r(k)); yr(k) + b * sin(theta_r(k))];
    zTilde_k = z_k - zr_k;
    
    Z(:,k) = z_k;
    Zr(:,k) = zr_k;
    ZTilde(:,k) = zTilde_k;

    % Compute U(theta)
    T_FL_theta_k = [cos(theta_k) sin(theta_k); -sin(theta_k)/b cos(theta_k)/b];
    H_theta_k = Hd * inv(T) * T_FL_theta_k;

    % Compute ur_k
    
    % Find i_min
    i_min = 0;
    for i = 1:Nradii1
        if isCovered(zTilde_k, radii.T1(i))
            i_min = i; 
            Imin(k) = i_min; 
            break;
        end
    end

    % Solve optimization
    if i_min > 1
        objFun = @(u) objective(u, A, B, zTilde_k, ur(:,k),u_pre,R_weight);
        conFun = @(u) constr(u, A, B, zTilde_k, ur(:,k), radii.T1(i_min-1));
        options = optimoptions('fmincon', 'Display', 'none');
        [u_opt, ~] = fmincon(objFun, [0; 0], ...
                          H_theta_k, ones(4,1), ... 
                          [], [], ... 
                          [], [], ... 
                          conFun, ...  
                          options);  
    else
        leftPart = ones(4,1) + H_theta_k * inv(B) * zTilde_k;
        options = optimoptions("quadprog","Display","off");
        [ur_k_hat, ~] = quadprog(2 * eye(2), -2 * ur(:,k), H_theta_k, leftPart, [], [],[],[],[],options);
        u_opt = -inv(B) * zTilde_k + ur_k_hat;
    end

    % Store u_opt
    U_opt(:,k) = u_opt;

    % Compute wR,wL
    command = inv(T) * T_FL_theta_k * u_opt;
    WR(k) = command(1);
    WL(k) = command(2);

    % Apply to robot
    velocities = T * command;
    vk = velocities(1);
    wk = velocities(2);

    X(:,k+1) = X(:,k) + [Ts * vk * cos(theta_k); Ts * vk * sin(theta_k); Ts * wk];
    
    u_pre = u_opt;
end

%% Plots
figure;
hold on;
title('Reference Trajectory vs Actual Trajectory')
referenceTraj = plot(xr, yr, 'b-.', 'LineWidth', 1.5);
actualTraj = plot(X(1,:),X(2,:), 'LineWidth', 1.5, 'Color', 'r');
xlabel('X Position (m)');
ylabel('Y Position (m)');
grid on;
xlim([-1 1]);
ylim([-1 1]);
axis square
legend([referenceTraj, actualTraj], {'Reference','Actual'});
hold off;

% ROSC sets
figure;
hold on;
title('ROSC sets $\mathcal{T}_i, i=1,\dots,N$', 'Interpreter','latex')
grid on;
theta = linspace(0,2*pi,360);
for i = 1:Nradii1
    r = radii.T1(i);
    x = r*cos(theta);
    y = r*sin(theta);
    if i == 1
        hRed = plot(x,y,'-','Color','r');
    else
        plot(x,y,'-','Color','b');
    end
    axis square;
end

for i = 1:Nradii2
    r = radii.T2(i);
    x = r*cos(theta);
    y = r*sin(theta);
    if i == 1
        hRed = plot(x,y,'-','Color','r');
    else
        plot(x,y,'-','Color','g');
    end
    axis square;
end

for i = 1:Nradii3
    r = radii.T3(i);
    x = r*cos(theta);
    y = r*sin(theta);
    if i == 1
        hRed = plot(x,y,'-','Color','r');
    else
        plot(x,y,'-','Color','m');
    end
    axis square;
end
hDot = plot(zTilde0(1),zTilde0(2),'LineStyle','none','Marker','hexagram','MarkerSize',8,'Color','black','MarkerFaceColor','black');
points = plot(ZTilde(1,:), ZTilde(2,:), 'Color','m','Marker','o');
legend( [hRed,hDot,points],{'Set $\mathcal{D}$','$\tilde{z}(0)$','$\tilde{z}(k)$'},'Interpreter','latex');
hold off

figure;
hold on;
title('$i(k)$ value over time', 'Interpreter','latex');
plot(1:N,Imin, 'LineWidth',1.5);
grid on;
axis tight
hold off

figure;
hold on;
subplot(3,3,[1 2 3])
hold on
title('Angular Velocities')
maxLine = yline(u1, 'LineStyle','--','Color','r');
plotWR = plot(1:N,WR,'LineWidth',1.5);
plotWL = plot(1:N,WL,'LineWidth',1.5);
legend([plotWL,plotWR,maxLine],{'$\omega_L$','$\omega_R$','$\bar{\omega}$'},'Interpreter','latex');
grid on;
hold off
subplot(3,3,[4 5 6])
hold on
title('Actual vs Reference Angle')
unwrapped_theta_r = unwrap(theta_r);
plotThetaRef = plot(1:N,unwrapped_theta_r,'LineStyle','--','Color','black','LineWidth',1.5);
plotThetaActual = plot(1:N, X(3,1:N), 'LineWidth', 1.5);
legend([plotThetaRef, plotThetaActual], {'$\theta_r$', '$\theta$'}, 'Interpreter', 'latex');
grid on;
hold off
subplot(3,3,[7 8 9])
hold on
title('Control Input')
u1r_plot = plot(1:N, ur(1,:), 'LineWidth', 1.5, 'LineStyle','--','Color','b');
u1_plot = plot(1:N, U_opt(1,:),'LineWidth', 1.5, 'Color','b');
u2r_plot = plot(1:N, ur(2,:), 'LineWidth', 1.5, 'LineStyle','--','Color','g');
u2_plot = plot(1:N, U_opt(2,:),'LineWidth', 1.5, 'Color','g');
legend([u1r_plot,u1_plot,u2r_plot,u2_plot],{'$u(1)_r$','$u(1)_{opt}$','$u(2)_r$','$u(2)_{opt}$'},'Interpreter','latex');
grid on;
hold off;
hold off;

%% Functions
function boolean = isCovered(zTilde0, r)
    dist = norm(zTilde0,2);
    boolean = dist <= r;
end

