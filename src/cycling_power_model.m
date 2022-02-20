% energy computation constants
W_0 = 14000;
D_CP = 100;
tau_W = 546 * exp(-0.01 * D_CP) + 316;

% Tokyo Olympic track elevation grades (by half-kilometer) and total
% distance
elevations = [-0.32, -0.1, -0.3, -0.3, -0.5, -0.44, -0.16, -0.5, 0.3, -0.1, -0.3, 0, 0.3, 0.6, 0.3, 0.2, 0.3, 0.4, 0.3, 0.5, 0.5, 0.7, 0.32, -0.22, -0.5, -0.4, -0.3, -0.36, -0.26, -0.18, -0.3, -0.6, -0.4, 0, 0.5, 0.5, 0.4, 0.4, -0.2, -0.3, -0.1, 0.1, -0.3, 0.5, 0.1, 0.1, 0.16, -0.04, -0.32, -0.1, -0.3, -0.3, -0.5, -0.44, -0.16, -0.5, 0.3, -0.1, -0.3, 0, 0.3, 0.6, 0.3, 0.2, 0.3, 0.4, 0.3, 0.5, 0.5, 0.7, 0.32, -0.22, -0.5, -0.4, -0.3, -0.36, -0.26, -0.18, -0.3, -0.6, -0.4, 0, 0.5, 0.5, 0.4, 0.4, -0.2, -0.3, -0.1, 0.1, -0.3, 0.5, 0.1, 0.1, 0.16, -0.04];
total = 42;

% model constants
num_iters = 100;
CP = calculate_CP(num_iters);

% initialize distributions of energy balance (in kilojoules)
W_bal_int = zeros(1, num_iters);
W_bal_ode = zeros(1, num_iters);
W_bal_ode(1) = W_0;
W_bal_int(1) = W_0;

% sample space for power curve parameters
x = 4000:50:4200;
y = 0:50:400;
A = cartesian(x, y, y, y);

% optimize parameters in power curve
min_v = zeros(1, 4);
min_time = Inf;
for k=1:size(A, 1)
    distance = 0;
    v = A(k, :);
    
    finished = false;
    exhausted = false;
    time = Inf;
    for j=2:num_iters

        % numerical integration method
        sum = 0;
        for i=1:j
            sum = sum + exp(-(j-i)/tau_W) * W_exp(i, v, CP);
        end

        W_bal_int(j) = W_0 - sum;

        if W_bal_int(j) <= 0 && W_bal_int(j-1) > 0
%             fprintf("Integral exhaustion occurred at iteration %d.\n", j);
        end

        % ODE method
        if P(i, v) > CP
            W_bal_ode(j) = W_bal_ode(j-1) - (P(i, v) - CP);
        else
            W_bal_ode(j) = W_0 - (W_0 - W_bal_ode(j-1)) * exp(-(CP - P(i, v)) / W_0);
        end

        if W_bal_ode(j) <= 0 && W_bal_ode(j-1) > 0
%             fprintf("ODE exhaustion occurred at iteration %d.\n", j);
            exhausted = true;
        end
        
        d = distance / 2;
        if floor(distance / 2) <= 0
            elevation = 0;
        else
            elevation = elevations(floor(distance / 2));
        end
        
        distance = distance + trapz(velocity(i-1:i, v, elevation) ./ 60);
        if distance >= total
            finished = true;
            time = j;
            break;
        end
    end
    
    if finished && ~exhausted && time <= min_time
        min_v = v;
        min_time = time;
    end
        
end

min_v
min_time

% plot results
% hold on
% plot(W_bal_int, '*g');
% plot(W_bal_ode, '*r');
% hold off

% plot power function and CP function
% figure;
% t = 1:100;
% hold on
plot(t, P(t, min_v));
% yline(CP);
% hold off

% plot velocity function and distance function
% figure;
% plot(t, velocity(t));

% figure;
% distance = zeros(1, num_iters);
% for i=2:num_iters
%     distance(i) = distance(i-1) + trapz(velocity(i-1:i));
% end
% plot(distance(t));

% power function
function power = P(t, v)
    power = v(1) + v(2) .* sin(2*t) + v(3) .* sin(t) + v(4) .* sin(3*t);
end

% computes CP given duration of ride
function CP = calculate_CP(total_time)
    exponential = exp(-0.2 * total_time + 4);
    CP = 10000 * exponential / (1 + exponential) + 4000;
end

% returns expended energy at time u
function output = W_exp(u, v, CP)
    if P(u, v) <= CP
        output = 0;
    else
        output = P(u, v) - CP;
    end
end

function velocity = velocity(t, w, elevation)
    % velocity computation constants
    rider_weight_lb = 165;
    bike_weight_lb = 17;
    frontal_area = 5.4788;
    drag_coefficient = 0.63;
    loss_dt = 2;
    G = elevation;
    head_wind = 0;
    C_rr = 0.005;
    rho = 0.0765;
    
    c1 = (1 - loss_dt / 100)^(-1);
    constant = (9.8067 * (bike_weight_lb + rider_weight_lb));
    constant = constant * sin(atan(G/100)) + C_rr * cos(atan(G/100)) * c1;
    velocity = ((P(t, w) - constant)/(0.5 * drag_coefficient * frontal_area * rho + head_wind));
    velocity = nthroot(velocity, 3);
end

function C = cartesian(varargin)
    args = varargin;
    n = nargin;

    [F{1:n}] = ndgrid(args{:});

    for i=n:-1:1
        G(:,i) = F{i}(:);
    end
    
    C = unique(G , 'rows');
end
