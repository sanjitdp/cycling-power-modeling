% energy computation constants
W_0 = 14000;
D_CP = 100;
tau_W = 546 * exp(-0.01 * D_CP) + 316;

% model constants
num_iters = 100;
CP = calculate_CP(num_iters);

% initialize distributions of energy balance (in kilojoules)
W_bal_int = zeros(1, num_iters);
W_bal_ode = zeros(1, num_iters);
W_bal_ode(1) = W_0;
W_bal_int(1) = W_0;

% fill up energy balance vectors using ODE and integral methods
for j=2:num_iters
    % numerical integration method
    sum = 0;
    for i=1:j
        sum = sum + exp(-(j-i)/tau_W) * W_exp(i, CP);
    end
    
    W_bal_int(j) = W_0 - sum;
    
    if W_bal_int(j) <= 0 && W_bal_int(j-1) > 0
        fprintf("Integral exhaustion occurred at iteration %d.\n", j);
    end
    
    % ODE method
    if P(i) > CP
        W_bal_ode(j) = W_bal_ode(j-1) - (P(i) - CP);
    else
        W_bal_ode(j) = W_0 - (W_0 - W_bal_ode(j-1)) * exp(-(CP - P(i)) / W_0);
    end
    
    if W_bal_ode(j) <= 0 && W_bal_ode(j-1) > 0
        fprintf("ODE exhaustion occurred at iteration %d.\n", j);
    end
end

% plot results
hold on
plot(W_bal_int, '*g');
plot(W_bal_ode, '*r');
hold off

% plot power function and CP function
figure;
t = 1:100;
hold on
plot(t, P(t));
yline(CP);
hold off

% plot velocity function and distance function
figure;
plot(t, velocity(t));

figure;
distance = zeros(1, num_iters);
for i=1:num_iters
    distance(i) = trapz(velocity(1:i));
end
plot(distance(t));

% power function
function power = P(t)
    power = 5000 - 5000 * mod(floor(t / 20), 2);
end

% computes CP given duration of ride
function CP = calculate_CP(total_time)
    exponential = exp(-0.2 * total_time + 4);
    CP = 10000 * exponential / (1 + exponential) + 4000;
end

% returns expended energy at time u
function output = W_exp(u, CP)
    if P(u) <= CP
        output = 0;
    else
        output = P(u) - CP;
    end
end

function velocity = velocity(t)
    % velocity computation constants
    rider_weight_lb = 165;
    bike_weight_lb = 17;
    frontal_area = 5.4788;
    drag_coefficient = 0.63;
    loss_dt = 2;
    G = 0;
    head_wind = 0;
    C_rr = 0.005;
    rho = 0.0765;
    
    c1 = (1 - loss_dt / 100)^(-1);
    constant = (9.8067 * (bike_weight_lb + rider_weight_lb));
    constant = constant * sin(atan(G/100)) + C_rr * cos(atan(G/100)) * c1;
    velocity = ((P(t) - constant)/(0.5 * drag_coefficient * frontal_area * rho + head_wind));
    velocity = nthroot(velocity, 3);
end
