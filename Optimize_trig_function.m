% energy computation constants
W_0 = 14000;
D_CP = 100;
tau_W = 546 * exp(-0.01 * D_CP) + 316;

% model constants
num_iters = 400;
CP = calculate_CP(num_iters);
a_ini = 1140;
b_ini = [1500 3.14159265/30 0 0];
Total_Dist = 8000;

%% Find Optimal Power Function
w_bal_ode_a = @(a, b) w_bal_ode(num_iters, CP, W_0, a, b);
%a = Find_A(w_bal_ode_a, Total_Dist, b_ini)
%b = b_ini;

%[W_bal_ode, ext_time] = w_bal_ode(num_iters, CP, W_0, a_ini, b_ini);
%T = Finish_Time(Total_Dist, num_iters, a_ini, b_ini);

[a, b, T_min] = Optimize(w_bal_ode_a, Total_Dist);
%[a, b, T_min] = Find_B(w_bal_ode_a, Total_Dist);


%% plot results
to_plot = true;

if to_plot
    %b = 0.8
    %a = Find_A(w_bal_ode_a, Total_Dist, b);
    [W_bal_ode, ext_time] = w_bal_ode(num_iters, CP, W_0, a, b);
    T = Finish_Time(Total_Dist, num_iters, a, b);

    hold on
    plot(W_bal_ode, '-r');
    hold off

% plot power function and CP function
    figure;
     t = 1:num_iters;
    hold on
    plot(t, P(t, a, b));
    yline(CP);
    hold off

% plot velocity function and distance function
    figure;
    plot(t, velocity(t, a, b));
    
    figure;
    distance = zeros(1, num_iters);
    for i=1:num_iters
        distance(i) = trapz(velocity(1:i, a, b));
    end
    plot(distance(t));
end


%% Functions
% actual power out
function power = P(t, a, b)
    power = b(1) * sin(b(2) * t) + b(3) * cos(b(4) * t) + a;
end

% computes CP given duration of ride
function CP = calculate_CP(total_time)
    exponential = exp(-0.2 * total_time + 4);
    CP = 10000 * exponential / (1 + exponential) + 4000;
end

% returns expended energy at time u
function output = W_exp(u, CP, a, b)
    if P(u, a, b) <= CP
        output = 0;
    else
        output = P(u, a, b) - CP;
    end
end

function velocity = velocity(t, a, b)         %%%% depend on P(t)
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
    velocity = ((P(t, a, b) - constant)/(0.5 * drag_coefficient * frontal_area * rho + head_wind));
    velocity = nthroot(velocity, 3);
end

%Time to finish the whole race
function time = Finish_Time(totalDist, num_iters, a, b)
    rDist = @(t) totalDist - Dist(t, a, b);
    for i = 0:num_iters
        if rDist(i) < 0
            break
        end
    end
    time = i;
end

%Distant traveled after time t
function dist = Dist(t, a, b)
    dist = trapz(velocity(1:t, a, b));
end

%%%%%%%  ODE MODEL FOR W_BAL
function [W_bal_ode, exhaust_time] = w_bal_ode(num_iters, CP, W_0, a, b)
    % initialize distributions of energy balance (in kilojoules)
    W_bal_ode = zeros(1, num_iters);
    W_bal_ode(1) = W_0;
    exhaust_time = -1;
    
    % fill up energy balance vectors using ODE
    for j=2:num_iters    
        % ODE method
        if P(j, a, b) > CP
            W_bal_ode(j) = W_bal_ode(j-1) - (P(j, a, b) - CP);
        else
            W_bal_ode(j) = W_0 - (W_0 - W_bal_ode(j-1)) * exp(-(CP - P(j, a, b)) / W_0);
        end
        
        if W_bal_ode(j) <= 0 && W_bal_ode(j-1) > 0
            %fprintf("ODE exhaustion occurred at iteration %d.\n", j);
            if exhaust_time < 0
                exhaust_time = j;
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% optimization%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Find a such that all the energy is used up when finishing the race
function a= Find_A(w_bal_ode_a, total_dist, b)
    a = 4000;
    step_l = 20;
    step_s = 1;
    iter_num = 400;
    a_hist = zeros(1, iter_num);
    time = Finish_Time(total_dist, iter_num, a, b);
    [W_bal_ode, exhaust_time] = w_bal_ode_a(a, b);
    for i = 1:iter_num
        a_hist(i) = a;
        if exhaust_time < 0
            a = a + step_l;
        elseif abs(time - exhaust_time) <= 1 || (i > 2 && a_hist(i)== a_hist(i-2))
            break;
        elseif time < exhaust_time - 5
            a = a + step_l;
        elseif time < exhaust_time
            a = a + step_s;
        elseif time > exhaust_time + 5
            a = a - step_l;
            if i > 1
                [W_bal_ode, exhaust_time] = w_bal_ode_a(a_hist(i-1), b);
                if time < exhaust_time
                    a = a_hist(i-1) + step_s;
                end
            end
        elseif time > exhaust_time
            a = a - step_s;
        end
        time = Finish_Time(total_dist, 400, a, b);
        [W_bal_ode, exhaust_time] = w_bal_ode_a(a, b);
    end
end

%%%% Optimize Power Output
function [a, b, T_min] = Optimize(w_bal_ode_a, total_dist)
    b = [2000 3.14159265/30 0 0];
    a = Find_A(w_bal_ode_a, total_dist, b);
    rng(0, 'twister');
    step = [20 0.0002 5 0.0001] .* (rand(1, 4)/2 + 0.5);
    iter_num = 100;
    b_last = b - 1;
    T_last = 0;
    T = Finish_Time(total_dist, 400, a, b);
    T_change = 1;
    b_change = [1 1 1 1];
    b_hist = zeros(iter_num, 4);
    b_hist(1, :) = b;
    for i = 1:iter_num
        T_last = T;
        b = b - step .* T_change ./ b_change
        a = Find_A(w_bal_ode_a, total_dist, b);
        T = Finish_Time(total_dist, 400, a, b)
        T_change = T - T_last;
        b_hist(i + 1, :) = b;
        b_change = b - b_hist(i, :);
        
        if abs(T_change) < 0.0005 && i > 10
            b = b_last;
            break;
        elseif Norm(b_change) == 0
            break;
        end
    end
    T_min = T;
end

function norm = Norm(b)
    norm = sqrt(b(1)^2 + b(2)^2 + b(3)^2 + b(4)^2);
end

%For testing minimum time %%%%% TO BE FIXED
function [a, b, T_min] = Find_B(w_bal_ode_a, total_dist)
    a = 4000;
    b = [2200 3.14159265/60 10 3.14159265/80];
    b_min = b;
    step = 0.02;
    iter_num = 100;
    T_min = 1000000;
    for i = 1:iter_num
        a = Find_A(w_bal_ode_a, total_dist, b);
        T = Finish_Time(total_dist, 400, a, b);
        if T < T_min
            T_min = T;
            b_min = b;
        end
    end
    b = b_min;
end
