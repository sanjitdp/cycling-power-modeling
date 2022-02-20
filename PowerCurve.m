%%%% Fitting Power Curve

rider_weight_lb = 165;
kg_to_lb = 2.20462;

%% Power Data

time = [5/60 1 5 60];

% [5sec 1min 5min 1h]   units: W/kg
All_Rounder = [19.85 9.55 5.74 4.71];
Sprinter = [21.03 9.09 4.81 3.91];
Pursuiter = [16.89 10.12 5.84 4.35];
Climber = [16.89 8.74 5.53 4.98];

% Data for CP,  units: W/lb
All_Rounder_w = All_Rounder  .* kg_to_lb;
Sprinter_w = Sprinter  .* kg_to_lb;
Pursuiter_w = Pursuiter  .* kg_to_lb;
Climber_w = Climber  .* kg_to_lb;


%% Fitting Power Curve

Dev = zeros(1, 4);
[All_Rounder_curve, Dev(1)] = fit_sig(time, All_Rounder, kg_to_lb);
[Sprinter_curve, Dev(2)] = fit_sig(time, Sprinter, kg_to_lb);
[Pursuiter_curve, Dev(3)] = fit_sig(time, Pursuiter, kg_to_lb);
[Climber_curve, Dev(4)] = fit_sig(time, Climber, kg_to_lb);

figure;
hold on
fplot(@(t) All_Rounder_curve(t), [0 60], '-r');
plot(time, All_Rounder_w, '*r');

fplot(@(t) Sprinter_curve(t), [0 60], '-b');
plot(time, Sprinter_w, '*b');

fplot(@(t) Pursuiter_curve(t), [0 60], '-g');
plot(time, Pursuiter_w, '*g');

fplot(@(t) Climber_curve(t), [0 60], '-k');
plot(time, Climber_w, '*k');
hold off


function [curve, dev] = fit_sig(time, data, scale)
    ft = fittype( @(a, b, c, d, x) a*exp(-b * x + c) ./ (1 + exp(-b*x + c)) + d, 'independent', {'x'}, 'dependent', 'y');
    [model, gof] = fit(time.', data.', ft, 'StartPoint', [-125, -1.4, 1.9, 130]);
    coef = coeffvalues(model)
    curve = @(t) coef(1) *exp(-coef(2) * t + coef(3)) ./ (1 + exp(-coef(2)*t + coef(3))) + coef(4);
    curve = @(t) curve(t) .* scale;
    dev = gof.rmse;
end
