%% Code for "Clarifying the Reliability Paradox: Poor Test-Retest Reliability 
% Attenuates Group Differences"
% Authors: Povilas Karvelis & Andreea Diaconescu (2025)
%
% This script generates synthetic data to analyze the impact of test–retest 
% reliability on observed effect sizes.

%% Set Random Seed for Reproducibility

rng(42);

%% Illustrative Figure Parameters (for the "Reliability Paradox" Figure)

e = 0.5;   % Standard deviation of measurement error
b = 0.5;   % Scaling factor for between–subject variance
w = 1;     % Condition difference (shift)
ss = 200;  % Sample size for illustrative figure

% Generate two independent "true" samples
x = randn([ss,1]);   % True latent variable sample
y = randn([ss,1]);   % Second true latent variable sample

% Generate observed measures for performance (x) and for a second measure (y)
xt = x*b + e*randn([ss,1]);        % Test sample for performance
xr = x*b + e*randn([ss,1]);        % Retest sample for performance    
xp = x*b + w + e*randn([ss,1]);      % Paired sample for performance (shifted)
xd = x*b + w + e*randn([ss,1]);      % Paired sample retest (not used further)

% %% Figure 1: Variance visualization
% figure('WindowStyle', 'docked')
% h = daboxplot([x * b, xt, xt + w], 'xtlabels', ...
%     {'\sigma_b','\sigma_b + \sigma_e', '\sigma_b + \sigma_e + \sigma_w'}, ...
%     'whiskers', 0, 'scatter', 1, 'scattersize', 25, 'scatteralpha', 0.6, ...
%     'withinlines', 1, 'outliers', 0);
% 
% set(gca, 'FontSize', 12)
% 
% %% Figure 2: Population variance vs raw group diffs
% figure('WindowStyle', 'docked')
% 
% % Identify groups based on threshold
% threshold = prctile(x, 50);
% %threshold = prctile(x, 90);
% 
% aid = x >= threshold;
% bid = x < threshold;
% 
% % Boxplot for group below threshold
% h = daboxplot([x(bid), x(bid) * 2, x(bid) * 2 + randn([sum(bid),1])] + 20, ...
%     'xtlabels', {'1', '2', '3'}, 'whiskers', 0, 'scatter', 1, 'scattersize', 25, ...
%     'scatteralpha', 0.6, 'withinlines', 1, 'jitter', 1, 'outliers', 0);
% 
% delete(h.bx); delete(h.md)
% 
% % Boxplot for group above threshold
% h = daboxplot([x(aid), x(aid) * 2, x(aid) * 2 + randn([sum(aid),1])] + 20, ...
%     'xtlabels', {'1', '2', '3'}, 'whiskers', 0, 'scatter', 1, 'scattersize', 25, ...
%     'scatteralpha', 0.6, 'withinlines', 1, 'jitter', 1, 'outliers', 0, ...
%     'scattercolors', {'m', 'w'});
% 
% delete(h.bx); delete(h.md)
% 
% ylabel('Mental health assessment')
% % ylabel('Trait assessment')
% xlabel('Population variance')
% xticklabels({'\sigma_b', '2\sigma_b', '2\sigma_b + \sigma_e'})
% set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'LineWidth', 2)
% ylim([13 26]); xlim([0.75 3.4])
% yticklabels('')

%% SIMULATION PARAMETERS

ss = 10000;          % Overall sample size for simulations
N = 20;              % Number of levels for variance parameters
vr = [0.3,2];        % Range for between–subject and error variances
bs = linspace(vr(1), vr(2), N);   % Scaling factors for between–subject variance
es = linspace(vr(1), vr(2), N);   % Scaling factors for error variance

w = 2;  % Constant shift applied to simulated measures

%% PERFORMANCE RELIABILITY SIMULATION
% In this section, we vary the reliability of the performance measure (x1t)
% while keeping the trait measure (tr) used for grouping fixed.

x1 = randn([ss,1]);    % Latent variable representing performance and traits
x2 = randn([ss,1]);    % Independent sample for paired comparisons

for i = 1:N
    for j = 1:N
        % Simulate the performance measure at test and retest:
        x1t = x1 * bs(i) + es(j)*randn([ss,1]);  % Test sample
        x1r = x1 * bs(i) + es(j)*randn([ss,1]);  % Retest sample    
        iccs(i,j) = ICC([x1t, x1r], 'A-1');      % Compute test–retest reliability (ICC)

        % Compute the one-sample effect size (shifted by constant w):
        dos(i,j) = cohensd(x1t + w, [], 'one-sample');   

        % Compute the paired-sample effect size:
        % (Using x2 to ensure near-zero correlation between repeated measures)
        x1p = x2 * bs(i) + w + es(j)*randn([ss,1]);  
        dps(i,j) = cohensd(x1p, x1t, 'paired-sample');

        % Generate a fixed trait measure for grouping (independent of the performance simulation):
        tr = x1 + std(x1)*randn([ss,1]);  

        % Compute the correlation between the performance measure and the fixed trait:
        r(i,j) = corr(x1t, tr);

        % Determine group splits based on the fixed trait (tr):
        xh = x1t(tr >= median(tr));      % High-trait subgroup (median split)
        xl = x1t(tr <  median(tr));       % Low-trait subgroup (median split)
        dms(i,j) = cohensd(xh, xl, 'two-sample');  % Compute median split effect size

        % Define alternative group splits ("patients" and "controls"):
        xpa = x1t(tr >= prctile(tr,85));   % "Patients" group (top 15%)
        xco = x1t(tr < prctile(tr,85));     % "Controls" group
        dss(i,j) = cohensd(xpa, xco, 'two-sample');  % Compute patient split effect size

        % Record the variance parameters for later plotting:
        B(i,j) = bs(i);
        E(i,j) = es(j);        
    end
end

%% TRAIT RELIABILITY SIMULATION
% In this section, the performance measure is held constant (with good reliability),
% while the reliability of the trait measure (used for grouping) is varied.

x1t_fixed = x1*1.0 + 1.0*randn([ss,1]);  % Fixed performance measure with high reliability

for i = 1:N
    for j = 1:N
        % Simulate the trait measure at test and retest:
        tr_test = x1 * bs(i) + es(j)*randn([ss,1]);
        tr_retest = x1 * bs(i) + es(j)*randn([ss,1]);
        iccs_traits(i,j) = ICC([tr_test, tr_retest], 'A-1');  % Compute trait reliability (ICC)

        % Simulate the trait measure for grouping:
        tr_sim = x1 * bs(i) + es(j)*randn([ss,1]);

        % Compute the correlation between the simulated trait and the fixed performance measure:
        r_traits(i,j) = corr(tr_sim, x1t_fixed);

        % Determine group splits based on the simulated trait:
        xh_traits = x1t_fixed(tr_sim >= median(tr_sim));
        xl_traits = x1t_fixed(tr_sim < median(tr_sim));
        dms_traits(i,j) = cohensd(xh_traits, xl_traits, 'two-sample');

        % Define alternative group splits based on the simulated trait:
        xpa_traits = x1t_fixed(tr_sim >= prctile(tr_sim,85));
        xco_traits = x1t_fixed(tr_sim < prctile(tr_sim,85));
        dss_traits(i,j) = cohensd(xpa_traits, xco_traits, 'two-sample');

        % Record the variance parameters for plotting:
        B_trait(i,j) = bs(i);
        E_trait(i,j) = es(j);
    end
end

%% RELIABILITY PARADOX FIGURE
% This figure reproduces the original 2×3 layout showing the "reliability paradox"
% (i.e., the relationship between test–retest reliability and observed effect sizes).

figure('WindowStyle','docked', 'Name', 'Reliability Paradox');

% Subplot 1: Test–retest reliability scatter plot
subplot(2,3,1)
scatter(xt, xr, 40, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.7);
hold on; set(gca, 'FontSize', 14);
axis([-3 3 -3 3]);
plot([-3 3], [-3 3], 'k--', 'LineWidth', 1.5);
h = lsline; set(h, "LineWidth", 2);
xlabel('Measurement at T1'); 
ylabel('Measurement at T2');
title('Test–retest reliability');

% Subplot 2: One–sample effects (using xp)
subplot(2,3,2)
h = daboxplot(xp, 'xtlabels', {'Condition 1'}, 'whiskers', 0, ...
    'scatter', 1, 'scattersize', 25, 'scatteralpha', 0.6, 'outliers', 0, ...
    'color', [0.9 0.9 0.9]);
title('One–sample effects');
set(gca, 'FontSize', 14);

% Subplot 3: Paired–sample effects (comparing xt and xp)
subplot(2,3,3)
h = daboxplot([xt, xp], 'xtlabels', {'Condition 1', 'Condition 2'}, 'whiskers', 0, ...
    'scatter', 1, 'scattersize', 25, 'scatteralpha', 0.6, 'withinlines', 1, 'outliers', 0, ...
    'color', [0.9 0.9 0.9]);
title('Paired–sample effects');
set(gca, 'FontSize', 14);

% Subplot 4: Surface plot of test–retest reliability (ICC)
subplot(2,3,4)
surf(B, E, iccs); set(gca, 'FontSize', 14);
xlabel('Between–subject variance'); 
ylabel('Error variance'); 
title('Test–retest reliability (ICC)');
colorbar; clim([0 1]);
xlim(vr); ylim(vr);
view(2);

% Subplot 5: Surface plot of one–sample Cohen''s d
subplot(2,3,5)
surf(B, E, dos); set(gca, 'FontSize', 14);
xlabel('Between–subject variance'); 
ylabel('Error variance'); 
title('One–sample Cohen''s d');
colorbar;
xlim(vr); ylim(vr);
view(2);

% Subplot 6: Surface plot of paired–sample Cohen''s d
subplot(2,3,6)
surf(B, E, dps); set(gca, 'FontSize', 14);
xlabel('Between–subject variance'); 
ylabel('Error variance');  
title('Paired–sample Cohen''s d');
colorbar;
xlim(vr); ylim(vr);
view(2);

%% GROUP DIFFERENCES VS. RELIABILITY (Combined)
% This figure presents three rows of plots in a 3×3 grid.
% The top two rows correspond to the performance simulation (illustrative and surface plots),
% and the bottom row shows the corresponding results for the trait simulation.

c = colororder;
c(1,:) = [0.7 0.7 0.7];

% Extract a subsample for illustrative plots (from the performance simulation)
trs = tr(1:200); 
x1s = x1t(1:200);    % Performance measure from the final iteration
xhs = xh(1:100); 
xls = xl(1:100);
xcos = xco(1:160); 
xpas = xpa(1:40);

% Create a new figure with a 3×3 grid layout
figure('WindowStyle','docked','Name','Group Differences vs. Reliability (Combined)');

% Top Row: Performance Simulation Illustrative Plots
subplot(3,3,1)
scatter(trs, x1s, 40, 'MarkerEdgeColor','w', 'MarkerFaceColor','k', 'MarkerFaceAlpha', 0.7);
hold on; lsline; set(gca, 'FontSize', 14);
xlabel('Traits/symptoms'); 
ylabel('Performance');
title('Traits/symptoms vs. performance');
ylim([-7 7]);
vline(median(tr), 'm--');
vline(prctile(trs,85), 'r--');

subplot(3,3,2)
gr = {xcos, xpas};
daboxplot(gr, 'xtlabels', {'Controls','Patients'}, 'whiskers', 0, 'outliers', 0, ...
    'scatter', 1, 'scattersize', 25, 'scatteralpha', 0.6, 'colors', c([1,2],:));
title('Patient split'); ylabel('Performance'); set(gca, 'FontSize', 14);

subplot(3,3,3)
gr = {xls, xhs};
daboxplot(gr, 'xtlabels', {'Lower traits','Higher traits'}, 'whiskers', 0, 'outliers', 0, ...
    'scatter', 1, 'scattersize', 25, 'scatteralpha', 0.6, 'colors', c([1,4],:));
title('Median split'); ylabel('Performance'); set(gca, 'FontSize', 14);

% Middle Row: Performance Simulation Surface Plots
subplot(3,3,4)
surf(B, E, r); set(gca, 'FontSize', 14);
xlabel('Performance between–subject var', 'FontSize', 15); 
ylabel('Performance error var', 'FontSize', 15);
title('Pearson''s r');
colorbar; clim([0.1 0.8]); xlim(vr); ylim(vr);
view(2);

subplot(3,3,5)
surf(B, E, dss); set(gca, 'FontSize', 14);
xlabel('Performance between–subject var', 'FontSize', 15); 
ylabel('Performance error var', 'FontSize', 15);
title('Patient split Cohen''s d*');
colorbar; clim([0 1.5]); xlim(vr); ylim(vr);
view(2);

subplot(3,3,6)
surf(B, E, dms); set(gca, 'FontSize', 14);
xlabel('Performance between–subject var', 'FontSize', 15); 
ylabel('Performance error var', 'FontSize', 15);
title('Median split Cohen''s d');
colorbar; clim([0 1.5]); xlim(vr); ylim(vr);
view(2);

% Bottom Row: Trait Simulation Surface Plots
subplot(3,3,7)
surf(B_trait, E_trait, r_traits); set(gca, 'FontSize', 14);
xlabel('Traits/symptoms between–subject var', 'FontSize', 15); 
ylabel('Traits/symptoms error var', 'FontSize', 15);
title('Pearson''s r');
colorbar; clim([0.1 0.8]); xlim(vr); ylim(vr);
view(2);

subplot(3,3,8)
surf(B_trait, E_trait, dss_traits); set(gca, 'FontSize', 14);
xlabel('Symptoms between–subject var', 'FontSize', 15); 
ylabel('Symptoms error var', 'FontSize', 15);
title('Patient split Cohen''s d*');
colorbar; clim([0 1.5]); xlim(vr); ylim(vr);
view(2);

subplot(3,3,9)
surf(B_trait, E_trait, dms_traits); set(gca, 'FontSize', 14);
xlabel('Traits between–subject var', 'FontSize', 15); 
ylabel('Traits error var', 'FontSize', 15);
title('Median split Cohen''s d');
colorbar; clim([0 1.5]); xlim(vr); ylim(vr);
view(2);

%% EFFECT SIZES AND P-VALUES AS A FUNCTION OF RELIABILITY
% This section examines how effect sizes and p-values change as a function 
% of reliability for different effects size metrics and statistical tests

ss = 1000000;           % Sample size for this analysis
x1 = randn([ss,1]);     % True sample for performance
x2 = randn([ss,1]);     % Independent sample for paired comparisons
b = 0.5;                % Scaling factor

N = 50;
vr = [0.01,3];
es = linspace(vr(1), vr(2), N);

for j = 1:2
    for i = 1:N
        % Simulate test and retest samples for performance:
        x1t = x1*b + es(i)*randn([ss,1]);  % Test sample
        x1r = x1*b + es(i)*randn([ss,1]);  % Retest sample    
        icce(i) = ICC([x1t, x1r], 'A-1');  % Compute reliability (ICC)
    
        % Generate a trait distribution (with different levels of true correlation)
        switch j
            case 1
                tr = x1 + 1.73*std(x1)*randn([ss,1]);
            case 2
                tr = x1 + 0.5*std(x1)*randn([ss,1]);
        end

        % Compute the correlation between the performance measure and the trait:
        rc(i,j) = corr(x1t, tr);  
        
        % Determine group splits based on the median of the trait:
        xh = x1t(tr >= median(tr));  % High-trait subgroup
        xl = x1t(tr < median(tr));   % Low-trait subgroup
        dmsc(i,j) = cohensd(xh, xl, 'two-sample');  % Compute effect size (Cohen's d)
    
        % Compute rank-biserial correlation (normalized z-value from ranksum test)
        [~,~,stats] = ranksum(xh, xl);
        urc(i,j) = stats.zval / sqrt(numel(x1t)); 
    end
end

% Plot the effect sizes as a function of reliability (ICC)
figure('WindowStyle','docked', 'Name', 'Effect Sizes and P-values vs. Reliability');

subplot(1,2,1)
plot(icce, rc(:,1) ./ max(rc(:,1)), 'LineWidth', 2); hold on;
plot(icce, dmsc(:,1) ./ max(dmsc(:,1)), 'LineWidth', 2);
plot(icce, urc(:,1) ./ max(urc(:,1)), 'LineWidth', 2);
plot(0:0.01:1, 1*sqrt(0:0.01:1), 'k--', 'LineWidth', 2);
set(gca, 'FontSize', 14);
xlabel('Reliability (ICC)'); 
ylabel('Normalized observed effect size'); 
legend({'Pearson''s r', 'Cohen''s d', 'Rank-biserial r_{rb}'}, 'FontSize', 14);
title('r_{true} = 0.5');

% Create an inset axes for additional data (r_true = 0.9)
axes('Position', [0.32, 0.24, 0.12, 0.3]);
plot(icce, rc(:,2) ./ max(rc(:,2)), 'LineWidth', 2); hold on;
plot(icce, dmsc(:,2) ./ max(dmsc(:,2)), 'LineWidth', 2);
plot(icce, urc(:,2) ./ max(urc(:,2)), 'LineWidth', 2);
plot(0:0.01:1, 1*sqrt(0:0.01:1), 'k--', 'LineWidth', 2);
set(gca, 'FontSize', 14, 'YTickLabel', [], 'XTickLabel', []);
title('r_{true} = 0.9');

%% Additional Analysis: p-values vs. Reliability (Small Sample)
ss = 60;               % Smaller sample size for p-value analysis
x1 = randn([ss,1]);    % True sample for performance
x2 = randn([ss,1]);    % Independent sample

b = 0.5;
N = 30;
vr = [0.01,1];
es = linspace(vr(1), vr(2), N);

for i = 1:N
    for j = 1:20000
        % Simulate test and retest samples:
        x1t = x1*b + es(i)*randn([ss,1]);  % Test sample
        x1r = x1*b + es(i)*randn([ss,1]);  % Retest sample    
        iccp(i,j) = ICC([x1t, x1r], 'A-1');  % Compute reliability
    
        % Generate a trait distribution:
        tr = x1 + 1.73*std(x1)*randn([ss,1]);
        [~, rp(i,j)] = corr(x1t, tr);  % Compute correlation between performance and trait
    
        % Determine group splits (median split):
        xh = x1t(tr >= median(tr));  % High-trait subgroup
        xl = x1t(tr < median(tr));   % Low-trait subgroup
        [~, dp(i,j)] = cohensd(xh, xl, 'two-sample');  % Compute effect size
        up(i,j) = ranksum(xh, xl);  % Compute ranksum test statistic (p-value proxy)
    end
end

% Plot p-values versus reliability
subplot(1,2,2)
plot(mean(iccp,2), mean(rp,2), 'LineWidth', 2); hold on;
plot(mean(iccp,2), mean(dp,2), 'LineWidth', 2); hold on;
plot(mean(iccp,2), mean(up,2), 'LineWidth', 2);
set(gca, 'FontSize', 14); ylim([0 0.15]);
xlabel('Reliability (ICC)'); 
ylabel('p-value'); 
title('r_{true} = 0.5; N = 60');