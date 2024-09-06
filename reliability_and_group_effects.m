% The code for reproducing results in "Clarifying the reliability paradox:
% poor test-retest reliability attenuates group differences"
% Povilas Karvelis & Andreea Diaconescu
% 2024


% set random seed for reproducibility
rng(42)

% the mapping between err std of two distributions and ICC
% xr(1,:,:) = [x1 + 3.00*randn([ss,1]), x1 + 3.00*randn([ss,1])]; % ICC = 0.1
% xr(2,:,:) = [x1 + 2.00*randn([ss,1]), x1 + 2.00*randn([ss,1])]; % ICC = 0.2
% xr(3,:,:) = [x1 + 1.52*randn([ss,1]), x1 + 1.52*randn([ss,1])]; % ICC = 0.3
% xr(4,:,:) = [x1 + 1.22*randn([ss,1]), x1 + 1.22*randn([ss,1])]; % ICC = 0.4
% xr(5,:,:) = [x1 + 1.00*randn([ss,1]), x1 + 1.00*randn([ss,1])]; % ICC = 0.5
% xr(6,:,:) = [x1 + 0.82*randn([ss,1]), x1 + 0.82*randn([ss,1])]; % ICC = 0.6
% xr(7,:,:) = [x1 + 0.66*randn([ss,1]), x1 + 0.66*randn([ss,1])]; % ICC = 0.7
% xr(8,:,:) = [x1 + 0.50*randn([ss,1]), x1 + 0.50*randn([ss,1])]; % ICC = 0.8
% xr(9,:,:) = [x1 + 0.33*randn([ss,1]), x1 + 0.33*randn([ss,1])]; % ICC = 0.9
% xr(10,:,:) = [x1 + 0.1*randn([ss,1]) , x1 + 0.1*randn([ss,1])]; % ICC = 1.0

%% plot illustration figure

e = 0.5; % error std
b = 0.5; % between-subject std
w = 1;   % condition diff
g = 1;   % group diff

ss = 200;

x = randn([ss,1]);   % true sample 
y = randn([ss,1]);   % true sample #2

xt = x*b + e*randn([ss,1]);        % test sample
xr = x*b + e*randn([ss,1]);        % retest sample    
xp = x*b + w + e*randn([ss,1]);    % paired sample
xd = x*b + w + e*randn([ss,1]);    % paired sample retest

yt = y*b - g + e*randn([ss,1]);    % test sample 2
yr = y*b*3 - g + e*randn([ss,1]);  % retest sample 2 
yp = y*b*3 - g + e*randn([ss,1]);  % paired sample 2 

% figure('WindowStyle','docked')
% h = daboxplot([x*b,xt,xt + w],'xtlabels', {'\sigma_b','\sigma_b + \sigma_e', '\sigma_b + \sigma_e + \sigma_w'},'whiskers',0,...
%     'scatter',1,'scattersize',25,'scatteralpha',0.6,'withinlines',1,'outliers',0);
% %title('b + e + w') %ylim([-25 20]);
% set(gca,'FontSize',12)
% 
% figure('WindowStyle','docked')
% h = daviolinplot([x,x*2,x*3] + 20 ,'xtlabels', {'1','2', '3'},'whiskers',0,...
%     'scatter',1,'scattersize',25,'scatteralpha',0.6,'withinlines',1,...
%     'jitter',1,'outliers',0,'violin','half2','color',[0.7 0.7 0.7]);
% xlim([0.5 3.5])
% delete(h.bx); delete(h.md)
% ylabel('Mental health assessment')
% xlabel('Population variance')
% set(gca,'FontSize',14, 'FontWeight','bold','LineWidth',2)
% 
% 
% figure('WindowStyle','docked')
% h = daboxplot([x,x*2,x*3] + 20 ,'xtlabels', {'1','2', '3'},'whiskers',0,...
%     'scatter',1,'scattersize',25,'scatteralpha',0.6,'withinlines',1,...
%     'jitter',1);
% xlim([0.5 3.5])
% delete(h.bx); delete(h.md)
% ylabel('Mental health assessment')
% xlabel('Population variance')
% set(gca,'FontSize',14, 'FontWeight','bold','LineWidth',2)


%% simulate observed within and between group effects

ss = 10000;            % sample size
x1 = randn([ss,1]);    % original sample
x2 = randn([ss,1]);    % another independent sample

w = 2; w2 = 1;
g = 1;

N = 20;
vr = [0.3,2];
bs = linspace(vr(1),vr(2),N);
es = linspace(vr(1),vr(2),N);

for i = 1:N
    for j = 1:N

        % test-retest
        x1t = x1*bs(i) + es(j)*randn([ss,1]);        % test sample
        x1r = x1*bs(i) + es(j)*randn([ss,1]);        % retest sample    
        iccs(i,j) = ICC([x1t,x1r],'A-1');            % compute reliability

        % one-sample effects
        dos(i,j) = cohensd(x1t+w,[],'one-sample');   % compute effect size

        % paired sample effects
        %x1p = x1*bs(i) + w + es(j)*randn([ss,1]);   % assume some corr
        x1p = x2*bs(i) + w + es(j)*randn([ss,1]);    % assume 0 corr
        dps(i,j) = cohensd(x1p,x1t,'paired-sample'); % compute effect size

        % correlation
        tr = x1+std(x1)*randn([ss,1]);      % trait distribution
        r(i,j) = corr(x1t,tr);   % correlation betwen traits and performance

        % median split
        xh = x1t(tr>=median(tr)); % subgroup with high traits
        xl = x1t(tr<median(tr));  % subgroup with low traits
        dms(i,j) = cohensd(xh,xl,'two-sample'); % compute effect size

        % patient split
        % xpa = x1t(tr>=mean(tr)+std(tr)); % patients
        % xco = x1t(tr<mean(tr)+std(tr));  % controls
        xpa = x1t(tr>=prctile(tr,85)); % patients
        xco = x1t(tr<prctile(tr,85));  % controls
        dss(i,j) = cohensd(xh,xl,'two-sample'); % compute effect size

        % record vars for plotting
        B(i,j) = bs(i);
        E(i,j) = es(j);        

    end
end

% the reliability paradox
figure('WindowStyle','docked')
subplot(2,3,1)
scatter(xt,xr,40,'MarkerEdgeColor','w', 'MarkerFaceColor', 'k',...
        'MarkerFaceAlpha',0.7); hold on; set(gca,'FontSize',14)
axis([-3 3 -3 3])
plot([-3 3], [-3 3], 'k--', 'LineWidth', 1.5)
h = lsline; set(h,"LineWidth",2)
xlabel('Measurement at T1'); 
ylabel('Measurement at T2')
title('Test-retest reliability')

subplot(2,3,2)
h = daboxplot(xp,'xtlabels', {'Condition 1'},'whiskers',0,...
    'scatter',1,'scattersize',25,'scatteralpha',0.6,'outliers',0,...
    'color', [0.9 0.9 0.9]);
title('One-sample effects') %ylim([-25 20]);
set(gca,'FontSize',14)

subplot(2,3,3)
h = daboxplot([xt,xp],'xtlabels', {'Condition 1','Condition 2'},'whiskers',0,...
    'scatter',1,'scattersize',25,'scatteralpha',0.6,'withinlines',1,'outliers',0,...
    'color', [0.9 0.9 0.9]);
title('Paired-sample effects') %ylim([-25 20]);
set(gca,'FontSize',14)

subplot(2,3,4)
surf(B,E,iccs); set(gca,'FontSize',14)
xlabel('Between-subject variance'); 
ylabel('Error variance'); 
title('Test-retest reliability (ICC)')
colorbar; clim([0 1]);
xlim(vr); ylim(vr)
view(2)

subplot(2,3,5)
surf(B,E,dos); set(gca,'FontSize',14)
xlabel('Between-subject variance'); 
ylabel('Error variance'); 
title('One-sample Cohen''s d')
colorbar
xlim(vr); ylim(vr)
view(2)

subplot(2,3,6)
surf(B,E,dps); set(gca,'FontSize',14)
xlabel('Between-subject variance'); 
ylabel('Error variance');  
title('Paired-sample Cohen''s d')
colorbar
xlim(vr); ylim(vr)
view(2)


% group differences vs reliability
figure('WindowStyle','docked')

c = colororder;
c(1,:) = [0.7 0.7 0.7]; 

% subsample the data for the illustrative plots
trs = tr(1:200); x1s = x1(1:200);    % traits/symptoms and performance
xhs = xh(1:100); xls = xl(1:100);    % high and low traits
xcos = xco(1:160); xpas = xpa(1:40); % controls and patients

% illustrative plot: traits/symptoms vs performance
subplot(2,3,1)
scatter(trs,x1s,40,'MarkerEdgeColor','w', 'MarkerFaceColor', 'k',...
        'MarkerFaceAlpha',0.7); hold on; set(gca,'FontSize',14)
h = lsline; set(h,"LineWidth",2)
xlabel('Traits/symptoms'); 
ylabel('Performance')
vline(median(tr),'m--');
vline(prctile(trs,85),'r--')
title('Traits/symptoms vs performance')

% illustrative plot: patient split
subplot(2,3,2)
gr{1} = xcos; gr{2} = xpas;
daboxplot(gr,'xtlabels', {'Controls','Patients'},'whiskers',0,'outliers',0,...
    'scatter',1,'scattersize',25,'scatteralpha',0.6,'colors',c([1,2],:))
title('Controls vs patients')
ylabel('Performance')
set(gca,'FontSize',14)

% illustrative plot: median split
subplot(2,3,3)
gr{1} = xls; gr{2} = xhs;
daboxplot(gr,'xtlabels', {'Lower traits','Higher traits'},'whiskers',0,'outliers',0,...
    'scatter',1,'scattersize',25,'scatteralpha',0.6,'colors',c([1,4],:))
title('Lower vs higher traits')
ylabel('Performance')
set(gca,'FontSize',14)

% pearson's correlation for traits/symptoms vs performance
subplot(2,3,4)
set(gca,'FontSize',14)
surf(B,E,r); set(gca,'FontSize',14)
xlabel('Between-subject variance'); 
ylabel('Error variance'); 
title('Pearson''s r')
colorbar; clim([0.1 0.8]);
xlim(vr); ylim(vr)
view(2)

% cohen's d for controls vs patients diffs
subplot(2,3,5)
surf(B,E,dss); set(gca,'FontSize',14)
xlabel('Between-subject variance'); 
ylabel('Error variance'); 
title('Patient split Cohen''s d')
colorbar; clim([0 1.5]);
xlim(vr); ylim(vr)
view(2)

% cohen's d for low vs high traits
subplot(2,3,6)
surf(B,E,dms); set(gca,'FontSize',14)
xlabel('Between-subject variance'); 
ylabel('Error variance'); 
title('Median split Cohen''s d')
colorbar; clim([0 1.5]);
xlim(vr); ylim(vr)
view(2)

%% how effect sizes and p-values change as a function of reliability

ss = 1000000;          % sample size
x1 = randn([ss,1]);    % original sample
x2 = randn([ss,1]);    % another independent sample

b = 0.5;

N = 50;
vr = [0.01,3];
es = linspace(vr(1),vr(2),N);

for j = 1:2
    for i = 1:N
    
        % test-retest
        x1t = x1*b + es(i)*randn([ss,1]);        % test sample
        x1r = x1*b + es(i)*randn([ss,1]);        % retest sample    
        icce(i) = ICC([x1t,x1r],'A-1');        % compute reliability
    
        % trait distribution 
        switch j
            case 1
                tr = x1 + 1.73*std(x1)*randn([ss,1]);  % 
            case 2
                tr = x1 + 0.5*std(x1)*randn([ss,1]);  % 
        end

        % compute correlation betwen traits and performance
        rc(i,j)= corr(x1t,tr);  
    
        % median split
        xh = x1t(tr >= median(tr)); % subgroup with high traits
        xl = x1t(tr < median(tr));  % subgroup with low traits
        dmsc(i,j)= cohensd(xh,xl,'two-sample'); % compute effect size
    
        [~,~,stats] = ranksum(xh,xl);
        urc(i,j) = stats.zval/sqrt(numel(x1t));
         
    
    end
end


figure('WindowStyle','docked')

subplot(1,2,1)
plot(icce,rc(:,1)./max(rc(:,1)),'LineWidth',2); hold on
plot(icce,dmsc(:,1)./max(dmsc(:,1)),'LineWidth',2);
plot(icce,urc(:,1)./max(urc(:,1)),'LineWidth',2);
plot(0:0.01:1, 1*sqrt(0:0.01:1),'k--', 'LineWidth',2);
set(gca,'FontSize',14)
xlabel('Reliability (ICC)'); 
ylabel('Observed effect size'); 
legend({'Pearson''s r', 'Cohen''s d','Mann-Whitney''s r'},'FontSize',16)


% create smaller axes in top right, and plot on it
axes('Position',[.32 .24 .12 .3])
plot(icce,rc(:,2)./max(rc(:,2)),'LineWidth',2); hold on
plot(icce,dmsc(:,2)./max(dmsc(:,2)),'LineWidth',2);
plot(icce,urc(:,2)./max(urc(:,2)),'LineWidth',2);
plot(0:0.01:1, 1*sqrt(0:0.01:1),'k--', 'LineWidth',2);
set(gca,'FontSize',14,'YTickLabel',[],'XTickLabel',[])

% xlabel('Reliability (ICC)'); 
% ylabel('Observed effect size');
% title('r_{true} = 0.8')

ss = 60;               % sample size
x1 = randn([ss,1]);    % original sample
x2 = randn([ss,1]);    % another independent sample

b = 0.5;

N = 30;
vr = [0.01,1];
es = linspace(vr(1),vr(2),N);

for i = 1:N
    for j = 1:20000

        % test-retest
        x1t = x1*b + es(i)*randn([ss,1]);        % test sample
        x1r = x1*b + es(i)*randn([ss,1]);        % retest sample    
        iccp(i,j) = ICC([x1t,x1r],'A-1');        % compute reliability
    
        % correlation
        tr = x1 + 1.73*std(x1)*randn([ss,1]);      % trait distribution
        [~, rp(i,j)] = corr(x1t,tr);   % correlation betwen traits and performance
    
        % median split
        xh = x1t(tr>=median(tr)); % subgroup with high traits
        xl = x1t(tr<median(tr));  % subgroup with low traits
        [~, dp(i,j)] = cohensd(xh,xl,'two-sample'); % compute effect size

        up(i,j) = ranksum(xh,xl);
         
    end
end

% p-value vs reliability
subplot(1,2,2)
plot(mean(iccp,2),mean(rp,2),'LineWidth',2); hold on
plot(mean(iccp,2),mean(dp,2),'LineWidth',2); hold on
plot(mean(iccp,2),mean(up,2),'LineWidth',2);
set(gca,'FontSize',14)
xlabel('Reliability (ICC)'); 
ylabel('p-value'); 