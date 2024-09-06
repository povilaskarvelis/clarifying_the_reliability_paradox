function [d, p] = cohensd(x,y,stype)

if strcmp(stype,'paired-sample')
    % d = abs(nanmean(x)-nanmean(y)) / std(x-y);
    [~,p,~,stats] = ttest(x,y);
    d = stats.tstat/sqrt(numel(x));

elseif strcmp(stype,'one-sample')

    [~,p,~,stats] = ttest(x);
    d = stats.tstat/sqrt(numel(x));

elseif strcmp(stype, 'two-sample')

    % std_p = sqrt(((numel(x)-1)*nanvar(x) +(numel(y)-1)*nanvar(y))/(numel(x) + numel(y) - 2));
    % 
    % d = abs(nanmean(x)-nanmean(y))/std_p;

    [~,p,~,stats] = ttest2(x,y,'vartype','unequal');
    d = stats.tstat*sqrt(1/numel(x) + 1/numel(y));

end

end

