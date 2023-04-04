function [m, m_se, m_var] = calc_mean(sig, p, p_set)

% calc_mean : calculates mean (and se, and variance) of variable 'sig'
% after grouping by variable p

% p must be a col (e.g. coh or dir); p_set is also a col
% sig can be a matrix, but only if rows correspond to p (mean is columnwise)

% CF: force p and p_set to be column vectors
if isrow(p_set)
    p_set = p_set';
end
if isrow(p)
    p = p';
end
if isrow(sig)
    sig = sig';
end

m = nan(size(p_set,1),size(sig,2));
m_se = nan(size(m));
m_var = nan(size(m));
if iscell(p_set)
    for i = 1 : length(p_set)
        I = ismember(p,p_set{i});
        m(i,:) = nanmean(sig(I,:), 1);
        if sum(I)>2
            m_se(i,:) = nanstd(sig(I,:))./sqrt(sum(~isnan(sig(I,:)),1));
            m_var(i,:) = nanvar(sig(I,:));
        else
            m_se(i,:) = NaN;
            m_var(i,:) = NaN;
        end
    end
else
    for i = 1 : length(p_set)
        I = ismember(p,p_set(i,:),'rows'); 
        m(i,:) = nanmean(sig(I,:), 1);
        if sum(I)>2
            m_se(i,:) = nanstd(sig(I,:)) ./ sqrt(sum(~isnan(sig(I,:)),1));
            m_var(i,:) = nanvar(sig(I,:));
        else
            m_se(i,:) = NaN;
            m_var(i,:) = NaN;
        end
    end
end

return