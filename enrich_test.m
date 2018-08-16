%two-sample rank test
function p_val = enrich_test(M_tot, bd_ind, start_ind)
n_max = size(M_tot,1);
inner_ind = bd_ind-start_ind+1; inner_n = inner_ind(2)-inner_ind(1)+1;
n1 = min(inner_ind(1)-1, floor(inner_n/2));
n2 = min(n_max-inner_ind(2), floor(inner_n/2));
outer_ind = [inner_ind(1)-n1, inner_ind(2)+n2];
outer = M_tot(outer_ind(1):outer_ind(2), outer_ind(1):outer_ind(2));
inner = M_tot(inner_ind(1):inner_ind(2), inner_ind(1):inner_ind(2));
inner_sample = zeros(1,size(inner,1)-1);
for i=1:(size(inner,1)-1)
    inner_sample(i) = mean(diag(inner,i));
end
outer((n1+1):(n1+inner_n), (n1+1):(n1+inner_n)) = NaN;
outer_sample = zeros(1,size(inner,1)-1);
for i=1:(size(inner,1)-1)
    outer_sample(i) = nanmean(diag(outer,i));
end
p_val = ranksum(inner_sample, outer_sample, 'tail','right');