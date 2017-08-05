raw_file = '/Volumes/G-DRIVE mobile USB/K562/K562/10kb_resolution_intrachromosomal/chr21/MAPQGE30/chr21_10kb.RAWobserved';
norm_file = '/Volumes/G-DRIVE mobile USB/K562/K562/10kb_resolution_intrachromosomal/chr21/MAPQGE30/chr21_10kb.KRnorm';
norm_sparse = read_matrix(raw_file, norm_file);

seg_ind = dlmread('/Users/rachelwang/Documents/data/K562_10kb_chr21_seg_ind.txt');
ctcf_ind = dlmread('/Users/rachelwang/Documents/data/K562_ctcf_chr21_relaxed.txt');
ctcf_ind = sort(ctcf_ind);

%longest segment
%seg_len = seg_ind(:,2)-seg_ind(:,1);
%temp = find(seg_len==max(seg_len));
%r1 = seg_ind(temp,1); r2 = seg_ind(temp, 2);
%M_tot = full(norm_sparse(r1:r2, r1:r2));
%thresh = quantile(M_tot(triu(true(size(M_tot)),0)), 0.9);

rst = [];
for i=3:4
    r1 = seg_ind(i,1); r2 = seg_ind(i, 2);
    M_tot = full(norm_sparse(r1:r2, r1:r2));
    Y_tot = ctcf_ind(ctcf_ind<=r2 & ctcf_ind>=r1)-r1+1;
    if length(Y_tot)>0
        [bd_start, bd_end, pi_est] = seg_proc(M_tot, Y_tot, 300, 250, 0.9);
        bd_start = bd_start+r1-1; bd_end = bd_end+r1-1;
        rst = [rst; [bd_start, bd_end, pi_est]];
    end
end

bd1=rst;
%enrichment test
pval=zeros(size(bd1,1),1);
for i=1:size(bd1,1)
    pval(i) = enrich_test(full(norm_sparse), bd1(i,1:2), 1);
end
%[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pval,0.1,'pdep','no'); %FDR
bd1=[bd1,pval];
bd1=bd1(bd1(:,4)<0.05,:);
%bd1=bd1(h,:);
%dlmwrite('/Users/rachelwang/Documents/MATLAB/chr21_results/bd1_lp_chr21_q9_q5_q5_relaxed.txt', bd1);