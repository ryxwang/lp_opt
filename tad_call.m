function [bd_start, bd_end, pval] = tad_call(hic_file, ctcf_file, bin_size, q)

%hic_file = './input/chr21_GM12878_example.txt';
%ctcf_file = './input/chr21_GM12878_ctcf_example.txt'; 
%bin_size = 10000;
%q=[0.9,0.5,0.5];
[data_sparse, start_ind] = read_data({hic_file}, bin_size);
data_sparse = data_sparse{1};
Y = dlmread(ctcf_file);

%%%%%level 1 TADs%%%%%%%%%
rst = [];
[bd_start, bd_end, pi_est] = seg_proc(data_sparse, Y, 300, 250, q(1));
rst = [rst; [bd_start, bd_end, pi_est]];

bd=rst;
%enrichment test
pval=zeros(size(bd,1),1);
for i=1:size(bd,1)
    pval(i) = enrich_test(data_sparse, bd(i,1:2), 1);
end
bd=[bd,pval];
bd=bd(bd(:,4)<0.05,:); %p-value <0.05

%sub levels
bd_tot_lp = bd; l=1;
%2 inner layers
while l<=2

bd_size = bd(:,2)-bd(:,1);

bd_start=[]; bd_end=[]; pi_est=[]; pval=[];
for i = 1:length(bd_size)
       r1 = bd(i,1); r2 = bd(i,2);
       M_tot = full(data_sparse(r1:r2, r1:r2));
       Y_tot = Y(Y<=r2 & Y>=r1)-r1+1;
       if length(Y_tot)>1
            A = binarize(M_tot,q(l+1));
            [pi_e,temp_start,temp_end,k_opt] = select_k(A,Y_tot);
            %enrichment test
                for j = 1:length(temp_start)
                    if ((temp_end(j)-temp_start(j)+1)==size(M_tot,1))
                        pval = [pval; 1];
                    else
                        pval = [pval; enrich_test(M_tot, [temp_start(j),temp_end(j)], 1)];            
                    end
                end
        
            bd_start = [bd_start(:); (temp_start(:) + r1 -1)];
            bd_end = [bd_end(:); (temp_end(:) + r1 -1)];
            pi_est = [pi_est; pi_e];
       end

end
temp = [bd_start, bd_end, pi_est, pval]; 
bd=temp(temp(:,4)<.05,:);
bd_tot_lp=[bd_tot_lp;bd];
l = l+1;
end

bd_start = bd_tot_lp(:,1)+start_ind-1; bd_end = bd_tot_lp(:,2)+start_ind-1;
pval=bd_tot_lp(:,4);
dlmwrite('./output/chr21_GM12878_tads.txt',[(bd_start-0.5)*bin_size, (bd_end-0.5)*bin_size], '\t');
end
