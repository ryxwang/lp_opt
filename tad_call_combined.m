function [bd_start, bd_end, pval] = tad_call_combined(hic_files, ctcf_files, bin_size, q, outfile)

%hic_file1 = './input/chr21_GM12878_example.txt';
%hic_file2 = './input/chr21_HUVEC_example.txt';
%ctcf_file1 = './input/chr21_GM12878_ctcf_example.txt';
%ctcf_file2 = './input/chr21_HUVEC_ctcf_example.txt';
%bin_size = 10000;
%q=[0.9,0.5,0.5];

%read in
%hic_files={hic_file1, hic_file2};
%ctcf_files={ctcf_file1, ctcf_file2};
[data_sparse, start_ind] = read_data(hic_files, bin_size);
L=length(hic_files);
for i = 1:L
    Y{i} = dlmread(ctcf_files{i});
end
%take intersection of CTCF sites
Y_mat=zeros(L,size(data_sparse,1));
for i = 1:L
    Y_mat(i, Y{i})=1;
end
Y_int=find(sum(Y_mat)==L)';

%%%%%level 1 TADs%%%%%%%%%
rst = [];
[bd_start, bd_end, pi_est] = seg_proc_combined(data_sparse, Y_int, 300, 250, q(1));
rst = [rst; [bd_start, bd_end, pi_est]];

bd=rst;
%enrichment test
pval=zeros(size(bd,1),L);
for i=1:size(bd,1)
    for j=1:L
        pval(i,j)=enrich_test(data_sparse{j}, bd(i,1:2), 1);
    end
end
bd=[bd,pval];
bd=bd(sum(pval<0.05,2)==L,:); %p-value <0.05

%sub levels
bd_tot_lp = bd; l=1;
%2 inner layers
while l<=2

bd_size = bd(:,2)-bd(:,1);

bd_start=[]; bd_end=[]; pi_est=[]; pval=[]; temp_pval=zeros(1,L);
for i = 1:length(bd_size)
       r1 = bd(i,1); r2 = bd(i,2);
       Y_tot = Y_int(Y_int<=r2 & Y_int>=r1)-r1+1;
       if length(Y_tot)>1
           A = [];  
           for k = 1:L
               M_tot{k} = full(data_sparse{k}(r1:r2, r1:r2));
               A(k,:,:) = binarize(M_tot{k},q(l+1));
           end
           [pi_e,temp_start,temp_end,k_opt] = select_k_combined(A,Y_tot);

           %enrichment test
           for j = 1:length(temp_start)
               if ((temp_end(j)-temp_start(j)+1)==size(M_tot{1},1))
                   pval = [pval; ones(1,L)];
               else
                   for k=1:L
                       temp_pval(k)=enrich_test(M_tot{k}, [temp_start(j),temp_end(j)], 1);
                   end
                   pval = [pval; temp_pval];
               end
           end            
           bd_start = [bd_start(:); (temp_start(:) + r1 -1)];
           bd_end = [bd_end(:); (temp_end(:) + r1 -1)];
           pi_est = [pi_est; pi_e];
       end

end
temp = [bd_start, bd_end, pi_est, pval];
bd=temp(sum(pval<0.05,2)==L,:); %p-value <0.05
bd_tot_lp=[bd_tot_lp;bd];
l = l+1;
end

bd_start = bd_tot_lp(:,1)+start_ind-1; bd_end = bd_tot_lp(:,2)+start_ind-1;
pval=bd_tot_lp(:,4:(4+L-1));
dlmwrite(outfile,[(bd_start-0.5)*bin_size, (bd_end-0.5)*bin_size]);
end
