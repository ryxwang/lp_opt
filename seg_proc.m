function [bd_start, bd_end, pi_est] = seg_proc(M_tot, Y_tot, n, w, q)
%read in 
%fileID = fopen('GM12878_primary_10kb_chr11_5mb_9mb.txt','r');
%M_tot = fscanf(fileID,'%f',[400,400]);
%fclose(fileID);
%thresh = quantile(M_tot(triu(true(size(M_tot)),0)), 0.9);
%fileID = fopen('GM12878_ctcf_10kb_chr11_5mb_9mb.txt','r');
%Y_tot = fscanf(fileID,'%f');
%Y_tot = sort(Y_tot);
%fclose(fileID);
n_max = size(M_tot,1);

%n=200; %segment length
%w=150; %shift

%options=optimoptions('quadprog');options.Display='off';options.TolX=1e-5;

start1=0;
end1=min(start1+n, n_max);
inds = (start1+1):end1;
M=M_tot(inds,inds);
thresh = quantile(M(triu(true(size(M)),0)), q);
A = (M>thresh)*1;
A(logical(eye(size(A)))) = 0;
Y = Y_tot(Y_tot<=max(inds)&Y_tot>=min(inds));
Y=Y-start1;
Y1 = max(Y);
[pi_e,bd1,bd2,k_opt] = select_k(A,Y);
bd_start1 = bd1+start1; bd_end1 = bd2+start1; pi_est1 = pi_e;

cnt = 0;
uind_col = [];
while end1<n_max 
    
start2=start1+w;
end2=min(start2+n, n_max);
inds=(start2+1):end2;
M=M_tot(inds,inds);
thresh = quantile(M(triu(true(size(M)),0)), q);
A = (M>thresh)*1;
A(logical(eye(size(A)))) = 0;
Y = Y_tot(Y_tot<=max(inds)&Y_tot>=min(inds));
Y2 = min(Y);
Y=Y-start2;
[pi_e,bd1,bd2,k_opt] = select_k(A,Y);
bd_start2 = bd1+start2; bd_end2 = bd2+start2; pi_est2=pi_e;

overlap = zeros(length(bd_start1), length(bd_start2));
for i = 1:length(bd_start1)
    for j = 1:length(bd_start2)
        c = max(bd_start1(i), bd_start2(j));
        d = min(bd_end1(i), bd_end2(j));
        overlap(i,j) = max(0,d-c);
    end
end

%collect nonoverlapping domains
uind_row = find(sum(overlap,2)==0);
uind = intersect(uind_col,uind_row);
if cnt==0
    bd_start = bd_start1(uind_row);
    bd_end = bd_end1(uind_row);
    pi_est = pi_est1(uind_row);
end
uind_col = find(sum(overlap)==0);
bd_start = [bd_start; bd_start1(uind)];
bd_end = [bd_end; bd_end1(uind)];
pi_est = [pi_est; pi_est1(uind)];

ind=find(overlap>0);
for i=ind'
    [ind1,ind2] = ind2sub(size(overlap),i);
    a1 = bd_start1(ind1); a2 = bd_start2(ind2);  
    b1 = bd_end1(ind1); b2 = bd_end2(ind2);
    if abs(a2-Y2)==0 & a1<a2  %first domain in second block starts at (or close) to the first site 
        a2 = a1;
    end
    if abs(b1-Y1)==0 & b2>b1
        b1=b2;
    end
    inter = min(b1,b2)-max(a1,a2);
    union = max(b1,b2)-min(a1,a2);
    jac = inter/(union+1e-8);
    if (a1<=a2 && b2<=b1) || (a2<=a1 && b1<=b2) %if nested take subset
        bd_start = [bd_start; max(a1,a2)];
        bd_end = [bd_end; min(b1,b2)];
        pi_est = [pi_est; max(pi_est1(ind1), pi_est2(ind2))];
    elseif jac>0.8 %merge to get the intersection
        bd_start = [bd_start; max(a1,a2)];
        bd_end = [bd_end; min(b1,b2)];
        pi_est = [pi_est; max(pi_est1(ind1), pi_est2(ind2))];
    end
end

bd_start1 = bd_start2; bd_end1 = bd_end2; pi_est1 = pi_est2;
Y1 = start2+max(Y); start1=start2; end1=end2;
cnt = cnt+1
end

bd_start = [bd_start; bd_start2(uind_col)];
bd_end = [bd_end; bd_end2(uind_col)];
pi_est = [pi_est; pi_est2(uind_col)];