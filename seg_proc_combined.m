function [bd_start, bd_end, pi_est] = seg_proc_combined(M_tot, Y_tot, n, w, q)

n_max = size(M_tot{1},1);
L=length(M_tot);

%options=optimoptions('quadprog');options.Display='off';

start1=0;
end1=min(start1+n, n_max);
inds = (start1+1):end1;
A = [];
for i=1:L
    A(i,:,:)=binarize(M_tot{i}(inds,inds), q);
end
Y = Y_tot(Y_tot<=max(inds)&Y_tot>=min(inds));
Y=Y-start1;
Y1 = max(Y);
[pi_e,bd1,bd2,k_opt] = select_k_combined(A,Y);
bd_start1 = bd1(:)+start1; bd_end1 = bd2(:)+start1; pi_est1 = pi_e;

cnt = 0;
uind_col = [];
while end1<n_max 
    ss=1;length_y=0;
    while length_y<=1
    start2=start1+ss*w;
    end2=min(start2+n, n_max);
    inds=(start2+1):end2;
    Y = Y_tot(Y_tot<=max(inds)&Y_tot>=min(inds));
    Y2 = min(Y);
    Y=Y-start2; length_y = length(Y);
    ss = ss+1;
    end
A = [];
for r=1:L
    A(r,:,:)=binarize(M_tot{r}(inds,inds), q);
end
[pi_e,bd1,bd2,k_opt] = select_k_combined(A,Y);
bd_start2 = bd1(:)+start2; bd_end2 = bd2(:)+start2; pi_est2=pi_e;

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
if length(ind)>0
for i=ind'
    [ind1,ind2] = ind2sub(size(overlap),i);
    a1 = bd_start1(ind1); a2 = bd_start2(ind2);  
    b1 = bd_end1(ind1); b2 = bd_end2(ind2);
    if abs(a2-Y2)<=2  
        a2 = a1;
    end
    if abs(b1-Y1)<=2
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
end

bd_start1 = bd_start2; bd_end1 = bd_end2; pi_est1 = pi_est2;
Y1 = start2+max(Y); start1=start2; end1=end2;
cnt = cnt+1
end

if n_max<=n
    bd_start = bd_start1; 
    bd_end = bd_end1;
    pi_est = pi_est1;
else
    bd_start = [bd_start; bd_start2(uind_col)];
    bd_end = [bd_end; bd_end2(uind_col)];
    pi_est = [pi_est; pi_est2(uind_col)];
end

