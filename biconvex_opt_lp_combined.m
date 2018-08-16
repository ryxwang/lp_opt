%find common TADs given binary adjacency matrices A, common CTCF vector Y, upper bound k
%on the number of TADs
function [pi_mat,x,pi_est,bd_start,bd_end]=biconvex_opt_lp_combined(A,Y,k)

n=size(A,2);
L = size(A,1);
A_sub = zeros(length(Y), length(Y));
n_sub = zeros(length(Y), length(Y));

alpha_hat = [];
for l=1:L
    for a=1:(length(Y)-1)
        for b=(a+1):length(Y)
            A_sub(a,b) = sum(sum(A(l, Y(a):Y(b), Y(a):Y(b))));
            n_sub(a,b) = (Y(b)-Y(a))*(Y(b)-Y(a)+1);
        end
    end
    A_sub = A_sub(find(~tril(ones(size(A_sub)))));
    n_sub = n_sub(find(~tril(ones(size(n_sub)))));

    A_sub_arr(l,:) = A_sub;
    n_sub_arr(l,:) = n_sub;
    alpha_hat(l,:) = A_sub./(n_sub+1e-8);
end
m = length(A_sub);
alpha_hat(alpha_hat==1) = 1-1e-8;
alpha_hat(alpha_hat==0) = 1e-8;

lin_ind = zeros(length(Y), length(Y));
lin_ind(logical(triu(ones(length(Y)),1))) = 1:m;
B = zeros(length(Y), m);
for i=1:length(Y)
    B(i, nonzeros(lin_ind(1:i, i:length(Y)))) = 1;
end

%updating Pi
%initial beta
beta_hat_new = [];
for l=1:L
    beta_hat_new(l,:) = sum(sum(A(l,:,:))) /(n*(n-1));
end
delta = 1;
iter = 0;
eps=1e-8;
while delta > 1e-5
    iter = iter+1;
    beta_hat = beta_hat_new;
    
    beta_hat_arr = repmat(beta_hat,1,m);
    t = 0.5*(A_sub_arr.*log((eps+alpha_hat)./(1-alpha_hat+eps).*(eps+1-beta_hat_arr)./(eps+beta_hat_arr)) + n_sub_arr.*log((1-alpha_hat+eps)./(eps+1-beta_hat_arr)));
    t2 = sum(t,1);
    

    A1(1:size(B,1),:)=B;
    A1(size(B,1)+1,:)=ones(1,length(t2));
    b(1:size(B,1))=1;
    b(size(B,1)+1)=k;
    x = linprog(-t2,A1,b,[],[],zeros(1,length(t2)),1+zeros(1,length(t2)));
    
    for l=1:L
        beta_hat_new(l,:) = (sum(sum(A(l,:,:))) - sum(x'.*A_sub_arr(l,:)))/(n*(n-1)-sum(x.*n_sub));
    end
    delta = norm(beta_hat_new-beta_hat)/norm(beta_hat)
    iter
end

pi_mat = zeros(length(Y), length(Y));
pi_mat(logical(triu(ones(length(Y)),1))) = x;
%imagesc(pi_mat)
%colorbar()

[I,J] = ind2sub(size(pi_mat), find(pi_mat>.95));
bd_start = Y(I); bd_end = Y(J);
pi_est = pi_mat(find(pi_mat>0.95));

%Aave = sum(A,1)/size(A,1);
%imagesc(reshape(Aave, size(Aave,2), size(Aave,3)))
%colorbar()
%hold on
%for i=1:length(bd_start)
%line([bd_start(i), bd_end(i)], [bd_start(i), bd_start(i)], 'Color','red', 'Linewidth',6)
%line([bd_end(i), bd_end(i)], [bd_start(i), bd_end(i)], 'Color','red', 'Linewidth',6)
%end
%hold off

end
