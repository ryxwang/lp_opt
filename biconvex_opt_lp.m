function [pi_mat,x,pi_est,bd_start,bd_end]=biconvex_opt_lp(A,Y,k)

n=size(A,1);

A_sub = zeros(length(Y), length(Y));
n_sub = zeros(length(Y), length(Y));


for a=1:(length(Y)-1)
    for b=(a+1):length(Y)
        A_sub(a,b) = sum(sum(A(Y(a):Y(b), Y(a):Y(b))));
        n_sub(a,b) = (Y(b)-Y(a)+1)*(Y(b)-Y(a)+1);
        if (A_sub(a,b)>n_sub(a,b))
            keyboard
        end
    end
end
A_sub = A_sub(find(~tril(ones(size(A_sub)))));
n_sub = n_sub(find(~tril(ones(size(n_sub)))));
m = length(A_sub);

alpha_hat = A_sub./(n_sub+1e-8);
alpha_hat(alpha_hat==1) = 1-1e-8;
alpha_hat(alpha_hat==0) = 1e-8;

%beta_hat = beta;
%alpha_hat = H(Y, Y);
%alpha_hat = nonzeros(triu(alpha_hat,1));
%alpha_hat = A_sub./n_sub;

lin_ind = zeros(length(Y), length(Y));
lin_ind(logical(triu(ones(length(Y)),1))) = 1:m;
B = zeros(length(Y), m);
for i=1:length(Y)
    B(i, nonzeros(lin_ind(1:i, i:length(Y)))) = 1;
end

%updating Pi
%initial beta
beta_hat_new = sum(sum(A)) /(n*(n-1));
delta = 1;
iter = 0;
eps=1e-5;
while delta > 1e-5
    iter = iter+1;
    beta_hat = beta_hat_new;
    t = log((eps+alpha_hat)./(1-alpha_hat+eps))-log((eps+beta_hat)/(eps+1-beta_hat));
    
    
    t2 = 0.5*(A_sub.*t + n_sub.*log((1-alpha_hat+eps)/(1-beta_hat+eps)));
   % keyboard
   
   %keyboard
    A1(1:size(B,1),:)=B;
    A1(size(B,1)+1,:)=ones(1,length(t2));
    b(1:size(B,1))=1;
    b(size(B,1)+1)=k;
    x = linprog(-t2,A1,b,[],[],zeros(1,length(t2)),1+zeros(1,length(t2)));
       
    beta_hat_new = (sum(sum(A)) - sum(x.*A_sub))/(n*(n-1)-sum(x.*n_sub));
    delta = abs(beta_hat_new-beta_hat)
    iter
end

pi_mat = zeros(length(Y), length(Y));
pi_mat(logical(triu(ones(length(Y)),1))) = x;
imagesc(pi_mat)
colorbar()

[I,J] = ind2sub(size(pi_mat), find(pi_mat>.9));
bd_start = Y(I); bd_end = Y(J);
pi_est = pi_mat(find(pi_mat>0.9));

imagesc(A)
colorbar()
hold on
for i=1:length(bd_start)
line([bd_start(i), bd_end(i)], [bd_start(i), bd_start(i)], 'Color','red', 'Linewidth',6)
line([bd_end(i), bd_end(i)], [bd_start(i), bd_end(i)], 'Color','red', 'Linewidth',6)
end
hold off

%x1=x;
%opt=t2'*x;
%opt=sum(t2'*x + epsilon*(entr(x)+entr(1-x)));
end
