%model selection
function [pi_est,bd_start,bd_end,k_opt] = select_k_combined(A,Y)
%options=optimoptions('quadprog');options.Display='off';options.TolX=1e-5;
delta=0;
k=1;
while delta < 1e-4 && k<length(Y)
    [pi_mat,x,temp_pi_est,temp_start,temp_end] = biconvex_opt_lp_combined(A,Y,k);
    delta = abs(sum(x)-k);
    if delta > 1e-4 && k>1
        break
    end
    pi_est=temp_pi_est;
    bd_start=temp_start;
    bd_end=temp_end;
    k = k+1;
end
k_opt=k-1;
