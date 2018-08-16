%binarize a symmetric matrix M at the q-th quantile of upper triangular
%entries
function A = binarize(M, q)
thresh = quantile(M(triu(true(size(M)),0)), q);
A = (M>thresh)*1;
A(logical(eye(size(A)))) = 0;

