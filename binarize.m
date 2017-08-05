function A = binarize(M, q)
thresh = quantile(M(triu(true(size(M)),0)), q);
A = (M>thresh)*1;
A(logical(eye(size(A)))) = 0;

