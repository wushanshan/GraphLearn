function grad = grad_mirror(w)
% Compute the gradient of the following function
% Phi(w) = e\ln(n)/p \sum_{i}^n ||w(i,:)||_2^p
% where p = 1+1/ln(n)
[n, k] = size(w);
p = 1 + 1/log(n);
w_norm = sqrt(w.^2*ones(k,1))+1e-9;   % 1e-9 is added to avoid zero norm
grad = diag(exp(1)*log(n)*w_norm.^(p-2))*w;
