function w_best = sparsitron(X, y, W1)
% Perform Sparsitron algorithm.
% Here x_i is the n-by-k-dim matrix, y = {-1,1}
% W1 is the l1 norm contraint on the weight vector.
% Transform X into vectors
[n, k, N] = size(X);
X = permute(X, [2 1 3]);
X = reshape(X, [n*k, N])';  % shape: [N, n*k]
X = [X, -X, zeros(N,1)];
y = (y+1)/2; % convert {-1,1} to {0,1}
s = max(100, floor(0.01*N));  % use the last 500 samples to select the w_best
X_tr = X(1:N-s, :); 
y_tr = y(1:N-s);
X_val = X(N-s+1:N,:); 
y_val = y(N-s+1:N,:);
d = size(X_tr,2);
beta = 1/(1+sqrt(log(d)/(N-s)));
w = ones(d,1)/d;
err_hat = (reshape(logsig(X_val*w), [1, s]) - reshape(y_val, [1,s])).^2;
err_hat = err_hat*ones(s,1)/s;
w_best = w;
best_err_hat = err_hat;
for t = 1:size(X_tr,1)
    Xs = X_tr(t, :);
    ys = y_tr(t);
    preds = logsig(Xs*w)-ys;
    grad = preds*Xs;
    w = w.*beta.^grad';
    w = w/sum(w)*W1;
    err_hat = (reshape(logsig(X_val*w), [1, s]) - reshape(y_val, [1,s])).^2;
    err_hat = err_hat*ones(s,1)/s;
    if err_hat<best_err_hat
        w_best = w;
    end
end
w_best = w_best(1:k*n) - w_best(k*n+1:2*n*k);
w_best = reshape(w_best, [k,n])';