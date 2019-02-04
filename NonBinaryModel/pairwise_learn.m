function [err, graph_recovery] = pairwise_learn(W, G, samples, theta, num_iter, l21_lr)
% Learning non-binary Ising model on a 2D grid graph.
% l21_lr: 0 or 1, whether to use l21_constrained logistic regression.
% l21_lr = 0 means to use the Sparsitron algorithm proposed by KM17.
% Output: err = |W-W^*|_{\infty}, graph_recovery is a boolean variable
% indicating whether the original graph is recovered or not.
n = size(G, 1);
k = floor(size(W,1)/n);
W_hat = zeros(n*k, n*k);
for s = 1:n
    U = zeros(k*k,(n-1)*k);
    for i = 1:k
        for j = (i+1):k
            sam = samples(samples(:,s)==i | samples(:,s)==j,:);
            if s>1 && s<n
                Xs = sam(:, [1:s-1, s+1:n]);
            elseif s==1
                Xs = sam(:, 2:n);
            else
                Xs = sam(:, 1:n-1);
            end
            y = sam(:, s);
            y(y==i) = 1.0;
            y(y==j) = -1.0;
            N = size(Xs,1);
            X = zeros(n, k, N);
            for m1 = 1:N
                for m2 = 1:n-1
                    X(m2, Xs(m1, m2), m1) = 1.0;
                end
                X(n, 1, m1) = 1.0;
            end
            if n <= 4
                W21 = 2*sqrt(k)*2*theta;
                W1 = 2*k*2*theta;
            else
                W21 = 2*sqrt(k)*4*theta;
                W1 = 2*k*4*theta;
            end
            if l21_lr == 1
                w_avg = l21_mirror_descent(X, y, num_iter, W21);
            else
                w_avg = sparsitron(X, y, W1);
            end
            w_avg = w_avg - w_avg*ones(k,k)/k;   % centering each row
            w_avg = reshape(w_avg', [1, n*k]);
            idx = (i-1)*k+j;
            U(idx, :) = w_avg(1:(n-1)*k);
            idx = (j-1)*k+i;
            U(idx, :) = -w_avg(1:(n-1)*k);
        end
    end
    for m = 1:k
        if s>1 && s<n
            idx = [1:k*(s-1), k*s+1:n*k];
        elseif s==1
            idx = (k+1):n*k;
        else
            idx = 1:k*(n-1);
        end
        v = zeros(1, k*k);
        v((m-1)*k+1:m*k) = 1.0/k;
        W_hat((s-1)*k+m, idx) = v*U;
    end
end
err = max(max(abs(W-W_hat)));
G_hat = zeros(n, n);
for i = 1:n
    for j = 1:n
        if i~=j
            sub_mat = W_hat((i-1)*k+1:i*k, (j-1)*k+1:j*k);
            if max(max(abs(sub_mat)))>theta/2
                G_hat(i,j)=1;
            end
        end
    end
end
if max(max(abs(G-G_hat)))>0.01
    graph_recovery = 0;
else
    graph_recovery = 1;
end