function [w_avg, g_norm]  = mirror_descent(X, y, T, W1)
% use mirror descent to solve a logistic regression
% problem s.t. the simplex constraint:
% x_{t+1} = x_t exp(-eta*grad), x_{t+1} = x_{t+1} / ||x_{t+1}||_1
d = size(X,2);
N = size(X,1);
eta = sqrt(2*log(d)/T)/(2*W1);
w = ones(d,1)/d;
w_sum = w;
for t = 1:T
    preds = logsig(X*w)-y;
    grad = (ones(1,N)/N)*(repmat(preds,1,d).*X); % to speed up, one can compute stochastic gradient
    w = w.*exp(-eta*grad');
    w = w/sum(w);
    w_sum = w_sum+w;
end
w_avg = w_sum/T;
% compute the gradient norm at w_avg
preds = logsig(X*w_avg)-y;
grad = (ones(1,N)/N)*(repmat(preds,1,d).*X);
g_norm = norm(grad);
