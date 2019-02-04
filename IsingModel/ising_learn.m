function error = ising_learn(graph_size, num_sample, theta, rand_sign, num_iter)
% graph_size: indicates the size of the graph; theta: weight of each edge
% if rand_sign = True, then the sign of each edge weight is random
% if rand_sign = False, then all edge weights = theta
% num_samples: number of samples; num_iter: number of mirror descent iterations
% construct the adjacency matrix and perform sampling
s = graph_size;
[A, Xs] = sampling_diamond(s,theta,rand_sign,num_sample);
W1 = 2*(s-2)*theta;
A_hat = zeros(s,s);
for i = 1:s
    y = (Xs(:,i)+1)/2;
    X = [Xs(:,1:(i-1)), Xs(:,(i+1):s)];
    X = [X,-X,zeros(num_sample,1)]*W1;
    [w, ~] = mirror_descent(X, y, num_iter, W1);
    w = (w(1:s-1) - w(s:2*s-2))*W1;
    A_hat(i,1:i-1) = w(1:i-1)/2;
    A_hat(i,i+1:s) = w(i:end)/2;
end
error = max(max(abs(A-A_hat)));
