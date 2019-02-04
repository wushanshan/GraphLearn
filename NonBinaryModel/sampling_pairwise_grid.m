function [W, G, samples] = sampling_pairwise_grid(graph_size, theta, alphabet_size, num_samples)
% sampling from a pairwise graphical model on a 2D grid
% construct the weight matrix
n = graph_size; 
k = alphabet_size;
W = zeros(n*k,n*k); % W is the weight matrix
G = zeros(n, n);  % G is the adjacency matrix
w1 = ones(k, k)*theta;
for i = 1:k
    s = mod(i,2)+1;
    w1(i, s:2:k) = -theta;
end
w2 = ones(k, k)*theta;
for i = 1:k
    s = mod(i+1,2)+1;
    w2(i, s:2:k) = -theta;
end
for i=1:sqrt(n)
    for j=1:sqrt(n)
        idx = (i-1)*sqrt(n)+j;
        if j>1
            G(idx, idx-1) = 1;
            G(idx-1, idx) = 1;
            r = randi([0,1]);
            if r == 0
                W((idx*k-k+1):idx*k,((idx-1)*k-k+1):(idx-1)*k) = w1;
                W(((idx-1)*k-k+1):(idx-1)*k, (idx*k-k+1):idx*k) = w1;
            else
                W((idx*k-k+1):idx*k,((idx-1)*k-k+1):(idx-1)*k) = w2;
                W(((idx-1)*k-k+1):(idx-1)*k, (idx*k-k+1):idx*k) = w2;
            end
        end
        if i>1
            G(idx, idx-sqrt(n)) = 1;
            G(idx-sqrt(n), idx) = 1;
            r = randi([0,1]);
            if r == 0
                W((idx*k-k+1):idx*k, (idx-sqrt(n))*k-k+1:(idx-sqrt(n))*k) = w1;
                W((idx-sqrt(n))*k-k+1:(idx-sqrt(n))*k, (idx*k-k+1):idx*k) = w1;
            else
                W((idx*k-k+1):idx*k, (idx-sqrt(n))*k-k+1:(idx-sqrt(n))*k) = w2;
                W((idx-sqrt(n))*k-k+1:(idx-sqrt(n))*k, (idx*k-k+1):idx*k) = w2;
            end
        end
    end
end
% sampling
samples = sampling_pairwise(W, alphabet_size, num_samples);