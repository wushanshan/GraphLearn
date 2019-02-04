function samples = sampling_pairwise(W, alphabet_size, num_samples)
% W is the weight matrix of the pairwise graphical model
% W has size n*k-by-n*k, where n is the graph size, k is the alphabet size
% Compute the probability of each configuration
n = floor(size(W,1)/alphabet_size);
probs = zeros(1, alphabet_size^n);
for i = 0:alphabet_size^n-1
    probs(i+1) = comp_prob(i, alphabet_size, n, W);
end
y = randsample(alphabet_size^n,num_samples,true,probs/sum(probs));
samples = (dec2base(y-1,alphabet_size,n)-'0')+1; % num_samples-by-n array