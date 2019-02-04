function prob = comp_prob(idx, alphabet_size, n, W)
% Compute the probability of a given configuration; used in the 'sampling_pariwise'.
x = ind2vec((dec2base(idx,alphabet_size,n)-'0')+1, alphabet_size);
x = reshape(full(x), [n*alphabet_size,1]);
prob = exp(x'*W*x/2);
