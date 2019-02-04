function samples = sampling_ising(A,N)
% A is the weight matrix, A_ij = A_ji, assuming the mean-field = 0
    s = size(A,1);
    probs = zeros(1,2^s);
    for i = 0:2^s-1
        x = (dec2bin(i,s)-'0')*2-1;
        probs(i+1) = exp(x*A*x'/2);
    end
    y = randsample(2^s,N,true,probs/sum(probs));
    samples = (dec2bin(y-1,s)-'0')*2-1;
end