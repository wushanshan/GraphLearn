function [A, samples] = sampling_diamond(s,theta,rand_sign,num_sample)
% sampling from an Ising model with diamond graph structure
% if rand_sign = True, then the sign of each edge weight is random
% if rand_sign = False, then all edge weights = theta
% construct the adjacency matrix
A = zeros(s,s);
A(1,2:s-1) = theta;
A(s,2:s-1) = theta;
A(2:s-1,1) = theta;
A(2:s-1,s) = theta;
if strcmp(rand_sign, 'True')==1
    sign_mat = triu(binornd(1,0.5,s,s)*2-1,1);
    A = A.*(sign_mat+sign_mat'); % need to make it symmetric
end
samples = sampling_ising(A, num_sample);
