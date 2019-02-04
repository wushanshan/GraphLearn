function [c,str] = proj_solve(p, s, eps)
% Input: p, positive scalar; s, n-dim positive array, eps: positive scalar
% Output: c, n-dim non-negative array
% We use Lagrange multiplier to solve the following constrained problem:
% min_c \sum_{i=1}^n c_i^p - s_i*c_i,  
%  s.t. \sum_{i=1}^n c_i \le 1, c_i >= 0 for all i=1,2...,n
% The Lagrange is: 
% L(c, \lambda, h)=\sum_i c_i^p - s_i*c_i-h_i*c_i + \lambda(\sum_i c_i-1)
% The KKT conditions are:
% pc_i^{p-1} - s_i - h_i + \lambda = 0;   
% \lambda >= 0, h_i >= 0,  h_i*c_i=0,  \lambda(\sum_i c_i - 1)=0
% ------------------------------------------------------------------
% Case 1: h_i^* = 0 and \lambda^* = 0
c = (s./p).^(1.0/(p-1));
if sum(c) <= 1
    str = 'case1';
    return
end
% Case 2: h_i^* = 0 and \lambda^* > 0
st = sort(s, 'descend');
lambda = st(end);
c = ((s-lambda)./p).^(1.0/(p-1));
if sum(c) <= 1
    % search between [0, st(end)] to find lambda^* s.t. sum(c)=1
    left = st(end);
    right = 0;
    lambda = (left+right)/2;
    c = ((s-lambda)./p).^(1.0/(p-1));
    while sum(c)-1<-eps || sum(c)-1>eps
        if sum(c)-1<-eps
            left = lambda;
        else
            right = lambda;
        end
        lambda = (left+right)/2;
        c = ((s-lambda)./p).^(1.0/(p-1));
    end
    str = 'case2';
    return
end
% Case 3: h_i^* > 0 and \lambda^* > 0
% Search for an index in [n] s.t. if lambda=st[i], then sum(c)<=1, 
% and if lambda=st[i+1], then sum(c)>=1. 
left = 1;
right = max(size(st));
index = floor((left+right)/2);
while right-left>1
    lambda = st(index);
    c = (max((s-lambda), 0)./p).^(1.0/(p-1));
    if sum(c) <= 1
        left = index;
    else
        right = index;
    end
    index = floor((left+right)/2);
end
% Search between [st(left), st(right)] to find lambda^* s.t. sum(c)=1
left = st(left);
right = st(right);
lambda = (left+right)/2;
c = (max((s-lambda),0)./p).^(1.0/(p-1));
while sum(c)-1<-eps || sum(c)-1>eps
    if sum(c)-1<-eps
        left = lambda;
    else
        right = lambda;
    end
    lambda = (left+right)/2;
    c = (max((s-lambda),0)./p).^(1.0/(p-1));
end
str = 'case3';