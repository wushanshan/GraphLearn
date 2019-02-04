This folder provides matlab files to reproduce Figure 2 of our paper.

The main function is 'pairwise_learn.m' file. It has a 'l21_lr' parameter, which decides whether to run l21-constrained logistic regression or the [Sparsitron algorithm](https://arxiv.org/abs/1706.06274).

Below is a snippet of running my code:

![pairwise_run](https://github.com/wushanshan/GraphLearn/blob/master/NonBinaryModel/pairwise_run.png)

Here 'err' measures the infinity norm of |W-W*| matrix, and 'succ' = 1 means the graph is exactly recovered, and 'succ' = 0 if not recovered. Note that for the case of k (the alphabet size) bigger than 2, to achieve 'succ'=1, 'err' does not necessarily be smaller than 0.1 in the above example.
