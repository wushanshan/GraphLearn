This folder provides matlab files that can be used to reproduce Figure 1 of our paper.

The main file is 'ising_learn.m'.

For example, if you have 4000 samples, and the diamond graph has 10 variables. Here is a snippet of running my code:

![ising_run](https://github.com/wushanshan/GraphLearn/IsingModel/ising_run.png)

The error is the infinity norm of matrix A-A_hat. In the above example, the error = 0.0801 < 0.2/2, which indicates successful recovery.
