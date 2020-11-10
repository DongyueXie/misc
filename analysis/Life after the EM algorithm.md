# Life after the EM algorithm

Explore EM-algorithm and variational Bayes.

## AN ALTERNATIVE VIEW OF THE EM ALGORITHM

The log-likelihood can be written as 
$$
\begin{equation}
\begin{split}
\log p(x;\theta) &= \log \frac{p(x,z;\theta)}{p(z|x;\theta)}
\\
&= \int q(z)\log \frac{p(x,z;\theta)}{p(z|x;\theta)}d(z)
\\
&= \int q(z)\log \frac{p(x,z;\theta)q(z)}{p(z|x;\theta)q(z)}d(z)
\\
&= \int q(z)\log \frac{p(x,z;\theta)}{q(z)}dz - \int q(z)\log\frac{p(z|x;\theta)}{q(z)}dz
\end{split}
\end{equation}
$$
The second term in the last line is the KL-divergence between $p(z|x;\theta)$ and $q(z)$, $KL(q(z)||p(z|x;\theta)) = - \int q(z)\log\frac{p(z|x;\theta)}{q(z)}dz$.

The first term in the last line is called evidence lower bound(ELBO), $F(q,\theta) = \int q(z)\log \frac{p(x,z;\theta)}{q(z)}dz$. Why the name? In Byaesian statistics, marginal likelihood is usually called evidence. The KL divergence is always non-negative.  So we have $\log p(x;\theta)\geq F(q,\theta)$. Obviously, $F(q,\theta)$ is a lower bound the evidence and the equality holds when $q(z) = p(z|x,\theta)$.

EM algorithm maximizes the lower bound $F(x,\theta)$: in E-step, it's maximized wrt $q(z)$ and in M-step, it's maximized wrt $\theta$. Note that the increase in the log-likelihood is greater than the increase in the lower bound(because KL is not zero since it uses old $\theta$). 

EM algorirhm requires $p(z|x;\theta)$ being explicitly known.

##  VARIATIONAL EM FRAMEWORK

To by pass the requirement of knowing $p(z|x;\theta)$ explicitly, we instead assuming an approximate of $q(z)$. 

The most popular approximation is the factorized one. The hidden variables $z$ are assumed to be partitioned into M partitions, $z_i, i=1,2,...,M$ and $q(z) = \prod_i^M q_i(z_i)$.

Then the evidence lower bound is 
$$
\begin{equation}
\begin{split}
F &= \int \prod_iq_i\log p(x,z;\theta)\prod_i dz_i - \sum_i \int \prod_l q_l\log(q_i)dz
\\
&= \int q_j\left( \int\log(p(x,z;\theta))\prod_{i\neq j}q_idz_i \right)dz_j - \int \prod_lq_l\log q_j dz-\sum_{i\neq j} \int \prod_l q_l\log(q_i)dz
\\
&= \int q_j\left( \int\log(p(x,z;\theta))\prod_{i\neq j}q_idz_i \right)dz_j - \int q_j\log q_j dz_j - \sum_{i\neq j} \int q_i\log(q_i)dz_i
\\
&= -KL(q_j||\exp(\int\log(p(x,z;\theta))\prod_{i\neq j}q_idz_i)) + constant.of.q_j
\end{split}
\end{equation}
$$
 

Hence, the ELBO is maximized when we take $q_j(z_j) = \frac{\exp(\int\log(p(x,z;\theta))\prod_{i\neq j}q_idz_i)}{\int \exp(\int\log(p(x,z;\theta))\prod_{i\neq j}q_idz_i) dz_j}$.

The variational EM is summarized as 

E-step: update $q_j$ for $j=1,2,...,M$. 

M-step: update $\theta = argmax F(q,\theta)$.







 