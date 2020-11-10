## Obtain Information matrix from EM algorithm

Let $\theta\in {\Theta}$ be the parameters to be estimated by MLE and $f(x|\theta)$ denote density function. We observe a measureble fucntion of $x$, denoted as $y$, and the MLE($\hat\theta$) is to be found using $y$.  

The information matrix is not immediatly available from EM interactions, unlike the traditional one derived from score function. 

Let 
$$
l(x,\theta) = \log f(x|\theta), 
\\
l(y,\theta) = \log f(y|\theta) = \log\int f(x|\theta)dz.
$$
Here, $x$ is the complete data, for example, $x=(y,z)$ where $z$ is unobserved. 

Then the gradient of $l(y,\theta)$ with respect to $\theta$ is
$$
s(y,\theta) = l'(y,\theta) = \frac{\int f'(x|\theta)dz}{\int f(x|\theta)dz}.
$$
Multiply the integrand by $f(x|\theta)$, we have 
$$
s(y,\theta) =\frac{\int f'(x|\theta)\frac{f(x|\theta)}{f(x|\theta)}dz}{\int f(x|\theta)dz} = \int f'(x|\theta)/f(x|\theta)*f(z|y,\theta)dz = E_{z|y,\theta}s(x,\theta),
$$
where $s(x,\theta) = l'(x,\theta) = f'(x|\theta)/f(x|\theta)$.

Information matrix is the variance of score function, or the negative expected second derivative of log-likelihood.

The second derivative of log-likelihood is 
$$
\begin{split}
s'(y,\theta) &= \frac{\int f''(x|\theta)dz}{\int f(x|\theta)dz} - s(y,\theta)s^T(y,\theta)
\\&= E_{z|y,\theta}(\frac{f''(x|\theta)}{f(x|\theta)}) - s(y,\theta)s^T(y,\theta)
\\
&= E_{z|y,\theta}(l(x,\theta)'')+E_{z|y,\theta}(s(x,\theta)s^T(x,\theta))-s(y,\theta)s^T(y,\theta)
\end{split}
$$
This is the same as the chain rule of information matrix,
$$
I_{y,z}(\theta) = I_y(\theta)+I_{z|y}(\theta).
$$


### Gaussian mixture model example

Let $y_i\sim \sum_k\pi_k N(0,V_k+S_i+\sigma^2I)$, assume all but $\sigma^2$ is known. We are interested in estimating the variance of $\hat\sigma^2$ from MLE. One way is to obtain the information of $\sigma^2$ and take inverse of it, based on asymptotic result of MLE. 

The latent vairables are 1. $z_i$, a vector indicating which mixture $y_i$ is from; 2. $\mu_i\sim N(0,V_k)$; 3. $a_i\sim N(0,\sigma^2 I)$

Tbd