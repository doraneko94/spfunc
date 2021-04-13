# Examples
## gamma.rs
Examples of how to use the family of the gamma functions.
## zeta.rs
Examples of how to use the family of the zeta functions.
## mle_gamma.rs
Example of obtaining the parameters of the gamma distribution 
by Maximum Likelihood Estimation (MLE).

In this example, the gamma distribution is defined as 

![\begin{align*}
f(x|a,b)\equiv\frac{x^{a-1}}{\Gamma(a)b^{a}}\exp{\left(-\frac{x}{b}\right)}
\end{align*}
](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cbegin%7Balign%2A%7D%0Af%28x%7Ca%2Cb%29%5Cequiv%5Cfrac%7Bx%5E%7Ba-1%7D%7D%7B%5CGamma%28a%29b%5E%7Ba%7D%7D%5Cexp%7B%5Cleft%28-%5Cfrac%7Bx%7D%7Bb%7D%5Cright%29%7D%0A%5Cend%7Balign%2A%7D%0A)

where `a` is the shape and `b` is the scale.
 
For a description of the algorithm used in this code, see

https://tminka.github.io/papers/minka-gamma.pdf