Brief overview of linear algebra
================
**Dmitry Kondrashov & Stefano Allesina**
Fundamentals of Biological Data Analysis – BIOS 26318

# Goals

  - Solving linear equations
  - Best-fit line through multiple data points
  - Concepts of linearity and vector spaces
  - Representation of vectors in multiple bases
  - Eigenvalues and eigenvectors of matrices

<!-- end list -->

``` r
library(tidyverse) # our friend the tidyverse
library(ggfortify) 
```

# Solving multivariate linear equations

Linear algebra is intimately tied to linear equations, that is, to
equations where all variables are multiplied by constant terms and added
together. Linear equations with one variable are easily solved by
division, for example:   
![
4x = 20
](https://latex.codecogs.com/png.latex?%0A4x%20%3D%2020%0A "
4x = 20
")  
is solved by dividing both sides by 4, obtaining the unique solution
![x=5](https://latex.codecogs.com/png.latex?x%3D5 "x=5").

The situation gets more interesting when multiple variables are
involved, with multiple equations, for example:   
![
\\begin{aligned}
4x - 3y &= 5 \\\\
\-x + 2y &= 10
\\end{aligned}
](https://latex.codecogs.com/png.latex?%0A%5Cbegin%7Baligned%7D%0A4x%20-%203y%20%26%3D%205%20%5C%5C%0A%20-x%20%2B%202y%20%26%3D%2010%0A%20%5Cend%7Baligned%7D%0A
"
\\begin{aligned}
4x - 3y &= 5 \\\\
 -x + 2y &= 10
 \\end{aligned}
")  

There are multiple ways to solve this, for example solving one equation
for one variable in terms of the other, then substituting it into the
second equation to obtain a one-variable problem. A more general
approach involves writing this problem in terms of a matrix ![\\mathbf
A](https://latex.codecogs.com/png.latex?%5Cmathbf%20A "\\mathbf A") that
contains the multiplicative constants of
![x](https://latex.codecogs.com/png.latex?x "x") and
![y](https://latex.codecogs.com/png.latex?y "y") and a vector ![\\vec
b](https://latex.codecogs.com/png.latex?%5Cvec%20b "\\vec b") that
contains the right-hand side constants:   
![
\\begin{aligned}
\\mathbf{A} = \\begin{pmatrix} 4 & -3 \\\\ -1 & 2 \\end{pmatrix}
\\;\\;\\; \\vec{b} = \\begin{pmatrix}5 \\\\ 10 \\end{pmatrix} &\\;\\;\\;
\\vec{x} = \\begin{pmatrix} x \\\\ y \\end{pmatrix} \\\\
\\mathbf{A} \\vec x = \\vec b &
\\end{aligned}
](https://latex.codecogs.com/png.latex?%0A%5Cbegin%7Baligned%7D%0A%5Cmathbf%7BA%7D%20%3D%20%5Cbegin%7Bpmatrix%7D%204%20%26%20-3%20%5C%5C%20-1%20%26%202%20%5Cend%7Bpmatrix%7D%20%5C%3B%5C%3B%5C%3B%20%5Cvec%7Bb%7D%20%3D%20%5Cbegin%7Bpmatrix%7D5%20%5C%5C%2010%20%5Cend%7Bpmatrix%7D%20%26%5C%3B%5C%3B%5C%3B%0A%5Cvec%7Bx%7D%20%3D%20%5Cbegin%7Bpmatrix%7D%20x%20%5C%5C%20y%20%5Cend%7Bpmatrix%7D%20%5C%5C%0A%5Cmathbf%7BA%7D%20%5Cvec%20x%20%3D%20%5Cvec%20b%20%26%0A%5Cend%7Baligned%7D%0A
"
\\begin{aligned}
\\mathbf{A} = \\begin{pmatrix} 4 & -3 \\\\ -1 & 2 \\end{pmatrix} \\;\\;\\; \\vec{b} = \\begin{pmatrix}5 \\\\ 10 \\end{pmatrix} &\\;\\;\\;
\\vec{x} = \\begin{pmatrix} x \\\\ y \\end{pmatrix} \\\\
\\mathbf{A} \\vec x = \\vec b &
\\end{aligned}
")  
This now looks analogous to the one-variable equation above, which we
solved by dividing both sides by the multiple of
![x](https://latex.codecogs.com/png.latex?x "x"). The difficulty is that
matrix operations are more complicated than scalar multiplication and
division. Matrix multiplication is used in the equation above to
multiply all the coefficients in the matrix by their respective
variables, which involves a relatively complicated [procedure in
general](https://en.wikipedia.org/wiki/Matrix_multiplication).

The “division” equivalent is called [matrix
inversion](https://en.wikipedia.org/wiki/Invertible_matrix) and it is
even more complicated. First, we need to define the identity matrix, or
the equivalent of the number 1 for matrix multiplication. The identity
matrix is defined only for square matrices (equal number of rows and
columns), so a size ![n](https://latex.codecogs.com/png.latex?n "n") by
![n](https://latex.codecogs.com/png.latex?n "n") idenitity matrix is
define to have all 1s on the diagonal and all zeros on the off-diagonal:
  
![
I = \\begin{pmatrix} 1 & 0 & \\dots & 0 \\\\ 0 & 1 &\\dots & 0 \\\\
\\vdots & \\vdots & \\ddots & \\vdots \\\\ 0 & 0 &\\dots
& 1\\end{pmatrix}
](https://latex.codecogs.com/png.latex?%0AI%20%3D%20%5Cbegin%7Bpmatrix%7D%201%20%26%200%20%26%20%5Cdots%20%26%200%20%5C%5C%200%20%26%201%20%20%26%5Cdots%20%26%200%20%5C%5C%20%5Cvdots%20%26%20%5Cvdots%20%26%20%5Cddots%20%26%20%5Cvdots%20%5C%5C%200%20%26%200%20%26%5Cdots%20%26%201%5Cend%7Bpmatrix%7D%0A
"
I = \\begin{pmatrix} 1 & 0 & \\dots & 0 \\\\ 0 & 1  &\\dots & 0 \\\\ \\vdots & \\vdots & \\ddots & \\vdots \\\\ 0 & 0 &\\dots & 1\\end{pmatrix}
")  
The identity matrix is special because multiplying any other matrix (of
compatible size) by it results in the same exact matrix (this is easy to
check on a couple of examples for 2 by 2 or 3 by 3 matrices):

  
![
I A = A I = A
](https://latex.codecogs.com/png.latex?%0AI%20A%20%3D%20A%20I%20%3D%20A%0A
"
I A = A I = A
")  
Then for an ![n](https://latex.codecogs.com/png.latex?n "n") by
![n](https://latex.codecogs.com/png.latex?n "n") matrix
![A](https://latex.codecogs.com/png.latex?A "A") its inverse
![A^{-1}](https://latex.codecogs.com/png.latex?A%5E%7B-1%7D "A^{-1}") is
defined to be the matrix multiplication by which results in the identity
matrix, that is:   
![
A^{-1} A = A A^{-1} = I
](https://latex.codecogs.com/png.latex?%0AA%5E%7B-1%7D%20A%20%3D%20A%20A%5E%7B-1%7D%20%20%3D%20I%0A
"
A^{-1} A = A A^{-1}  = I
")  
Defining the inverse is one task, but calculating it for any given
matrix, especially of large size, is quite laborious. We will not
describe the algorithms here, but you can read about [Gauss-Jordan
elimination](https://en.wikipedia.org/wiki/Gaussian_elimination#Finding_the_inverse_of_a_matrix),
which is one classic example. One important point is that not all
matrices are invertible, and for some no inverse matrix exists,
analogous to zero for real number division. The difference is that there
are infinitely many matrices for which this is the case, called
*singular* matrices.

In the cases in which the inverse matrix exists, the linear system of
equations can be solved by multiplying both sides by the inverse matrix,
like this:   
![
\\vec x = \\mathbf{A}^{-1}\\vec b
](https://latex.codecogs.com/png.latex?%0A%20%5Cvec%20x%20%3D%20%5Cmathbf%7BA%7D%5E%7B-1%7D%5Cvec%20b%0A
"
 \\vec x = \\mathbf{A}^{-1}\\vec b
")  

*Example:* Take the linear 2 by 2 system of equations of above and solve
it using matrix inversion.The R function solve() calculates the inverse
and multiplies it by the constant vector b:

``` r
A <- matrix(c(4,-1,-3,2), nrow = 2)
b <- c(5,10)
solve(A,b)
```

    # [1] 8 9

## Fitting a line to data

One geometric application of solving multiple linear equations is to
find the coefficients of a line that passes through two points in the
2-dimensional plane (or of a plane that passes through three points in
three-dimensional space, but we won’t go there.) In that case, the
coordinates of the points are the data, and the unknown variables are
the parameters slope ![m](https://latex.codecogs.com/png.latex?m "m")
and intercept ![b](https://latex.codecogs.com/png.latex?b "b") of the
line that we want to find.

**Example:** If the data set consists of two points ![(3,5),
(6, 2)](https://latex.codecogs.com/png.latex?%283%2C5%29%2C%20%286%2C%202%29
"(3,5), (6, 2)"), then finding the best fit values of
![m](https://latex.codecogs.com/png.latex?m "m") and
![b](https://latex.codecogs.com/png.latex?b "b") means solving the
following two equations:

  
![
\\begin{aligned}
3m + b &= 5 \\\\
6m + b &= 2
\\end{aligned}
](https://latex.codecogs.com/png.latex?%0A%5Cbegin%7Baligned%7D%0A%203m%20%2B%20b%20%26%3D%205%20%5C%5C%0A%206m%20%2B%20b%20%26%3D%202%0A%5Cend%7Baligned%7D%0A
"
\\begin{aligned}
 3m + b &= 5 \\\\
 6m + b &= 2
\\end{aligned}
")  

These equations have a solution for the slope and intercept, which can
be calculated in R using solve() and then plot the line with the
parameters from the solution vector beta:

``` r
xs <- c(3,6)
ys <- c(5,2)
A <- matrix(c(xs[1],xs[2],1,1), nrow = 2)
b <- c(ys[1],ys[2])
beta <- solve(A,b)
data1 <- tibble(xs,ys)
ggplot(data = data1) + aes(x=xs, y= ys) + geom_point() + geom_abline(slope=beta[1], intercept=beta[2])
```

<img src="linalg_basics_files/figure-gfm/unnamed-chunk-3-1.png" style="display: block; margin: auto;" />

However, a data set with two points is very small and nobody would
accept these values as reasonable estimates. Let us add one more data
point, to increase our sample size to three: ![(3,5, (6, 2),
(9, 1)](https://latex.codecogs.com/png.latex?%283%2C5%2C%20%286%2C%202%29%2C%20%289%2C%201%29
"(3,5, (6, 2), (9, 1)"). How do you find the best fit slope and
intercept?

**Bad idea:** take two points and find a line, that is the slope and the
intercept, that passes through the two. It should be clear why this is a
bad idea: we are arbitrarily ignoring some of the data, while perfectly
fitting two points. So how do we use all the data? Let us write down the
equations that a line with slope
![m](https://latex.codecogs.com/png.latex?m "m") and intercept
![b](https://latex.codecogs.com/png.latex?b "b") have to satisfy in
order to fit our data points:

  
![
\\begin{aligned}
3m + b = 5 \\\\
6m + b = 2\\\\
9m + b = 1
\\end{aligned}
](https://latex.codecogs.com/png.latex?%0A%5Cbegin%7Baligned%7D%0A3m%20%2B%20b%20%3D%205%20%5C%5C%0A6m%20%2B%20b%20%3D%202%5C%5C%0A9m%20%2B%20b%20%3D%201%0A%5Cend%7Baligned%7D%0A
"
\\begin{aligned}
3m + b = 5 \\\\
6m + b = 2\\\\
9m + b = 1
\\end{aligned}
")  

This system has no exact solution, since there are three equations and
only two unknowns. We need to find
![m](https://latex.codecogs.com/png.latex?m "m") and
![b](https://latex.codecogs.com/png.latex?b "b") such that they are a
“best fit” to the data, not the perfect solution.

## Least-squares line

Let us write the equation in matrix form as follows:

  
![
\\begin{aligned}
\\mathbf{A} = \\begin{pmatrix} 3 & 1 \\\\ 6 & 1 \\\\ 9 & 1
\\end{pmatrix} \\;\\;\\; \\vec{b} = \\begin{pmatrix}5 \\\\ 2 \\\\ 1
\\end{pmatrix} \\;\\;\\;
\\vec{\\beta} = \\begin{pmatrix} m \\\\ b \\end{pmatrix} & \\\\
\\mathbf{A} \\vec \\beta = \\vec b &
\\end{aligned}
](https://latex.codecogs.com/png.latex?%0A%5Cbegin%7Baligned%7D%0A%5Cmathbf%7BA%7D%20%3D%20%5Cbegin%7Bpmatrix%7D%203%20%26%201%20%5C%5C%206%20%26%201%20%5C%5C%209%20%26%201%20%5Cend%7Bpmatrix%7D%20%5C%3B%5C%3B%5C%3B%20%5Cvec%7Bb%7D%20%3D%20%5Cbegin%7Bpmatrix%7D5%20%5C%5C%202%20%5C%5C%201%20%5Cend%7Bpmatrix%7D%20%5C%3B%5C%3B%5C%3B%0A%5Cvec%7B%5Cbeta%7D%20%3D%20%5Cbegin%7Bpmatrix%7D%20m%20%5C%5C%20b%20%5Cend%7Bpmatrix%7D%20%26%20%5C%5C%0A%5Cmathbf%7BA%7D%20%5Cvec%20%5Cbeta%20%3D%20%5Cvec%20b%20%26%0A%5Cend%7Baligned%7D%0A
"
\\begin{aligned}
\\mathbf{A} = \\begin{pmatrix} 3 & 1 \\\\ 6 & 1 \\\\ 9 & 1 \\end{pmatrix} \\;\\;\\; \\vec{b} = \\begin{pmatrix}5 \\\\ 2 \\\\ 1 \\end{pmatrix} \\;\\;\\;
\\vec{\\beta} = \\begin{pmatrix} m \\\\ b \\end{pmatrix} & \\\\
\\mathbf{A} \\vec \\beta = \\vec b &
\\end{aligned}
")  

Mathematically, the problem is that one cannot invert a non-square
matrix. However, there is a way of turning the matrix into a square one,
by multiplying it by its own transpose (same matrix with rows and
columns reversed):

  
![
\\mathbf{A}^T \\mathbf{A} \\vec \\beta = \\mathbf{A}^T \\vec b
](https://latex.codecogs.com/png.latex?%0A%5Cmathbf%7BA%7D%5ET%20%5Cmathbf%7BA%7D%20%5Cvec%20%5Cbeta%20%3D%20%5Cmathbf%7BA%7D%5ET%20%5Cvec%20b%0A
"
\\mathbf{A}^T \\mathbf{A} \\vec \\beta = \\mathbf{A}^T \\vec b
")  
*Exercise:* Carry out the matrix multiplications to verify that
![\\mathbf{A}^T
\\mathbf{A}](https://latex.codecogs.com/png.latex?%5Cmathbf%7BA%7D%5ET%20%5Cmathbf%7BA%7D
"\\mathbf{A}^T \\mathbf{A}") is a 2 by 2 matrix and ![\\mathbf{A}^T
\\vec
b](https://latex.codecogs.com/png.latex?%5Cmathbf%7BA%7D%5ET%20%5Cvec%20b
"\\mathbf{A}^T \\vec b") is 2 by 1 vector.

Now we can solve this equation with a square matrix ![\\mathbf{A}^T
\\mathbf{A}](https://latex.codecogs.com/png.latex?%5Cmathbf%7BA%7D%5ET%20%5Cmathbf%7BA%7D
"\\mathbf{A}^T \\mathbf{A}") by multiplying both sides by the inverse\!
In general, for an ![n](https://latex.codecogs.com/png.latex?n
"n")-dimensional data set consisting of a bunch of values of
![x](https://latex.codecogs.com/png.latex?x "x") and
![y](https://latex.codecogs.com/png.latex?y "y"), the process loooks
like this:

  
![
\\vec Y = \\begin{pmatrix} y\_1\\\\ y\_2\\\\ \\vdots \\\\ y\_n
\\end{pmatrix} \\;\\;\\; 
\\mathbf{X} = \\begin{pmatrix} 1 & x\_1\\\\ 1 & x\_2\\\\ \\vdots &
\\vdots \\\\ 1 & x\_n \\end{pmatrix}
\\;\\;\\; 
\\mathbf{\\beta} = \\begin{pmatrix} m \\\\ b \\end{pmatrix} \\\\
\\beta = (\\mathbf{X}^{T} \\mathbf{X})^{-1} \\mathbf{X}^{T}\\vec Y
](https://latex.codecogs.com/png.latex?%0A%5Cvec%20Y%20%3D%20%5Cbegin%7Bpmatrix%7D%20y_1%5C%5C%20y_2%5C%5C%20%5Cvdots%20%5C%5C%20y_n%20%5Cend%7Bpmatrix%7D%20%5C%3B%5C%3B%5C%3B%20%0A%5Cmathbf%7BX%7D%20%3D%20%5Cbegin%7Bpmatrix%7D%201%20%26%20x_1%5C%5C%201%20%26%20x_2%5C%5C%20%5Cvdots%20%26%20%5Cvdots%20%5C%5C%201%20%26%20x_n%20%5Cend%7Bpmatrix%7D%0A%20%5C%3B%5C%3B%5C%3B%20%0A%5Cmathbf%7B%5Cbeta%7D%20%3D%20%5Cbegin%7Bpmatrix%7D%20m%20%5C%5C%20b%20%5Cend%7Bpmatrix%7D%20%5C%5C%0A%5Cbeta%20%3D%20%28%5Cmathbf%7BX%7D%5E%7BT%7D%20%5Cmathbf%7BX%7D%29%5E%7B-1%7D%20%5Cmathbf%7BX%7D%5E%7BT%7D%5Cvec%20Y%0A
"
\\vec Y = \\begin{pmatrix} y_1\\\\ y_2\\\\ \\vdots \\\\ y_n \\end{pmatrix} \\;\\;\\; 
\\mathbf{X} = \\begin{pmatrix} 1 & x_1\\\\ 1 & x_2\\\\ \\vdots & \\vdots \\\\ 1 & x_n \\end{pmatrix}
 \\;\\;\\; 
\\mathbf{\\beta} = \\begin{pmatrix} m \\\\ b \\end{pmatrix} \\\\
\\beta = (\\mathbf{X}^{T} \\mathbf{X})^{-1} \\mathbf{X}^{T}\\vec Y
")  
*Example:* Let us see the best-fit line for the 3-point data set above:

``` r
xs <- c(3,6,9)
ys <- c(5,2,1)
A <- matrix(c(xs[1],xs[2],xs[3],1,1,1), nrow = 3)
b <- c(ys[1],ys[2], ys[3])
beta <- solve(t(A) %*% A,t(A) %*% b)
data1 <- tibble(xs,ys)
ggplot(data = data1) + aes(x=xs, y= ys) + geom_point() + geom_abline(slope=beta[1], intercept=beta[2])
```

<img src="linalg_basics_files/figure-gfm/unnamed-chunk-4-1.png" style="display: block; margin: auto;" />

Let us use the classic data set of Karl Pearson’s from 1903 containing
the height of fathers and sons, which we will return to next week when
we tackle linear regression
properly:

``` r
heights <- read_tsv("http://www.randomservices.org/random/data/Pearson.txt")
pl <- ggplot(data = heights) + aes(x = Father, y = Son) + geom_point() + coord_equal()
pl
```

<img src="linalg_basics_files/figure-gfm/unnamed-chunk-5-1.png" style="display: block; margin: auto;" />

*Exercise:* Let’s try to find the best fit line to this data set (the
hard way) using the same process as above for the three - point data
set:

Of course, `R` can do this calculation for you with just one command:

``` r
best_beta_easy <- lm(Son ~ Father, data = heights)
best_beta_easy
```

    # 
    # Call:
    # lm(formula = Son ~ Father, data = heights)
    # 
    # Coefficients:
    # (Intercept)       Father  
    #      33.893        0.514

But it feels good to know that this is not black magic\! In fact,
plotting it on top of the data does not even require computing the
coefficients:

``` r
pl + geom_smooth(method = "lm") # lm stands for linear model
```

<img src="linalg_basics_files/figure-gfm/unnamed-chunk-8-1.png" style="display: block; margin: auto;" />

# Linearity and vector spaces

We have dealt with linear models in various guises, so now would be a
good time to define properly what linearity means. The word comes from
the shape of graphs of linear functions of one variable, e.g. ![f(x) =
ax +
b](https://latex.codecogs.com/png.latex?f%28x%29%20%3D%20ax%20%2B%20b
"f(x) = ax + b"), but the algebraic meaning rests on the following two
general properties:

**Definition.** A *linear transformation* or *linear operator* is a
mapping ![L](https://latex.codecogs.com/png.latex?L "L") between two
sets of vectors with the following properties:

1.  *(scalar multiplication)* ![L(c \\vec v) = c L(\\vec
    v)](https://latex.codecogs.com/png.latex?L%28c%20%5Cvec%20v%29%20%3D%20c%20L%28%5Cvec%20v%29
    "L(c \\vec v) = c L(\\vec v)"); where
    ![c](https://latex.codecogs.com/png.latex?c "c") is a scalar and
    ![\\vec v](https://latex.codecogs.com/png.latex?%5Cvec%20v
    "\\vec v") is a vector
2.  *(additive)* ![L(\\vec v\_1 + \\vec v\_2) = L(\\vec v\_1) + L(\\vec
    v\_2)](https://latex.codecogs.com/png.latex?L%28%5Cvec%20v_1%20%2B%20%5Cvec%20v_2%29%20%3D%20L%28%5Cvec%20v_1%29%20%2B%20L%28%5Cvec%20v_2%29
    "L(\\vec v_1 + \\vec v_2) = L(\\vec v_1) + L(\\vec v_2)"); where
    ![\\vec v\_1](https://latex.codecogs.com/png.latex?%5Cvec%20v_1
    "\\vec v_1") and ![\\vec
    v\_2](https://latex.codecogs.com/png.latex?%5Cvec%20v_2 "\\vec v_2")
    are vectors

Here we have two types of objects: vectors and transformations/operators
that act on those vectors. The basic example of this are vectors and
matrices, because a matrix multiplied by a vector (on the right) results
another vector, provided the number of columns in the matrix is the same
as the number of rows in the vector. This can be interpreted as the
matrix transforming the vector ![\\vec
v](https://latex.codecogs.com/png.latex?%5Cvec%20v "\\vec v") into
another one: $ A v = u$.

**Example:** Let us multiply the following matrix and vector (specially
chosen to make a point):

``` r
A <- matrix(c(2, 2, 1, 3), nrow = 2)
vec1 <- c(1,-1)
vec2 <- A %*% vec1
print(vec1)
print(vec2)
```

    # [1]  1 -1
    #      [,1]
    # [1,]    1
    # [2,]   -1

We see that this particular vector
![(1,-1)](https://latex.codecogs.com/png.latex?%281%2C-1%29 "(1,-1)") is
unchanged when multiplied by this matrix, or we can say that the matrix
multiplication is equivalent to multiplication by 1. Here is another
such vector for the same matrix:

``` r
vec1 <- c(1,2)
vec2 <- A %*% vec1
print(vec1)
print(vec2)
```

    # [1] 1 2
    #      [,1]
    # [1,]    4
    # [2,]    8

In this case, the vector is changed, but only by multiplication by a
constant (4). Thus the geometric direction of the vector remained
unchanged.

The notion of linearity leads to the important idea of combining
different vectors:

**Definition:** A *linear combination* of
![n](https://latex.codecogs.com/png.latex?n "n") vectors ![\\{ \\vec
v\_i
\\}](https://latex.codecogs.com/png.latex?%5C%7B%20%5Cvec%20v_i%20%5C%7D
"\\{ \\vec v_i \\}") is a weighted sum of these vectors with any real
numbers
![\\{a\_i\\}](https://latex.codecogs.com/png.latex?%5C%7Ba_i%5C%7D
"\\{a_i\\}"):   
![ a\_1 \\vec v\_1+ a\_2 \\vec v\_2... + a\_n \\vec
v\_n](https://latex.codecogs.com/png.latex?%20a_1%20%5Cvec%20v_1%2B%20a_2%20%5Cvec%20v_2...%20%2B%20a_n%20%5Cvec%20v_n
" a_1 \\vec v_1+ a_2 \\vec v_2... + a_n \\vec v_n")  

Linear combinations arise naturally from the notion of linearity,
combining the additive property and the scalar multiplication property.
Speaking intuitively, a linear combination of vectors produces a new
vector that is related to the original set. Linear combinations give a
simple way of generating new vectors, and thus invite the following
definition for a collection of vectors closed under linear combinations:

**Definition.** A *vector space* is a collection of vectors such that a
linear combination of any ![n](https://latex.codecogs.com/png.latex?n
"n") vectors is contained in the vector space.

The most common examples are the spaces of all real-valued vectors of
dimension ![n](https://latex.codecogs.com/png.latex?n "n"), which are
denoted by
![\\mathbb{R}^n](https://latex.codecogs.com/png.latex?%5Cmathbb%7BR%7D%5En
"\\mathbb{R}^n"). For instance,
![\\mathbb{R}^2](https://latex.codecogs.com/png.latex?%5Cmathbb%7BR%7D%5E2
"\\mathbb{R}^2") (pronounced “r two”) is the vector space of two
dimensional real-valued vectors such as
![(1,3)](https://latex.codecogs.com/png.latex?%281%2C3%29 "(1,3)") and
![(\\pi,
-\\sqrt{17})](https://latex.codecogs.com/png.latex?%28%5Cpi%2C%20-%5Csqrt%7B17%7D%29
"(\\pi, -\\sqrt{17})"); similarly,
![\\mathbb{R}^3](https://latex.codecogs.com/png.latex?%5Cmathbb%7BR%7D%5E3
"\\mathbb{R}^3") is the vector space consisting of three dimensional
real-valued vectors such as
![(0.1,0,-5.6)](https://latex.codecogs.com/png.latex?%280.1%2C0%2C-5.6%29
"(0.1,0,-5.6)"). You can convince yourself, by taking linear
combinations of vectors, that these vector spaces contain all the points
in the usual Euclidean plane and three-dimensional space. The real
number line can also be thought of as the vector space
![\\mathbb{R}^1](https://latex.codecogs.com/png.latex?%5Cmathbb%7BR%7D%5E1
"\\mathbb{R}^1").

## Linear independence and basis vectors

How can we describe a vector space without trying to list all of its
elements? We know that one can generate an element by taking linear
combinations of vectors. It turns out that it is possible to generate
(or “span”) a vector space by taking linear combinations of a subset of
its vectors. The challenge is to find a minimal subset of subset that is
not redundant. In order to do this, we first introduce a new concept:

**Definition:** A set of vectors ![\\{ \\vec v\_i
\\}](https://latex.codecogs.com/png.latex?%5C%7B%20%5Cvec%20v_i%20%5C%7D
"\\{ \\vec v_i \\}") is called *linearly independent* if the only linear
combination involving them that equals the zero vector is if all the
coefficients are zero. ( ![a\_1 \\vec v\_1 + a\_2 \\vec v\_2 + ... +
a\_n \\vec v\_n
= 0](https://latex.codecogs.com/png.latex?a_1%20%5Cvec%20v_1%20%2B%20a_2%20%5Cvec%20v_2%20%2B%20...%20%2B%20a_n%20%5Cvec%20v_n%20%3D%200
"a_1 \\vec v_1 + a_2 \\vec v_2 + ... + a_n \\vec v_n = 0") only if
![a\_i = 0](https://latex.codecogs.com/png.latex?a_i%20%3D%200
"a_i = 0") for all ![i](https://latex.codecogs.com/png.latex?i "i").)

In the familiar Euclidean spaces, e.g.
![\\mathbb{R}^2](https://latex.codecogs.com/png.latex?%5Cmathbb%7BR%7D%5E2
"\\mathbb{R}^2"), linear independence has a geometric meaning: two
vectors are linearly independent if the segments from the origin to the
endpoint do not lie on the same line. But it can be shown that any set
of three vectors in the plane is linearly dependent, because there are
only two dimensions in the vector space. This brings us to the key
definition of this section:

**Definition:** A *basis* of a vector space is a linearly independent
set of vectors that generate (or span) the vector space. The number of
vectors (cardinality) in such a set is called the *dimension* of the
vector space.

A vector space generally has many possible bases, as illustrated in
figure. In the case of
![\\mathbb{R}^2](https://latex.codecogs.com/png.latex?%5Cmathbb%7BR%7D%5E2
"\\mathbb{R}^2"), the usual (canonical) basis set is ![\\{(1,0);
(0,1)\\}](https://latex.codecogs.com/png.latex?%5C%7B%281%2C0%29%3B%20%280%2C1%29%5C%7D
"\\{(1,0); (0,1)\\}") which obviously generates any point on the plane
and is linearly independent. But any two linearly independent vectors
can generate any vector in the plane.

**Example:** The vector ![\\vec r =
(2,1)](https://latex.codecogs.com/png.latex?%5Cvec%20r%20%3D%20%282%2C1%29
"\\vec r = (2,1)") can be represented as a linear combination of the two
canonical vectors: ![\\vec r = 2\\times (1,0)+1\\times
(0,1)](https://latex.codecogs.com/png.latex?%5Cvec%20r%20%3D%202%5Ctimes%20%281%2C0%29%2B1%5Ctimes%20%280%2C1%29
"\\vec r = 2\\times (1,0)+1\\times (0,1)"). Let us choose another basis
set, say ![\\{(1,1);
(-1,1)\\}](https://latex.codecogs.com/png.latex?%5C%7B%281%2C1%29%3B%20%28-1%2C1%29%5C%7D
"\\{(1,1); (-1,1)\\}") (this is the canonical basis vectors rotated by
![\\pi/2](https://latex.codecogs.com/png.latex?%5Cpi%2F2 "\\pi/2").) The
same vector can be represented by a linear combination of these two
vectors, with coefficients
![1.5](https://latex.codecogs.com/png.latex?1.5 "1.5") and
![-0.5](https://latex.codecogs.com/png.latex?-0.5 "-0.5"): ![\\vec r
= 1.5\\times (1,1) - 0.5 \\times
(-1,1)](https://latex.codecogs.com/png.latex?%5Cvec%20r%20%3D%201.5%5Ctimes%20%281%2C1%29%20-%200.5%20%5Ctimes%20%28-1%2C1%29
"\\vec r = 1.5\\times (1,1) - 0.5 \\times (-1,1)"). If we call the first
basis ![C](https://latex.codecogs.com/png.latex?C "C") for canonical and
the second basis ![D](https://latex.codecogs.com/png.latex?D "D") for
different, we can write the same vector using different sets of
coordinates for each basis:   
![ 
\\vec r\_{C} = (2,1); \\; \\vec r\_D = (1.5, -0.5)
](https://latex.codecogs.com/png.latex?%20%0A%5Cvec%20r_%7BC%7D%20%3D%20%282%2C1%29%3B%20%5C%3B%20%5Cvec%20r_D%20%3D%20%281.5%2C%20-0.5%29%0A
" 
\\vec r_{C} = (2,1); \\; \\vec r_D = (1.5, -0.5)
")  

## Projections and changes of basis

The representation of an arbitrary vector (point) in a vector space as a
linear combination of a given basis set is called the  of the point in
terms of the basis, which gives the coordinates for the vector in terms
of each basis vector. The decomposition of a point in terms of a
particular basis is very useful in high-dimensional spaces, where a
clever choice of a basis can allow a description of a set of points
(such as a data set) in terms of contributions of only a few basis
vectors, if the data set primarily extends only in a few dimensions.

To obtain the coefficients of the basis vectors in a decomposition of a
vector ![\\vec r](https://latex.codecogs.com/png.latex?%5Cvec%20r
"\\vec r"), we need to perform what is termed a *projection* of the
vector onto the basis vectors. Think of shining a light perpendicular to
the basis vector, and measuring the length of the shadow cast by the
vector ![\\vec r](https://latex.codecogs.com/png.latex?%5Cvec%20r
"\\vec r") onto ![\\vec
v\_i](https://latex.codecogs.com/png.latex?%5Cvec%20v_i "\\vec v_i"). If
the vectors are parallel, the shadow is equal to the length of ![\\vec
r](https://latex.codecogs.com/png.latex?%5Cvec%20r "\\vec r"); if they
are orthogonal, the shadow is nonexistent. To find the length of the
shadow, use the inner product of ![\\vec
r](https://latex.codecogs.com/png.latex?%5Cvec%20r "\\vec r") and
![\\vec v](https://latex.codecogs.com/png.latex?%5Cvec%20v "\\vec v"),
which as you recall corresponds to the cosine of the angle between the
two vectors multiplied by their norms: $r, v=rv() $. We do not care
about the length of the vector ![\\vec
v](https://latex.codecogs.com/png.latex?%5Cvec%20v "\\vec v") we are
projecting onto, thus we divide the inner product by the square norm of
![\\vec v](https://latex.codecogs.com/png.latex?%5Cvec%20v "\\vec v"),
and then multiply the vector ![\\vec
v](https://latex.codecogs.com/png.latex?%5Cvec%20v "\\vec v") by this
projection coefficient:   
![ 
Proj(\\vec r ; \\vec v) = \\frac{ \\langle \\vec r , \\vec v \\rangle }
{\\langle \\vec v , \\vec v \\rangle } \\vec v = \\frac{ \\langle \\vec
r , \\vec v \\rangle } {\\vert \\vec v \\vert^2} \\vec v= \\frac{
\\vert\\vec r\\vert \\cos(\\theta) } {\\vert \\vec v \\vert}\\vec v
](https://latex.codecogs.com/png.latex?%20%0AProj%28%5Cvec%20r%20%3B%20%5Cvec%20v%29%20%3D%20%5Cfrac%7B%20%5Clangle%20%5Cvec%20r%20%2C%20%5Cvec%20v%20%5Crangle%20%20%7D%20%7B%5Clangle%20%5Cvec%20v%20%2C%20%5Cvec%20v%20%5Crangle%20%7D%20%5Cvec%20v%20%3D%20%5Cfrac%7B%20%5Clangle%20%5Cvec%20r%20%2C%20%20%5Cvec%20v%20%5Crangle%20%20%7D%20%7B%5Cvert%20%5Cvec%20v%20%5Cvert%5E2%7D%20%5Cvec%20v%3D%20%5Cfrac%7B%20%20%5Cvert%5Cvec%20r%5Cvert%20%5Ccos%28%5Ctheta%29%20%7D%20%7B%5Cvert%20%5Cvec%20v%20%5Cvert%7D%5Cvec%20v%0A
" 
Proj(\\vec r ; \\vec v) = \\frac{ \\langle \\vec r , \\vec v \\rangle  } {\\langle \\vec v , \\vec v \\rangle } \\vec v = \\frac{ \\langle \\vec r ,  \\vec v \\rangle  } {\\vert \\vec v \\vert^2} \\vec v= \\frac{  \\vert\\vec r\\vert \\cos(\\theta) } {\\vert \\vec v \\vert}\\vec v
")  

This formula gives the projection of the vector ![\\vec
r](https://latex.codecogs.com/png.latex?%5Cvec%20r "\\vec r") onto
![\\vec v](https://latex.codecogs.com/png.latex?%5Cvec%20v "\\vec v"),
the result is a new vector in the direction of ![\\vec
v](https://latex.codecogs.com/png.latex?%5Cvec%20v "\\vec v"), with the
scalar coefficient ![a = \\ \\langle \\vec r , \\vec v \\rangle /\\vert
\\vec v
\\vert^2](https://latex.codecogs.com/png.latex?a%20%3D%20%5C%20%5Clangle%20%5Cvec%20r%20%2C%20%5Cvec%20v%20%5Crangle%20%2F%5Cvert%20%5Cvec%20v%20%5Cvert%5E2
"a = \\ \\langle \\vec r , \\vec v \\rangle /\\vert \\vec v \\vert^2").

**Example:** Here is how one might calculate the projection of the point
![(2,1)](https://latex.codecogs.com/png.latex?%282%2C1%29 "(2,1)") onto
the basis set ![\\{(1,1);
(-1,1)\\}](https://latex.codecogs.com/png.latex?%5C%7B%281%2C1%29%3B%20%28-1%2C1%29%5C%7D
"\\{(1,1); (-1,1)\\}"):

``` r
v1 <- c(1,1)
v2 <- c(-1,1)
u <- c(2,1)
ProjMat <- matrix(cbind(v1,v2), byrow = T, nrow = 2)
print(ProjMat)
ProjMat%*%u
```

    #      [,1] [,2]
    # [1,]    1    1
    # [2,]   -1    1
    #      [,1]
    # [1,]    3
    # [2,]   -1

This is not quite right: the projection coefficients are off by a factor
of two compared to the correct values in the example above. This is
because we have neglected to *normalize* the basis vectors, so we should
modify the script as follows:

``` r
v1 <- c(1,1)
v1 <- v1/(sum(v1^2))
v2 <- c(-1,1)
v2 <- v2/(sum(v2^2))
u <- c(2,1)
ProjMat <- matrix(cbind(v1,v2), byrow = T, nrow = 2)
print(ProjMat)
print(ProjMat%*%u)
```

    #      [,1] [,2]
    # [1,]  0.5  0.5
    # [2,] -0.5  0.5
    #      [,1]
    # [1,]  1.5
    # [2,] -0.5

This is an example of how to convert a vector/point from representation
in one basis set to another. The new basis vectors, expressed in the
original basis set, are arranged in a matrix by row, scaled by their
norm squared, and multiplied by the vector that one wants to express in
the new basis. The resulting vector contains the coordinates in the new
basis.

# Matrices as linear operators

## matrices transform vectors

In this section we will learn to characterize square matrices by finding
special numbers and vectors associated with them. At the core of this
analysis lies the concept of a matrix as an *operator* that transforms
vectors by multiplication. To be clear, in this section we take as
default that the matrices ![A](https://latex.codecogs.com/png.latex?A
"A") are square, and that vectors ![\\vec
v](https://latex.codecogs.com/png.latex?%5Cvec%20v "\\vec v") are column
vectors, and thus will multiply the matrix on the right: ![A \\times
\\vec v](https://latex.codecogs.com/png.latex?A%20%5Ctimes%20%5Cvec%20v
"A \\times \\vec v").

A matrix multiplied by a vector produces another vector, provided the
number of columns in the matrix is the same as the number of rows in the
vector. This can be interpreted as the matrix transforming the vector
![\\vec v](https://latex.codecogs.com/png.latex?%5Cvec%20v "\\vec v")
into another one: $ A v = u$. The resultant vector ![\\vec
u](https://latex.codecogs.com/png.latex?%5Cvec%20u "\\vec u") may or may
not resemble ![\\vec v](https://latex.codecogs.com/png.latex?%5Cvec%20v
"\\vec v"), but there are special vectors for which the transformation
is very simple.

**Example.** Let us multiply the following matrix and vector (specially
chosen to make a point):   
![
\\left(\\begin{array}{cc}2 & 1
\\\\ 2& 3\\end{array}\\right)\\left(\\begin{array}{c}1 \\\\ -1
\\end{array}\\right) = \\left(\\begin{array}{c}2 -1 \\\\ 2 - 3
\\end{array}\\right) = \\left(\\begin{array}{c} 1 \\\\ -1
\\end{array}\\right)
](https://latex.codecogs.com/png.latex?%0A%5Cleft%28%5Cbegin%7Barray%7D%7Bcc%7D2%20%26%201%20%5C%5C%202%26%203%5Cend%7Barray%7D%5Cright%29%5Cleft%28%5Cbegin%7Barray%7D%7Bc%7D1%20%5C%5C%20-1%20%5Cend%7Barray%7D%5Cright%29%20%3D%20%5Cleft%28%5Cbegin%7Barray%7D%7Bc%7D2%20-1%20%5C%5C%202%20-%203%20%5Cend%7Barray%7D%5Cright%29%20%3D%20%20%5Cleft%28%5Cbegin%7Barray%7D%7Bc%7D%201%20%5C%5C%20-1%20%5Cend%7Barray%7D%5Cright%29%0A
"
\\left(\\begin{array}{cc}2 & 1 \\\\ 2& 3\\end{array}\\right)\\left(\\begin{array}{c}1 \\\\ -1 \\end{array}\\right) = \\left(\\begin{array}{c}2 -1 \\\\ 2 - 3 \\end{array}\\right) =  \\left(\\begin{array}{c} 1 \\\\ -1 \\end{array}\\right)
")  
We see that this particular vector is unchanged when multiplied by this
matrix, or we can say that the matrix multiplication is equivalent to
multiplication by 1. Here is another such vector for the same matrix:   
![
\\left(\\begin{array}{cc}2 & 1
\\\\ 2& 3\\end{array}\\right)\\left(\\begin{array}{c}1 \\\\ 2
\\end{array}\\right) = \\left(\\begin{array}{c}2 +2 \\\\ 2 + 6
\\end{array}\\right) = \\left(\\begin{array}{c} 4 \\\\ 8
\\end{array}\\right)
](https://latex.codecogs.com/png.latex?%0A%5Cleft%28%5Cbegin%7Barray%7D%7Bcc%7D2%20%26%201%20%5C%5C%202%26%203%5Cend%7Barray%7D%5Cright%29%5Cleft%28%5Cbegin%7Barray%7D%7Bc%7D1%20%5C%5C%202%20%5Cend%7Barray%7D%5Cright%29%20%3D%20%5Cleft%28%5Cbegin%7Barray%7D%7Bc%7D2%20%2B2%20%5C%5C%202%20%2B%206%20%5Cend%7Barray%7D%5Cright%29%20%3D%20%20%5Cleft%28%5Cbegin%7Barray%7D%7Bc%7D%204%20%5C%5C%208%20%5Cend%7Barray%7D%5Cright%29%0A
"
\\left(\\begin{array}{cc}2 & 1 \\\\ 2& 3\\end{array}\\right)\\left(\\begin{array}{c}1 \\\\ 2 \\end{array}\\right) = \\left(\\begin{array}{c}2 +2 \\\\ 2 + 6 \\end{array}\\right) =  \\left(\\begin{array}{c} 4 \\\\ 8 \\end{array}\\right)
")  
In this case, the vector is changed, but only by multiplication by a
constant (4). Thus the geometric direction of the vector remained
unchanged.

Generally, a square matrix has an associated set of vectors for which
multiplication by the matrix is equivalent to multiplication by a
constant. This can be written down as a definition:

An *eigenvector* of a square matrix
![A](https://latex.codecogs.com/png.latex?A "A") is a vector ![\\vec
v](https://latex.codecogs.com/png.latex?%5Cvec%20v "\\vec v") for which
matrix multiplication by ![A](https://latex.codecogs.com/png.latex?A
"A") is equivalent to multiplication by a constant. This constant
![\\lambda](https://latex.codecogs.com/png.latex?%5Clambda "\\lambda")
is called its *eigenvalue* of
![A](https://latex.codecogs.com/png.latex?A "A") corresponding the the
eigenvector ![\\vec v](https://latex.codecogs.com/png.latex?%5Cvec%20v
"\\vec v"). The relationship is summarized in the following equation:   
![
A \\times \\vec v = \\lambda \\vec v
](https://latex.codecogs.com/png.latex?%0AA%20%20%5Ctimes%20%20%5Cvec%20v%20%3D%20%5Clambda%20%5Cvec%20v%0A
"
A  \\times  \\vec v = \\lambda \\vec v
")  

Note that this equation combines a matrix
(![A](https://latex.codecogs.com/png.latex?A "A")), a vector (![\\vec
v](https://latex.codecogs.com/png.latex?%5Cvec%20v "\\vec v")) and a
scalar ![\\lambda](https://latex.codecogs.com/png.latex?%5Clambda
"\\lambda"), and that both sides of the equation are column vectors.

The definition does not specify how many such eigenvectors and
eigenvalues can exist for a given matrix
![A](https://latex.codecogs.com/png.latex?A "A"). There are usually as
many such vectors ![\\vec
v](https://latex.codecogs.com/png.latex?%5Cvec%20v "\\vec v") and
corresponding numbers
![\\lambda](https://latex.codecogs.com/png.latex?%5Clambda "\\lambda")
as the number of rows or columns of the square matrix
![A](https://latex.codecogs.com/png.latex?A "A"), so a 2 by 2 matrix has
two eigenvectors and two eigenvalues, a 5x5 matrix has 5 of each, etc.
One ironclad rule is that there cannot be more distinct eigenvalues than
the matrix dimension. Some matrices possess fewer eigenvalues than the
matrix dimension, those are said to have a *degenerate* set of
eigenvalues, and at least two of the eigenvectors share the same
eigenvalue.

The situation with eigenvectors is trickier. There are some matrices for
which any vector is an eigenvector, and others which have a limited set
of eigenvectors. What is difficult about counting eigenvectors is that
an eigenvector is still an eigenvector when multiplied by a constant.
You can show that for any matrix, multiplication by a constant is
commutative: $cA = Ac $, where
![A](https://latex.codecogs.com/png.latex?A "A") is a matrix and
![c](https://latex.codecogs.com/png.latex?c "c") is a constant. This
leads us to the important result that if ![\\vec
v](https://latex.codecogs.com/png.latex?%5Cvec%20v "\\vec v") is an
eigenvector with eigenvalue
![\\lambda](https://latex.codecogs.com/png.latex?%5Clambda "\\lambda"),
then any scalar multiple ![c \\vec
v](https://latex.codecogs.com/png.latex?c%20%5Cvec%20v "c \\vec v") is
also an eigenvector with the same eigenvalue. The following demonstrates
this algebraically:   
![A \\times (c \\vec v) = c A \\times \\vec v = c \\lambda \\vec v =
\\lambda (c \\vec
v)](https://latex.codecogs.com/png.latex?A%20%20%5Ctimes%20%20%28c%20%5Cvec%20v%29%20%3D%20c%20A%20%20%5Ctimes%20%20%5Cvec%20v%20%3D%20c%20%5Clambda%20%5Cvec%20v%20%3D%20%20%5Clambda%20%28c%20%5Cvec%20v%29
"A  \\times  (c \\vec v) = c A  \\times  \\vec v = c \\lambda \\vec v =  \\lambda (c \\vec v)")  
This shows that when the vector ![c \\vec
v](https://latex.codecogs.com/png.latex?c%20%5Cvec%20v "c \\vec v") is
multiplied by the matrix ![A](https://latex.codecogs.com/png.latex?A
"A"), it results in its being multiplied by the same number
![\\lambda](https://latex.codecogs.com/png.latex?%5Clambda "\\lambda"),
so by definition it is an eigenvector. Therefore, an eigenvector ![\\vec
v](https://latex.codecogs.com/png.latex?%5Cvec%20v "\\vec v") is not
unique, as any constant multiple ![c \\vec
v](https://latex.codecogs.com/png.latex?c%20%5Cvec%20v "c \\vec v") is
also an eigenvector. It is more useful to think not of a single
eigenvector ![\\vec v](https://latex.codecogs.com/png.latex?%5Cvec%20v
"\\vec v"), but of a **collection of vectors that can be interconverted
by scalar multiplication** that are all essentially the same
eigenvector. Another way to represent this, if the eigenvector is real,
is that an eigenvector as a **direction that remains unchanged by
multiplication by the matrix**, such as direction of the vector
![v](https://latex.codecogs.com/png.latex?v "v") in figure . As
mentioned above, this is true only for real eigenvalues and
eigenvectors, since complex eigenvectors cannot be used to define a
direction in a real space.

To summarize, eigenvalues and eigenvectors of a matrix are a set of
numbers and a set of vectors (up to scalar multiple) that describe the
action of the matrix as a multiplicative operator on vectors.
“Well-behaved” square ![n](https://latex.codecogs.com/png.latex?n "n")
by ![n](https://latex.codecogs.com/png.latex?n "n") matrices have
![n](https://latex.codecogs.com/png.latex?n "n") distinct eigenvalues
and ![n](https://latex.codecogs.com/png.latex?n "n") eigenvectors
pointing in distinct directions. In a deep sense, the collection of
eigenvectors and eigenvalues defines a matrix
![A](https://latex.codecogs.com/png.latex?A "A"), which is why an older
name for them is characteristic vectors and values.

### Calculating eigenvalues

Finding the eigenvalues and eigenvectors analytically, that is on paper,
is quite laborious even for 3 by 3 or 4 by 4 matrices and for larger
ones there is no analytical solution. In practice, the task is
outsourced to a computer, and MATLAB has a number of functions for this
purpose. Nevertheless, it is useful to go through the process in 2
dimensions in order to gain an understanding of what is involved. From
the definition of eigenvalues and eigenvectors, the condition can be
written in terms of the four elements of a 2 by 2 matrix:   
![
\\left(\\begin{array}{cc}a & b \\\\c &
d\\end{array}\\right)\\left(\\begin{array}{c}v\_1 \\\\ v\_2
\\end{array}\\right) = \\left(\\begin{array}{c}av\_1 +b v\_2\\\\ cv\_1+
dv\_2 \\end{array}\\right) = \\lambda \\left(\\begin{array}{c}v\_1 \\\\
v\_2 \\end{array}\\right)
](https://latex.codecogs.com/png.latex?%0A%5Cleft%28%5Cbegin%7Barray%7D%7Bcc%7Da%20%26%20b%20%5C%5Cc%20%26%20d%5Cend%7Barray%7D%5Cright%29%5Cleft%28%5Cbegin%7Barray%7D%7Bc%7Dv_1%20%5C%5C%20v_2%20%5Cend%7Barray%7D%5Cright%29%20%3D%20%5Cleft%28%5Cbegin%7Barray%7D%7Bc%7Dav_1%20%2Bb%20v_2%5C%5C%20cv_1%2B%20dv_2%20%5Cend%7Barray%7D%5Cright%29%20%3D%20%5Clambda%20%5Cleft%28%5Cbegin%7Barray%7D%7Bc%7Dv_1%20%5C%5C%20v_2%20%5Cend%7Barray%7D%5Cright%29%0A
"
\\left(\\begin{array}{cc}a & b \\\\c & d\\end{array}\\right)\\left(\\begin{array}{c}v_1 \\\\ v_2 \\end{array}\\right) = \\left(\\begin{array}{c}av_1 +b v_2\\\\ cv_1+ dv_2 \\end{array}\\right) = \\lambda \\left(\\begin{array}{c}v_1 \\\\ v_2 \\end{array}\\right)
")  
This is now a system of two linear algebraic equations, which we can
solve by substitution. First, let us solve for
![v\_1](https://latex.codecogs.com/png.latex?v_1 "v_1") in the first
row, to get   
![v\_1 =
\\frac{-bv\_2}{a-\\lambda}](https://latex.codecogs.com/png.latex?v_1%20%3D%20%5Cfrac%7B-bv_2%7D%7Ba-%5Clambda%7D
"v_1 = \\frac{-bv_2}{a-\\lambda}")  
Then we substitute this into the second equation and get:   
![
\\frac{-bcv\_2}{a-\\lambda} +(d-\\lambda)v\_2 = 0
](https://latex.codecogs.com/png.latex?%0A%5Cfrac%7B-bcv_2%7D%7Ba-%5Clambda%7D%20%2B%28d-%5Clambda%29v_2%20%3D%200%0A
"
\\frac{-bcv_2}{a-\\lambda} +(d-\\lambda)v_2 = 0
")  

Since ![v\_2](https://latex.codecogs.com/png.latex?v_2 "v_2") multiplies
both terms, and is not necessarily zero, we require that its
multiplicative factor be zero. Doing a little algebra, we obtain the
following, known as the *characteristic equation* of the matrix:   
![-bc +(a-\\lambda)(d-\\lambda) = \\lambda^2-(a+d)\\lambda +ad-bc
= 0](https://latex.codecogs.com/png.latex?-bc%20%2B%28a-%5Clambda%29%28d-%5Clambda%29%20%3D%20%5Clambda%5E2-%28a%2Bd%29%5Clambda%20%2Bad-bc%20%3D%200
"-bc +(a-\\lambda)(d-\\lambda) = \\lambda^2-(a+d)\\lambda +ad-bc = 0")  
This equation can be simplified by using two quantities we defined at
the beginning of the section: the sum of the diagonal elements called
the trace ![\\tau =
a+d](https://latex.codecogs.com/png.latex?%5Ctau%20%3D%20a%2Bd
"\\tau = a+d"), and the determinant ![\\Delta =
ad-bc](https://latex.codecogs.com/png.latex?%5CDelta%20%3D%20ad-bc
"\\Delta = ad-bc"). The quadratic equation has two solutions, dependent
solely on ![\\tau](https://latex.codecogs.com/png.latex?%5Ctau "\\tau")
and ![\\Delta](https://latex.codecogs.com/png.latex?%5CDelta "\\Delta"):
  
![\\lambda = \\frac{\\tau \\pm \\sqrt{\\tau^2-4\\Delta}}{2}
\\label{eq:2D\_eig}](https://latex.codecogs.com/png.latex?%5Clambda%20%3D%20%5Cfrac%7B%5Ctau%20%5Cpm%20%5Csqrt%7B%5Ctau%5E2-4%5CDelta%7D%7D%7B2%7D%0A%5Clabel%7Beq%3A2D_eig%7D
"\\lambda = \\frac{\\tau \\pm \\sqrt{\\tau^2-4\\Delta}}{2}
\\label{eq:2D_eig}")  
This is the general expression for a 2 by 2 matrix, showing there are
two possible eigenvalues. Note that if
![\\tau^2-4\\Delta\>0](https://latex.codecogs.com/png.latex?%5Ctau%5E2-4%5CDelta%3E0
"\\tau^2-4\\Delta\>0"), the eigenvalues are real, if
![\\tau^2-4\\Delta\<0](https://latex.codecogs.com/png.latex?%5Ctau%5E2-4%5CDelta%3C0
"\\tau^2-4\\Delta\<0"), they are complex (have real and imaginary
parts), and if
![\\tau^2-4\\Delta=0](https://latex.codecogs.com/png.latex?%5Ctau%5E2-4%5CDelta%3D0
"\\tau^2-4\\Delta=0"), there is only one eigenvalue. This situation is
known as degenerate, because two eigenvectors share the same eigenvalue.

**Example.** Let us take the same matrix we looked at in the previous
subsection:   
![A = \\left(\\begin{array}{cc}2 & 1
\\\\ 2& 3\\end{array}\\right)](https://latex.codecogs.com/png.latex?A%20%3D%20%5Cleft%28%5Cbegin%7Barray%7D%7Bcc%7D2%20%26%201%20%5C%5C%202%26%203%5Cend%7Barray%7D%5Cright%29
"A = \\left(\\begin{array}{cc}2 & 1 \\\\ 2& 3\\end{array}\\right)")  
The trace of this matrix is ![\\tau = 2+3
=5](https://latex.codecogs.com/png.latex?%5Ctau%20%3D%202%2B3%20%3D5
"\\tau = 2+3 =5") and the determinant is ![\\Delta = 6 - 2
= 4](https://latex.codecogs.com/png.latex?%5CDelta%20%3D%206%20-%202%20%3D%204
"\\Delta = 6 - 2 = 4"). Then by our formula, the eigenvalues are:   
![\\lambda = \\frac{5 \\pm \\sqrt{5^2-4 \\times 4}}{2} = \\frac{5
\\pm 3}{2}
= 4, 1](https://latex.codecogs.com/png.latex?%5Clambda%20%3D%20%5Cfrac%7B5%20%5Cpm%20%5Csqrt%7B5%5E2-4%20%5Ctimes%204%7D%7D%7B2%7D%20%20%3D%20%20%5Cfrac%7B5%20%5Cpm%203%7D%7B2%7D%20%20%3D%204%2C%201
"\\lambda = \\frac{5 \\pm \\sqrt{5^2-4 \\times 4}}{2}  =  \\frac{5 \\pm 3}{2}  = 4, 1")  
These are the multiples we found in the example above, as expected. Of
course `R` has functions to calculate this instead of doing this by
hand:

``` r
A <- matrix(c(2,2, 1, 3), nrow =2)
eigs <- eigen(A)
eigs$values
eigs$vectors
```

    # [1] 4 1
    #            [,1]       [,2]
    # [1,] -0.4472136 -0.7071068
    # [2,] -0.8944272  0.7071068

**Note:** a real-valued matrix can have complex eigenvalues and
eigenvectors, but whenever it acts on a real vector, the result is still
real. This is because the complex numbers cancel each other’s imaginary
parts.
