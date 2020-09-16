Multidimensional Scaling and Clustering
================
**Dmitry Kondrashov & Stefano Allesina**
Fundamentals of Biological Data Analysis – BIOS 26318

# Goal

Introduce Multidimensional Scaling (MDS). Given a set of “distances”
between samples, MDS attempts to arrange the samples in a
![k](https://latex.codecogs.com/png.latex?k "k")-dimensional space, such
that distances are preserved. In practice, MDS is not picky on the
notion of a distance (i.e., needs not to be a metric).

Introduce the notion of clustering, and show how a simple clustering
algorithm (k-means) work.

Let’s import some libraries:

``` r
library(reshape2) # more data massaging
library(tidyverse) # our friend the tidyverse
library(vegan) # for procrustes analysis
```

# Multidimensional Scaling (MDS)

## Mathematical approach

The input is the matrix of dissimilarities
![D](https://latex.codecogs.com/png.latex?D "D"), potentially
representing distances ![d\_{ij} = d(x\_i,
x\_j)](https://latex.codecogs.com/png.latex?d_%7Bij%7D%20%3D%20d%28x_i%2C%20x_j%29
"d_{ij} = d(x_i, x_j)"). A distance function is “metric” if:

  - ![d(x\_i, x\_j)
    \\geq 0](https://latex.codecogs.com/png.latex?d%28x_i%2C%20x_j%29%20%5Cgeq%200
    "d(x_i, x_j) \\geq 0") (non-negativity)
  - ![d(x\_i, x\_j)
    = 0](https://latex.codecogs.com/png.latex?d%28x_i%2C%20x_j%29%20%3D%200
    "d(x_i, x_j) = 0") only if ![x\_i =
    x\_j](https://latex.codecogs.com/png.latex?x_i%20%3D%20x_j
    "x_i = x_j") (identity)
  - ![d(x\_i, x\_j) = d(x\_j,
    x\_i)](https://latex.codecogs.com/png.latex?d%28x_i%2C%20x_j%29%20%3D%20d%28x_j%2C%20x_i%29
    "d(x_i, x_j) = d(x_j, x_i)") (symmetry)
  - ![d(x\_i, x\_k) \\leq d(x\_i, x\_j) + d(x\_j,
    x\_k)](https://latex.codecogs.com/png.latex?d%28x_i%2C%20x_k%29%20%5Cleq%20d%28x_i%2C%20x_j%29%20%2B%20d%28x_j%2C%20x_k%29
    "d(x_i, x_k) \\leq d(x_i, x_j) + d(x_j, x_k)") (triangle inequality)

Given a set of dissimilarities, we can therefore ask whether they are
distances, and particularly whether they represent Euclidean distances.

## Goal of MDS

Given the ![n \\times
n](https://latex.codecogs.com/png.latex?n%20%5Ctimes%20n "n \\times n")
matrix ![D](https://latex.codecogs.com/png.latex?D "D"), find a set of
coordinates ![x\_i, \\ldots x\_n \\in \\mathbb
R^p](https://latex.codecogs.com/png.latex?x_i%2C%20%5Cldots%20x_n%20%5Cin%20%5Cmathbb%20R%5Ep
"x_i, \\ldots x_n \\in \\mathbb R^p"), such that ![d\_{ij} \\approx
\\lVert x\_i - x\_j
\\rVert\_2](https://latex.codecogs.com/png.latex?d_%7Bij%7D%20%5Capprox%20%5ClVert%20x_i%20-%20x_j%20%5CrVert_2
"d_{ij} \\approx \\lVert x_i - x_j \\rVert_2") (as close as possible).
The operator ![\\lVert \\cdot
\\rVert\_2](https://latex.codecogs.com/png.latex?%5ClVert%20%5Ccdot%20%5CrVert_2
"\\lVert \\cdot \\rVert_2") is the Euclidean norm, measuring Euclidean
distance.

As such, if we can find a perfect solution, then the dissimilarities can
be mapped into Euclidean distances in a
![k](https://latex.codecogs.com/png.latex?k "k")-dimensional space.

## Classic MDS

Suppose that the elements of ![D](https://latex.codecogs.com/png.latex?D
"D") measure Euclidean distances between
![n](https://latex.codecogs.com/png.latex?n "n") points, each of which
has ![k](https://latex.codecogs.com/png.latex?k "k") coordinates:

  
![
X = \\begin{bmatrix}
x\_{11} & x\_{12} & \\dots & x\_{1k} \\\\
x\_{21} & x\_{22} & \\dots & x\_{2k} \\\\
\\vdots & \\vdots & \\ddots & \\vdots \\\\
x\_{n1} & x\_{n2} & \\dots & x\_{nk}
\\end{bmatrix}
](https://latex.codecogs.com/png.latex?%0AX%20%3D%20%5Cbegin%7Bbmatrix%7D%0A%20%20%20%20x_%7B11%7D%20%26%20x_%7B12%7D%20%26%20%20%5Cdots%20%20%26%20x_%7B1k%7D%20%5C%5C%0A%20%20%20%20x_%7B21%7D%20%26%20x_%7B22%7D%20%26%20%20%5Cdots%20%20%26%20x_%7B2k%7D%20%5C%5C%0A%20%20%20%20%5Cvdots%20%26%20%5Cvdots%20%26%20%20%5Cddots%20%26%20%5Cvdots%20%5C%5C%0A%20%20%20%20x_%7Bn1%7D%20%26%20x_%7Bn2%7D%20%26%20%20%5Cdots%20%20%26%20x_%7Bnk%7D%0A%5Cend%7Bbmatrix%7D%0A
"
X = \\begin{bmatrix}
    x_{11} & x_{12} &  \\dots  & x_{1k} \\\\
    x_{21} & x_{22} &  \\dots  & x_{2k} \\\\
    \\vdots & \\vdots &  \\ddots & \\vdots \\\\
    x_{n1} & x_{n2} &  \\dots  & x_{nk}
\\end{bmatrix}
")  
We consider the centered coordinates:

  
![
\\sum\_i x\_{ij} = 0
](https://latex.codecogs.com/png.latex?%0A%5Csum_i%20x_%7Bij%7D%20%3D%200%0A
"
\\sum_i x_{ij} = 0
")  
And the matrix ![B = X
X^t](https://latex.codecogs.com/png.latex?B%20%3D%20X%20X%5Et
"B = X X^t"), whose coefficients are ![B\_{ij} = \\sum\_k x\_{ik}
x\_{jk}](https://latex.codecogs.com/png.latex?B_%7Bij%7D%20%3D%20%5Csum_k%20x_%7Bik%7D%20x_%7Bjk%7D
"B_{ij} = \\sum_k x_{ik} x_{jk}"). We can write the square of the
distance between point ![i](https://latex.codecogs.com/png.latex?i "i")
and ![j](https://latex.codecogs.com/png.latex?j "j") as:

  
![ d\_{ij}^2 = \\sum\_k (x\_{ik} - x\_{jk})^2 = \\sum\_k x\_{ik}^2 +
\\sum\_k x\_{jk}^2 -2 \\sum\_k x\_{ik} x\_{jk} = B\_{ii} + B\_{jj} - 2
B\_{ij}](https://latex.codecogs.com/png.latex?%20d_%7Bij%7D%5E2%20%3D%20%5Csum_k%20%28x_%7Bik%7D%20-%20x_%7Bjk%7D%29%5E2%20%20%3D%20%5Csum_k%20x_%7Bik%7D%5E2%20%2B%20%5Csum_k%20x_%7Bjk%7D%5E2%20-2%20%5Csum_k%20x_%7Bik%7D%20x_%7Bjk%7D%20%3D%20B_%7Bii%7D%20%2B%20B_%7Bjj%7D%20-%202%20B_%7Bij%7D
" d_{ij}^2 = \\sum_k (x_{ik} - x_{jk})^2  = \\sum_k x_{ik}^2 + \\sum_k x_{jk}^2 -2 \\sum_k x_{ik} x_{jk} = B_{ii} + B_{jj} - 2 B_{ij}")  

Note that, because of the centering:

  
![
\\sum\_i B\_{ij} = \\sum\_i \\sum\_k x\_{ik} x\_{jk} = \\sum\_k x\_{jk}
\\sum\_i x\_{ik} = 0
](https://latex.codecogs.com/png.latex?%0A%5Csum_i%20B_%7Bij%7D%20%3D%20%5Csum_i%20%5Csum_k%20x_%7Bik%7D%20x_%7Bjk%7D%20%3D%20%5Csum_k%20x_%7Bjk%7D%20%5Csum_i%20x_%7Bik%7D%20%3D%200%0A
"
\\sum_i B_{ij} = \\sum_i \\sum_k x_{ik} x_{jk} = \\sum_k x_{jk} \\sum_i x_{ik} = 0
")  

Now we compute:

  
![
\\sum\_i d\_{ij}^2 = \\sum\_i (B\_{ii} + B\_{jj} - 2 B\_{ij}) = \\sum\_i
B\_{ii} + \\sum\_i B\_{jj} - 2 \\sum\_i B\_{ij} = \\text{Tr}(B) + n
B\_{jj} 
](https://latex.codecogs.com/png.latex?%0A%5Csum_i%20d_%7Bij%7D%5E2%20%3D%20%5Csum_i%20%28B_%7Bii%7D%20%2B%20B_%7Bjj%7D%20-%202%20B_%7Bij%7D%29%20%3D%20%5Csum_i%20B_%7Bii%7D%20%2B%20%5Csum_i%20B_%7Bjj%7D%20-%202%20%5Csum_i%20B_%7Bij%7D%20%3D%20%5Ctext%7BTr%7D%28B%29%20%2B%20n%20B_%7Bjj%7D%20%0A
"
\\sum_i d_{ij}^2 = \\sum_i (B_{ii} + B_{jj} - 2 B_{ij}) = \\sum_i B_{ii} + \\sum_i B_{jj} - 2 \\sum_i B_{ij} = \\text{Tr}(B) + n B_{jj} 
")  

Similarly (distances are symmetric):

  
![
\\sum\_j d\_{ij}^2 = \\text{Tr}(B) + n B\_{ii} 
](https://latex.codecogs.com/png.latex?%0A%5Csum_j%20d_%7Bij%7D%5E2%20%3D%20%5Ctext%7BTr%7D%28B%29%20%2B%20n%20B_%7Bii%7D%20%0A
"
\\sum_j d_{ij}^2 = \\text{Tr}(B) + n B_{ii} 
")  
And, finally:

  
![
\\sum\_i \\sum\_j d\_{ij}^2 = 2 n \\text{Tr}(B)
](https://latex.codecogs.com/png.latex?%0A%5Csum_i%20%5Csum_j%20d_%7Bij%7D%5E2%20%3D%202%20n%20%5Ctext%7BTr%7D%28B%29%0A
"
\\sum_i \\sum_j d_{ij}^2 = 2 n \\text{Tr}(B)
")  

From these three equations, we obtain:

  
![
B\_{ii} = \\frac{\\sum\_j d\_{ij}^2}{n} - \\frac{\\sum\_i \\sum\_j
d\_{ij}^2 }{2 n^2}
](https://latex.codecogs.com/png.latex?%0AB_%7Bii%7D%20%3D%20%5Cfrac%7B%5Csum_j%20d_%7Bij%7D%5E2%7D%7Bn%7D%20-%20%5Cfrac%7B%5Csum_i%20%5Csum_j%20d_%7Bij%7D%5E2%20%7D%7B2%20n%5E2%7D%0A
"
B_{ii} = \\frac{\\sum_j d_{ij}^2}{n} - \\frac{\\sum_i \\sum_j d_{ij}^2 }{2 n^2}
")  
and

  
![
B\_{jj} = \\frac{\\sum\_i d\_{ij}^2}{n} - \\frac{\\sum\_i \\sum\_j
d\_{ij}^2 }{2 n^2}
](https://latex.codecogs.com/png.latex?%0AB_%7Bjj%7D%20%3D%20%5Cfrac%7B%5Csum_i%20d_%7Bij%7D%5E2%7D%7Bn%7D%20-%20%5Cfrac%7B%5Csum_i%20%5Csum_j%20d_%7Bij%7D%5E2%20%7D%7B2%20n%5E2%7D%0A
"
B_{jj} = \\frac{\\sum_i d_{ij}^2}{n} - \\frac{\\sum_i \\sum_j d_{ij}^2 }{2 n^2}
")  

Therefore:

  
![ 
B\_{ij} = -\\frac{1}{2}(d\_{ij}^2 - B\_{ii} - B\_{jj}) =
-\\frac{1}{2}\\left(d\_{ij}^2 - \\frac{\\sum\_i d\_{ij}^2}{n} -
\\frac{\\sum\_j d\_{ij}^2}{n} + \\frac{\\sum\_i \\sum\_j d\_{ij}^2
}{n^2} \\right)
](https://latex.codecogs.com/png.latex?%20%0AB_%7Bij%7D%20%3D%20-%5Cfrac%7B1%7D%7B2%7D%28d_%7Bij%7D%5E2%20-%20B_%7Bii%7D%20-%20B_%7Bjj%7D%29%20%3D%20-%5Cfrac%7B1%7D%7B2%7D%5Cleft%28d_%7Bij%7D%5E2%20-%20%5Cfrac%7B%5Csum_i%20d_%7Bij%7D%5E2%7D%7Bn%7D%20-%20%5Cfrac%7B%5Csum_j%20d_%7Bij%7D%5E2%7D%7Bn%7D%20%20%2B%20%5Cfrac%7B%5Csum_i%20%5Csum_j%20d_%7Bij%7D%5E2%20%7D%7Bn%5E2%7D%20%5Cright%29%0A
" 
B_{ij} = -\\frac{1}{2}(d_{ij}^2 - B_{ii} - B_{jj}) = -\\frac{1}{2}\\left(d_{ij}^2 - \\frac{\\sum_i d_{ij}^2}{n} - \\frac{\\sum_j d_{ij}^2}{n}  + \\frac{\\sum_i \\sum_j d_{ij}^2 }{n^2} \\right)
")  

With some algebra, one can show that this is equivalent to:

  
![B = -\\frac{1}{2} C D^{(2)}
C](https://latex.codecogs.com/png.latex?B%20%3D%20-%5Cfrac%7B1%7D%7B2%7D%20C%20D%5E%7B%282%29%7D%20C
"B = -\\frac{1}{2} C D^{(2)} C")  

Where ![D^{(2)}](https://latex.codecogs.com/png.latex?D%5E%7B%282%29%7D
"D^{(2)}") is the matrix of squared distances, and
![C](https://latex.codecogs.com/png.latex?C "C") is the centering matrix
![C = 1 - \\frac{1}{n}\\mathcal
O](https://latex.codecogs.com/png.latex?C%20%3D%201%20-%20%5Cfrac%7B1%7D%7Bn%7D%5Cmathcal%20O
"C = 1 - \\frac{1}{n}\\mathcal O") (and ![\\mathcal
O](https://latex.codecogs.com/png.latex?%5Cmathcal%20O "\\mathcal O") is
the matrix of all ones). Thus, we can obtain
![B](https://latex.codecogs.com/png.latex?B "B") directly from the
distance matrix. Once we’ve done this,
![X](https://latex.codecogs.com/png.latex?X "X") can be found by taking
the eigenvalue decomposition:

  
![
B = X X^t = Q \\Lambda Q^t
](https://latex.codecogs.com/png.latex?%0AB%20%3D%20X%20X%5Et%20%3D%20Q%20%5CLambda%20Q%5Et%0A
"
B = X X^t = Q \\Lambda Q^t
")  

(where ![Q](https://latex.codecogs.com/png.latex?Q "Q") is the matrix of
eigenvectors of ![B](https://latex.codecogs.com/png.latex?B "B"), and
![\\Lambda](https://latex.codecogs.com/png.latex?%5CLambda "\\Lambda") a
diagonal matrix of the eigenvalues of
![B](https://latex.codecogs.com/png.latex?B "B")). Therefore:

  
![ X = Q
\\Lambda^{\\frac{1}{2}}](https://latex.codecogs.com/png.latex?%20X%20%3D%20Q%20%5CLambda%5E%7B%5Cfrac%7B1%7D%7B2%7D%7D
" X = Q \\Lambda^{\\frac{1}{2}}")  

## Reconstructing the map of Chicago

To test MDS, I have built a matrix expressing the distances between the
604 Divvy bikes station in Chicago.

The distances are measured in degrees of latitude/longitude. We
introduce some noise to see how robust our estimate is.

``` r
# load distance matrix 
load("data/divvy_stations_distances.RData")
# add some noise to make it more fun
# we change each distance of +/- 20%
n <- nrow(distance_matrix)
distance_matrix <- distance_matrix * 
  matrix(runif(n * n, min = 0.8, max = 1.2), n, n)
# plot distances
distance_matrix %>% reshape2::melt() %>% 
  ggplot(aes(x = Var1, y = Var2, fill = value)) + geom_tile()
```

<img src="multidimensional_scaling_files/figure-gfm/unnamed-chunk-2-1.png" style="display: block; margin: auto;" />

Now we perform classical MDS, and plot the coordinates we’ve recovered:

``` r
# classical MDS
mds_fit <- cmdscale(distance_matrix, k = 2) # k is the dimension of the embedding
mds_fit <- tibble(id = rownames(distance_matrix), 
                  longitude = mds_fit[,1], latitude = mds_fit[,2])
mds_fit %>% ggplot() + aes(x = latitude, y = longitude) + geom_point()
```

<img src="multidimensional_scaling_files/figure-gfm/unnamed-chunk-3-1.png" style="display: block; margin: auto;" />

And now let’s do it the “hard way”:

``` r
D2 <- distance_matrix^2
CM <- diag(rep(1, n)) - 1 / n
B <- -(1/2) * CM %*% D2 %*% CM
eB <- eigen(B)
mds_hard <- Re(eB$vectors) %*% diag(sqrt(Re(abs(eB$values))))
ggplot(data = tibble(latitude = mds_hard[,2], longitude = mds_hard[,1])) +  
  aes(x = latitude, y = longitude) + 
  geom_point()
```

<img src="multidimensional_scaling_files/figure-gfm/unnamed-chunk-4-1.png" style="display: block; margin: auto;" />

As we were saying, distances are invariant to rotation, translation and
reflection. As such, the coordinates can be rotated and centered
arbitrarily. To find the best matching rotation, we use [procrustes
analysis](https://en.wikipedia.org/wiki/Procrustes_analysis):

``` r
# this is the true location of the stations
actual_locations <- read_csv("data/divvy_stations.csv")
# aligning coordinates using Procrustes rotation and centering
procr <- procrustes(actual_locations %>% 
                      select(latitude, longitude),
                    mds_fit %>% 
                      select(latitude, longitude), 
                    scale = FALSE)
new_coord <- t(t(as.matrix(mds_fit %>% 
                             select(latitude, longitude)) %*% procr$rotation) + procr$translation[1,])
ggplot(actual_locations) + aes(x = longitude, y = latitude) + geom_point() + 
  geom_point(aes(x = new_coord[,2], y = new_coord[,1]), colour = "red", shape = 1)
```

<img src="multidimensional_scaling_files/figure-gfm/unnamed-chunk-5-1.png" style="display: block; margin: auto;" />

You can see that we’ve been doing really well. The largest discrepancies
are around the borders, where we have less information. To show how much
should we shift the points to recover the actual data, you can
use

``` r
plot(procr)
```

<img src="multidimensional_scaling_files/figure-gfm/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

# Clustering

## Goal of clustering

Given a set of ![n](https://latex.codecogs.com/png.latex?n "n")
observations, divide them into
![k](https://latex.codecogs.com/png.latex?k "k") disjoint sets ![S =
\\{S\_1, S\_2, \\ldots,
S\_k\\}](https://latex.codecogs.com/png.latex?S%20%3D%20%5C%7BS_1%2C%20S_2%2C%20%5Cldots%2C%20S_k%5C%7D
"S = \\{S_1, S_2, \\ldots, S_k\\}") such that observations in the same
set are “similar” to each other and “differ” from observations in other
sets. The notion of similarity is context-dependent. Clustering (or
classification) is a very active area of research, with hundreds of
algorithms and techniques developed for different problems.

## K-means

This is a simple technique to cluster
![n](https://latex.codecogs.com/png.latex?n "n") observations into
![k](https://latex.codecogs.com/png.latex?k "k") clusters. The number of
clusters ![k](https://latex.codecogs.com/png.latex?k "k") has to be
specified. The input is a ![d-](https://latex.codecogs.com/png.latex?d-
"d-")dimensional vector ![(\\mathbf x\_1, \\mathbf x\_2, \\ldots,
\\mathbf
x\_n)](https://latex.codecogs.com/png.latex?%28%5Cmathbf%20x_1%2C%20%5Cmathbf%20x_2%2C%20%5Cldots%2C%20%5Cmathbf%20x_n%29
"(\\mathbf x_1, \\mathbf x_2, \\ldots, \\mathbf x_n)"), where ![\\mathbf
x\_i \\in \\mathbb
R^d](https://latex.codecogs.com/png.latex?%5Cmathbf%20x_i%20%5Cin%20%5Cmathbb%20R%5Ed
"\\mathbf x_i \\in \\mathbb R^d"), and a number of clusters
![k](https://latex.codecogs.com/png.latex?k "k"). We want to assign each
vector ![\\matbb
x\_i](https://latex.codecogs.com/png.latex?%5Cmatbb%20x_i "\\matbb x_i")
to one of ![k](https://latex.codecogs.com/png.latex?k "k") sets
![S\_j](https://latex.codecogs.com/png.latex?S_j "S_j"), so that the
vectors (points in ![\\mathbb
R^d](https://latex.codecogs.com/png.latex?%5Cmathbb%20R%5Ed
"\\mathbb R^d")) are closest to each other. Mathematically, we want to
find the assignments to sets such that the variance within each set is
minimized:

  
![
\\underset{S}{\\mathrm{arg\\, min}} \\sum\_{i=i}^{k}\\sum\_{\\mathbf{x}
\\in S\_i} \\Vert \\mathbf{x} - \\mathbf{\\mu}\_i \\Vert^2 
](https://latex.codecogs.com/png.latex?%0A%5Cunderset%7BS%7D%7B%5Cmathrm%7Barg%5C%2C%20min%7D%7D%20%5Csum_%7Bi%3Di%7D%5E%7Bk%7D%5Csum_%7B%5Cmathbf%7Bx%7D%20%5Cin%20S_i%7D%20%5CVert%20%5Cmathbf%7Bx%7D%20-%20%5Cmathbf%7B%5Cmu%7D_i%20%5CVert%5E2%20%0A
"
\\underset{S}{\\mathrm{arg\\, min}} \\sum_{i=i}^{k}\\sum_{\\mathbf{x} \\in S_i} \\Vert \\mathbf{x} - \\mathbf{\\mu}_i \\Vert^2 
")  
where ![\\mathbf
\\mu\_i](https://latex.codecogs.com/png.latex?%5Cmathbf%20%5Cmu_i
"\\mathbf \\mu_i") is the vector of the means of points in set
![i](https://latex.codecogs.com/png.latex?i "i") (one mean for each
dimension). Equivalently:

  
![
\\underset{S}{\\mathrm{arg\\, min}} \\vert S\_i \\vert \\text{Var}
(S\_i)
](https://latex.codecogs.com/png.latex?%0A%5Cunderset%7BS%7D%7B%5Cmathrm%7Barg%5C%2C%20min%7D%7D%20%5Cvert%20S_i%20%5Cvert%20%5Ctext%7BVar%7D%20%28S_i%29%0A
"
\\underset{S}{\\mathrm{arg\\, min}} \\vert S_i \\vert \\text{Var} (S_i)
")  

where ![\\vert S\_i
\\vert](https://latex.codecogs.com/png.latex?%5Cvert%20S_i%20%5Cvert
"\\vert S_i \\vert") is the size of set
![i](https://latex.codecogs.com/png.latex?i "i") and ![\\text{Var}
(S\_i)](https://latex.codecogs.com/png.latex?%5Ctext%7BVar%7D%20%28S_i%29
"\\text{Var} (S_i)") is the variance of the coordinates in set
![i](https://latex.codecogs.com/png.latex?i "i").

It can be proved that this problem is in general NP-hard (i.e., finding
the best solution requires a number of operations that scales
non-polynomially with the size of the data). Fortunately, there are good
heuristics, and several `R` packages that implement different
algorithms.

The simplest algorithm starts by choosing at random (or more smartly)
![k](https://latex.codecogs.com/png.latex?k "k") observations to serve
as initial “means” of the clusters; then:

1.  Each point is assigned to the cluster with the closest mean
2.  Means are re-computed
3.  Go back to 1 until the algorithm has converged (i.e., membership has
    not changed)

Note that the quality of the solution depends on which points are used
to intialize the means; as such, if these are chosen at random, the
algorithm will possibly converge to different solutions when different
initial conditions are chosen.

## Example

We are going to use a simple 2-d synthetic data set taken from
[here](http://cs.joensuu.fi/sipu/datasets/)

``` r
dt <- read_delim("data/s1.txt", 
                 col_names = c("x", "y"), 
                 delim = " ", 
                 trim_ws = TRUE)
ggplot(dt) + aes(x = x, y = y) + geom_point()
```

<img src="multidimensional_scaling_files/figure-gfm/unnamed-chunk-7-1.png" style="display: block; margin: auto;" />

There are 15 clusters in the data. Now we’re going to try to detect them
using k-means

``` r
# with three clusters
cl <- as.factor(kmeans(dt, centers = 3)$cluster)
ggplot(dt %>% add_column(cl = cl)) + 
  aes(x = x, y = y, colour = cl) + 
  geom_point()
```

<img src="multidimensional_scaling_files/figure-gfm/unnamed-chunk-8-1.png" style="display: block; margin: auto;" />

And now with 15 clusters

``` r
# with three clusters
set.seed(1) #set a random starting point
cl <- as.factor(kmeans(dt, centers = 15)$cluster)
ggplot(dt %>% add_column(cl = cl)) + 
  aes(x = x, y = y, colour = cl) + 
  geom_point()
```

<img src="multidimensional_scaling_files/figure-gfm/unnamed-chunk-9-1.png" style="display: block; margin: auto;" />

**Exercise:** Now try again but choosing a different seed: do you get
the same result? Write code that repeats the same analysis several times
and takes the best result (note that `kmeans(...)$tot.withinss` is the
total sum of squares within the clusters—i.e., the quantity you want to
minimize).

### Exercise: PCA sommelier

The file `Wine.csv` contains several measures made on 178 wines from
Piedmont, produced using three different grapes (column `Grape`, with 1
= Barolo, 2 = Grignolino, 3 = Barbera). Use the 13 measured variables
(i.e., all but `Grape`) to perform a PCA, and use k-means to cluster the
data in PC space. Can you recover the right classification of grapes?

``` r
dt <- read_csv("data/Wine.csv")
# make into a matrix for PCA
mat <- dt %>% select(-Grape) %>% as.matrix()
# perform PCA by scaling and centering all 13 variables
pca <- prcomp(mat, center = TRUE, scale. = TRUE)
# pca$x contains the coordinates
# let's see how are we doing
ggplot(cbind(dt, as_tibble(pca$x))) + 
  aes(x = PC1, y = PC2, colour = factor(Grape)) + 
  geom_point()
```

<img src="multidimensional_scaling_files/figure-gfm/unnamed-chunk-10-1.png" style="display: block; margin: auto;" />

Now let’s apply k-means to divide the data into three clusters

``` r
cl <- kmeans(pca$x, centers = 3)
ggplot(dt %>% add_column(cluster = cl$cluster)) + 
  aes(x = Grape, y = cluster, colour = factor(Grape)) + 
  geom_jitter()
```

<img src="multidimensional_scaling_files/figure-gfm/unnamed-chunk-11-1.png" style="display: block; margin: auto;" />

You can see that we correctly classify all the wines for grapes 1 and 3,
while some of the ones from grape 2 are misclassified.

``` r
# find convex hull for each cluster
# i.e. minimal convex polygon containing all points
hull <- pca$x %>% as_tibble() %>% 
  select(PC1, PC2) %>% 
  add_column(cluster = cl$cluster) %>% 
  group_by(cluster) %>% 
  slice(chull(PC1, PC2))

# plot the convex hulls
pl <- hull %>% ggplot(aes(x = PC1, y = PC2, group = cluster, fill = factor(cluster))) + geom_polygon(alpha = 0.5)
  
show(pl)
```

<img src="multidimensional_scaling_files/figure-gfm/unnamed-chunk-12-1.png" style="display: block; margin: auto;" />

``` r
# now add the points
pl + geom_point(data = pca$x %>% as_tibble() %>% 
         select(PC1, PC2) %>% 
         add_column(Grape = dt$Grape, cluster = cl$cluster), 
  aes(x = PC1, y = PC2, colour = factor(Grape), shape = factor(Grape)))
```

<img src="multidimensional_scaling_files/figure-gfm/unnamed-chunk-12-2.png" style="display: block; margin: auto;" />

You can see that we are misclassifying the grapes that are at the border
of the cluster or nestled among the points belonging to a different
grape. But the outcome is quite good: we would have misclassified only 6
wines out of 178 (3.3%).
