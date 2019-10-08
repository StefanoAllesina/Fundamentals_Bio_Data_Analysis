Distributions and their properties
================
**Dmitry Kondrashov & Stefano Allesina**
Fundamentals of Biological Data Analysis -- BIOS 26318

-   [Objectives:](#objectives)
-   [Independence](#independence)
    -   [Conditional probability](#conditional-probability)
    -   [Independence](#independence-1)
    -   [Usefulness of independence](#usefulness-of-independence)
-   [Probability distribution examples (discrete)](#probability-distribution-examples-discrete)
    -   [Uniform](#uniform)
    -   [Binomial](#binomial)
    -   [Geometric](#geometric)
    -   [Poisson](#poisson)
-   [Probability distribution examples (continuous)](#probability-distribution-examples-continuous)
    -   [Uniform](#uniform-1)
    -   [Exponential](#exponential)
    -   [Normal distribution](#normal-distribution)
-   [Application of normal distribution: confidence intervals:](#application-of-normal-distribution-confidence-intervals)

Objectives:
===========

-   Apply concepts of conditional probability to practical scenarios and questions
-   Describe independence as a concept and apply to data sets
-   Use random number generators to simulate various distributions
-   Be familiar with the shape of several common distributions and describe the role of their parameters

Independence
============

Conditional probability
-----------------------

In the basic definitions of probability, we considered the probabilities of each outcome and events separately. Let us consider how information about one event affects the probability of another event. The concept is that if one event (let's call it ![B](https://latex.codecogs.com/png.latex?B "B")) is true, unless the event is the entire space, it rules out some other outcomes. This may affect the probability of other events (e.g., ![A](https://latex.codecogs.com/png.latex?A "A")) in the sample space, because knowledge of ![B](https://latex.codecogs.com/png.latex?B "B") may rule out some of the outcomes in ![A](https://latex.codecogs.com/png.latex?A "A") as well. Here is the formal definition:

> **Definition**: For two events ![A](https://latex.codecogs.com/png.latex?A "A") and ![B](https://latex.codecogs.com/png.latex?B "B") in a sample space ![\\Omega](https://latex.codecogs.com/png.latex?%5COmega "\Omega") with a probability measure ![P](https://latex.codecogs.com/png.latex?P "P"), the probability of ![A](https://latex.codecogs.com/png.latex?A "A") given ![B](https://latex.codecogs.com/png.latex?B "B"), called the **conditional probability**, defined as:

![P(A \\vert B) = \\frac{P(A \\cap B)}{P(B)}](https://latex.codecogs.com/png.latex?P%28A%20%5Cvert%20B%29%20%3D%20%5Cfrac%7BP%28A%20%5Ccap%20B%29%7D%7BP%28B%29%7D "P(A \vert B) = \frac{P(A \cap B)}{P(B)}")

> where ![A \\cap B](https://latex.codecogs.com/png.latex?A%20%5Ccap%20B "A \cap B") or ![A, B](https://latex.codecogs.com/png.latex?A%2C%20B "A, B") is the intersection of events ![A](https://latex.codecogs.com/png.latex?A "A") and ![B](https://latex.codecogs.com/png.latex?B "B"), also known as "![A](https://latex.codecogs.com/png.latex?A "A") and ![B](https://latex.codecogs.com/png.latex?B "B")"---the event consisting of all outcomes that are in both ![A](https://latex.codecogs.com/png.latex?A "A") and ![B](https://latex.codecogs.com/png.latex?B "B").

In words, given the knowledge that an event ![B](https://latex.codecogs.com/png.latex?B "B") occurs, the sample space is restricted to the subset ![B](https://latex.codecogs.com/png.latex?B "B"), which is why the denominator in the definition is ![P(B)](https://latex.codecogs.com/png.latex?P%28B%29 "P(B)"). The numerator encompasses all the outcomes we are interested in, (i.e., ![A](https://latex.codecogs.com/png.latex?A "A")), but since we are now restricted to ![B](https://latex.codecogs.com/png.latex?B "B"), the numerator consists of all the outcomes of ![A](https://latex.codecogs.com/png.latex?A "A") which are also in ![B](https://latex.codecogs.com/png.latex?B "B"), or ![A \\cap B](https://latex.codecogs.com/png.latex?A%20%5Ccap%20B "A \cap B"). The definition makes sense in two extreme cases: if ![A = B](https://latex.codecogs.com/png.latex?A%20%3D%20B "A = B") and if ![A](https://latex.codecogs.com/png.latex?A "A") and ![B](https://latex.codecogs.com/png.latex?B "B") are mutually exclusive:

![P(B\\vert B) = P(B \\cap B) /P(B) = P(B)/P(B) = 1](https://latex.codecogs.com/png.latex?P%28B%5Cvert%20B%29%20%3D%20P%28B%20%5Ccap%20B%29%20%2FP%28B%29%20%3D%20P%28B%29%2FP%28B%29%20%3D%201 "P(B\vert B) = P(B \cap B) /P(B) = P(B)/P(B) = 1")

If ![P(A\\cap B) = 0](https://latex.codecogs.com/png.latex?P%28A%5Ccap%20B%29%20%3D%200 "P(A\cap B) = 0"), then ![P(A\\vert B) = 0/P(B) = 0](https://latex.codecogs.com/png.latex?P%28A%5Cvert%20B%29%20%3D%200%2FP%28B%29%20%3D%200 "P(A\vert B) = 0/P(B) = 0")

**Important note:** one common source of confusion about conditional probability is the difference between the probability of ![A](https://latex.codecogs.com/png.latex?A "A") and ![B](https://latex.codecogs.com/png.latex?B "B") and the probability of ![A](https://latex.codecogs.com/png.latex?A "A") given ![B](https://latex.codecogs.com/png.latex?B "B"). This is a result of the discrepancy between everyday word usage and mathematical terminology, because the statement "what are the odds of finding a tall person who also likes tea?" is hard to distinguish from "what are the odds that a person who is tall likes tea?" The critical difference between these two statements is that in the former you start out with no information and are picking out a person from the entire population, while is in the latter you start out with the knowledge that a person is tall.

**Example:** In the classic Mendelian pea experiment, each diploid organism carries two alleles. The allele ![A](https://latex.codecogs.com/png.latex?A "A") is dominant and results in pink flowers, while ![a](https://latex.codecogs.com/png.latex?a "a") is recessive and results in white flowers. There are three possible genotypes (![AA](https://latex.codecogs.com/png.latex?AA "AA"), ![Aa](https://latex.codecogs.com/png.latex?Aa "Aa"), ![aa](https://latex.codecogs.com/png.latex?aa "aa")) and two phenotypes (Pink or White). For the questions below, assume that two heterozygous pea plants (each having genotype ![Aa](https://latex.codecogs.com/png.latex?Aa "Aa")) are crossed, producing the following table of genotypes with equal probabilites in each cell:

| parent | A         | a          |
|--------|-----------|------------|
| A      | AA (pink) | Aa (pink)  |
| a      | Aa (pink) | aa (white) |

1.  What is the probability that a plant with pink flowers has genotype ![AA](https://latex.codecogs.com/png.latex?AA "AA")? Write this down in terms of conditional probability and explain how it's different from the probability of a plant having both pink flower and genotype ![AA](https://latex.codecogs.com/png.latex?AA "AA").

2.  What is the probability that a plant with genotype ![AA](https://latex.codecogs.com/png.latex?AA "AA") has pink flowers? Again, write down the conditional probability and explain how it's different from the probability of a plant having both pink flower and genotype ![AA](https://latex.codecogs.com/png.latex?AA "AA").

**Lesson:** in general,

![P(X \\vert  Y ) \\neq P(Y \\vert  X)](https://latex.codecogs.com/png.latex?P%28X%20%5Cvert%20%20Y%20%29%20%5Cneq%20P%28Y%20%5Cvert%20%20X%29 "P(X \vert  Y ) \neq P(Y \vert  X)")

Independence
------------

Independence is a fundamental concept in probability that may be misinterpreted without careful thinking. Intuitively, two events (or random variables) are independent if one does not influence the other. More precisely, it means that the probability of one event is the same regarless of whether the other one happens or not. This is expressed precisely using conditional probabilities:

> **Definition**: Two events ![A](https://latex.codecogs.com/png.latex?A "A") and ![B](https://latex.codecogs.com/png.latex?B "B") in a sample space ![\\Omega](https://latex.codecogs.com/png.latex?%5COmega "\Omega") with a probability measure ![P](https://latex.codecogs.com/png.latex?P "P") are **independent** if ![P(A\\vert B) = P(A)](https://latex.codecogs.com/png.latex?P%28A%5Cvert%20B%29%20%3D%20P%28A%29 "P(A\vert B) = P(A)"), or equivalently if ![P(B\\vert A) = P(B)](https://latex.codecogs.com/png.latex?P%28B%5Cvert%20A%29%20%3D%20P%28B%29 "P(B\vert A) = P(B)").

Independence is not a straightforward concept. It may be confused with mutual exclusivity, as one might surmise that if ![A](https://latex.codecogs.com/png.latex?A "A") and ![B](https://latex.codecogs.com/png.latex?B "B") have no overlap, then they are independent. That however, is false by definition, since ![P(A\\vert B)](https://latex.codecogs.com/png.latex?P%28A%5Cvert%20B%29 "P(A\vert B)") is 0 for two mutually exclusive events. The confusion stems from thinking that if ![A](https://latex.codecogs.com/png.latex?A "A") and ![B](https://latex.codecogs.com/png.latex?B "B") are non-overlapping, then they do not influence each other. But the notion of influence in this definition is about information; so if ![A](https://latex.codecogs.com/png.latex?A "A") and ![B](https://latex.codecogs.com/png.latex?B "B") are mutually exclusive, the knowledge that one of them occurs has an influence of the probability of the other one occurring, specifically it rules the other one out.

**Example:** In the sample space of weather phenomena, are the events of snowing and hot weather independent?

**Example:** A slighly more subtle example, the lifetime risk of breast cancer is about 1 in 8 for women and about 1 in 1000 for men. Are sex and breast cancer independent?

Usefulness of independence
--------------------------

Independence is a mathematical abstraction, and reality rarely provides us with perfectly independent variables. But it's a very useful abstraction in that it enables calculations that would be difficult or impossible to carry out without this assumption.

First, independence allows for calculating the probability of two events or two random variables simultaneously. This is a straightforward consequence of the definition conditional probability (first equality) and independence (second equality):

![\\frac{P(A \\cap B)}{P(B)}= P(A\\vert B) = P(A)](https://latex.codecogs.com/png.latex?%5Cfrac%7BP%28A%20%5Ccap%20B%29%7D%7BP%28B%29%7D%3D%20P%28A%5Cvert%20B%29%20%3D%20P%28A%29 "\frac{P(A \cap B)}{P(B)}= P(A\vert B) = P(A)")

 Multiplying both sides by ![P(B)](https://latex.codecogs.com/png.latex?P%28B%29 "P(B)"), we get the **product rule** of independence, perhaps the most widely used formula in applied probability:

![P(A \\cap B) = P(A)P(B)](https://latex.codecogs.com/png.latex?P%28A%20%5Ccap%20B%29%20%3D%20P%28A%29P%28B%29 "P(A \cap B) = P(A)P(B)")

**Example:** The probability that two randomly selected individuals have red hair--assuming that the occurence of this trait is independent--is the square of the probability of red hair in one individual. (Note that this is never exactly the case for a finite population---why?)

**Example:** The probability of two alleles of two separate genes (call them A and B) occurring on the same gamete may be independent or may be linked. In population genetics, the concept of *linkage disequilibrium* describes the extent of such linkage; for example, alleles that are located on separate chromosomes (in eukaryotes) are usually not linked and their occurrence is independent. The *coefficient of linkage disequilibrium* is defined as the difference between what is expected from independence and the actual probability of both alleles being present:

![D\_{AB} = P(A \\cap B) - P(A)P(B) ](https://latex.codecogs.com/png.latex?D_%7BAB%7D%20%3D%20P%28A%20%5Ccap%20B%29%20-%20P%28A%29P%28B%29%20 "D_{AB} = P(A \cap B) - P(A)P(B) ")

 ![P(A)](https://latex.codecogs.com/png.latex?P%28A%29 "P(A)") and ![P(B)](https://latex.codecogs.com/png.latex?P%28B%29 "P(B)") are the frequencies of the two respective alleles (haplotypes) in the population, while ![P(A \\cap B)](https://latex.codecogs.com/png.latex?P%28A%20%5Ccap%20B%29 "P(A \cap B)") is the frequency of the haplotypes occurring together in the same copy of the genome (that is, on the same gamete). For two independent loci, ![D\_{AB} = 0](https://latex.codecogs.com/png.latex?D_%7BAB%7D%20%3D%200 "D_{AB} = 0"), while for loci that usually occur together the coefficient will be positive, and its magnitude is influenced both by physical proximity of the loci on a chromosome, the evolutionary history of the species, and other factors.

Another important consequence of independence has to do with the sum of two independent random variables. The expectation of the sum of any random variables is linear, which can be demonstrated using some work with sums, starting from the definition of expectation (the same can be shown for continuous random variables, using integrals instead of sums):

![E(X + Y) = \\sum\_i \\sum\_j (x\_i + y\_j) P(x\_i, y\_j) =](https://latex.codecogs.com/png.latex?E%28X%20%2B%20Y%29%20%3D%20%5Csum_i%20%5Csum_j%20%28x_i%20%2B%20y_j%29%20P%28x_i%2C%20y_j%29%20%3D "E(X + Y) = \sum_i \sum_j (x_i + y_j) P(x_i, y_j) =")

![= \\sum\_i \\sum\_j x\_iP(x\_i, y\_j) + \\sum\_i \\sum\_j y\_j P(x\_i, y\_j) = \\sum\_i x\_i \\sum\_j P(x\_i, y\_j) + \\sum\_j y\_j \\sum\_i  P(x\_i, y\_j) = ](https://latex.codecogs.com/png.latex?%3D%20%5Csum_i%20%5Csum_j%20x_iP%28x_i%2C%20y_j%29%20%2B%20%5Csum_i%20%5Csum_j%20y_j%20P%28x_i%2C%20y_j%29%20%3D%20%5Csum_i%20x_i%20%5Csum_j%20P%28x_i%2C%20y_j%29%20%2B%20%5Csum_j%20y_j%20%5Csum_i%20%20P%28x_i%2C%20y_j%29%20%3D%20 "= \sum_i \sum_j x_iP(x_i, y_j) + \sum_i \sum_j y_j P(x_i, y_j) = \sum_i x_i \sum_j P(x_i, y_j) + \sum_j y_j \sum_i  P(x_i, y_j) = ")

 Summing up a joint probability distribution over all values of one variable removes that variable, ![\\sum\_j P(x\_i, y\_j) = P(x\_i)](https://latex.codecogs.com/png.latex?%5Csum_j%20P%28x_i%2C%20y_j%29%20%3D%20P%28x_i%29 "\sum_j P(x_i, y_j) = P(x_i)") ![\\sum\_i P(x\_i, y\_j) = P(y\_j)](https://latex.codecogs.com/png.latex?%5Csum_i%20P%28x_i%2C%20y_j%29%20%3D%20P%28y_j%29 "\sum_i P(x_i, y_j) = P(y_j)"), so this leave us with the two separate expected values:

![= \\sum\_i x\_i P(x\_i) + \\sum\_j y\_j P(y\_j) = E(X) + E(Y)](https://latex.codecogs.com/png.latex?%3D%20%5Csum_i%20x_i%20P%28x_i%29%20%2B%20%5Csum_j%20y_j%20P%28y_j%29%20%3D%20E%28X%29%20%2B%20E%28Y%29 "= \sum_i x_i P(x_i) + \sum_j y_j P(y_j) = E(X) + E(Y)")

 However, this is not the case for the variance in general (using ![E\_X](https://latex.codecogs.com/png.latex?E_X "E_X") and ![E\_Y](https://latex.codecogs.com/png.latex?E_Y "E_Y") to indicate the expected values of ![X](https://latex.codecogs.com/png.latex?X "X") and ![Y](https://latex.codecogs.com/png.latex?Y "Y") to reduce the number of parentheses):

![\\text{Var}(X+Y) = E \\left\[ (X+Y)-(E\_X+E\_Y) \\right\]^2 = ](https://latex.codecogs.com/png.latex?%5Ctext%7BVar%7D%28X%2BY%29%20%3D%20E%20%5Cleft%5B%20%28X%2BY%29-%28E_X%2BE_Y%29%20%5Cright%5D%5E2%20%3D%20 "\text{Var}(X+Y) = E \left[ (X+Y)-(E_X+E_Y) \right]^2 = ")

![=E\[ (X-E\_X)^2 +(Y-E\_Y)^2 - 2(X-E\_X)(Y-E\_Y)\] =  ](https://latex.codecogs.com/png.latex?%3DE%5B%20%28X-E_X%29%5E2%20%2B%28Y-E_Y%29%5E2%20-%202%28X-E_X%29%28Y-E_Y%29%5D%20%3D%20%20 "=E[ (X-E_X)^2 +(Y-E_Y)^2 - 2(X-E_X)(Y-E_Y)] =  ")

![=E (X-E\_X)^2 +  E(Y-E\_Y)^2 - 2 E\[(X-E\_X)(Y-E\_Y)\]  = ](https://latex.codecogs.com/png.latex?%3DE%20%28X-E_X%29%5E2%20%2B%20%20E%28Y-E_Y%29%5E2%20-%202%20E%5B%28X-E_X%29%28Y-E_Y%29%5D%20%20%3D%20 "=E (X-E_X)^2 +  E(Y-E_Y)^2 - 2 E[(X-E_X)(Y-E_Y)]  = ")

 The first two terms are the respective variances, while the third term is called the *covariance* of ![X](https://latex.codecogs.com/png.latex?X "X") and ![Y](https://latex.codecogs.com/png.latex?Y "Y"):

![= \\text{Var}(X) + \\text{Var}(Y) - 2 \\text{Cov}(X,Y) ](https://latex.codecogs.com/png.latex?%3D%20%5Ctext%7BVar%7D%28X%29%20%2B%20%5Ctext%7BVar%7D%28Y%29%20-%202%20%5Ctext%7BCov%7D%28X%2CY%29%20 "= \text{Var}(X) + \text{Var}(Y) - 2 \text{Cov}(X,Y) ")

 Covariance describes how much two random variables vary together, or more precisely, how much they deviate from their respective means in the same direction. Thus it should be reasonable to think that two indepdendent random variables have covariance 0, which is demonstrated as follows:

![E\[(X-E\_X)(Y-E\_Y)\] = E(XY) - E\_Y E\_X - E\_YE\_X + E\_XE\_Y = E(XY) - E\_X E\_Y  ](https://latex.codecogs.com/png.latex?E%5B%28X-E_X%29%28Y-E_Y%29%5D%20%3D%20E%28XY%29%20-%20E_Y%20E_X%20-%20E_YE_X%20%2B%20E_XE_Y%20%3D%20E%28XY%29%20-%20E_X%20E_Y%20%20 "E[(X-E_X)(Y-E_Y)] = E(XY) - E_Y E_X - E_YE_X + E_XE_Y = E(XY) - E_X E_Y  ")

We can write the expression for the expectation of the random variable comprised of all pairs of values of ![X](https://latex.codecogs.com/png.latex?X "X") and ![Y](https://latex.codecogs.com/png.latex?Y "Y"), using the fact that for two independent random variables, ![P(x\_i,y\_j) = P(x\_i)P(y\_j)](https://latex.codecogs.com/png.latex?P%28x_i%2Cy_j%29%20%3D%20P%28x_i%29P%28y_j%29 "P(x_i,y_j) = P(x_i)P(y_j)") for all values ![x\_i](https://latex.codecogs.com/png.latex?x_i "x_i") and ![y\_j](https://latex.codecogs.com/png.latex?y_j "y_j"):

![E(XY) = \\sum\_i \\sum\_j x\_iy\_j P(x\_i,y\_j) = \\sum\_i x\_i P(x\_i) \\sum\_j y\_j P(y\_j) = E\_X E\_Y](https://latex.codecogs.com/png.latex?E%28XY%29%20%3D%20%5Csum_i%20%5Csum_j%20x_iy_j%20P%28x_i%2Cy_j%29%20%3D%20%5Csum_i%20x_i%20P%28x_i%29%20%5Csum_j%20y_j%20P%28y_j%29%20%3D%20E_X%20E_Y "E(XY) = \sum_i \sum_j x_iy_j P(x_i,y_j) = \sum_i x_i P(x_i) \sum_j y_j P(y_j) = E_X E_Y")

 The calculation for two continuous random variables is analogous, only with integrals instead of sums.

This demonstrates that the covariance of two independent random variables is 0, and thus that the variance of a sum of two independent random variables is the sum of the two separate variables.

**Example:** This property of variance is often used in analysis of noise or error in data. It is commonly assumed in least squares fitting that noise in data is independent of the signal or model underlying the data. This is the foundation for statements like "this linear regression explains 80% of the variance in the data."

Probability distribution examples (discrete)
============================================

The following are examples of distributions of random variables with discrete values. The first two have finite support (finitely many values) while the second two have infinite support.

Uniform
-------

The simplest probability distribution in which every value has the same probability (and one which is sometimes called "purely random" even though any random variable with any distribution is just as random). The probability distribution for a uniform random variable with ![n](https://latex.codecogs.com/png.latex?n "n") values is ![P(x) = 1/n](https://latex.codecogs.com/png.latex?P%28x%29%20%3D%201%2Fn "P(x) = 1/n") for any value ![x](https://latex.codecogs.com/png.latex?x "x").

``` r
low <- 0 # minimum value
high <- 10 # maximum value
values <- low:high # vector of discrete values of the RV
num <- length(values)
probs <- rep(1 / num, num) # uniform mass function vector
barplot(probs, names.arg = values, xlab = 'values', ylab = 'probability',
        main = paste("uniform  distribution on integers from ", low, "to ", high))
```

<img src="distributions_files/figure-markdown_github/unnamed-chunk-1-1.png" style="display: block; margin: auto;" />

``` r
unif.exp <- sum(values*probs)
paste("The expected value of uniform distribution is", unif.exp)
unif.var <- sum((unif.exp - values)^2*probs)
paste("The variance of uniform distribution is", unif.var)
```

    # [1] "The expected value of uniform distribution is 5"
    # [1] "The variance of uniform distribution is 10"

*Exercise:* experiment with the low and high values to see how the expectation and variance depend on them. Can you postulate a relationship without looking it up?

Binomial
--------

Binary or Bernoulli trials have two discrete outcomes (mutant/wildtype, win/lose, etc.). The number of "successes" out of a sequence of ![n](https://latex.codecogs.com/png.latex?n "n") independent binary trials with probability of success ![p](https://latex.codecogs.com/png.latex?p "p") is described by the binomial distribution.

``` r
n <- 10 # the number of trials
p <- 0.3 # the probability of success in one trial
values <- 0:n # vector of discrete values of the binomial
probs <- dbinom(values, n, p)
barplot(probs, names.arg = values, xlab = 'values', ylab = 'probability',
        main = paste("binomial distribution with n=", n, "and p=", p))
```

<img src="distributions_files/figure-markdown_github/unnamed-chunk-3-1.png" style="display: block; margin: auto;" />

``` r
bin.exp <- sum(values*probs)
paste("The expected value of binomial distribution is", bin.exp)
bin.var <- sum((bin.exp - values)^2*probs)
paste("The variance of binomial distribution is", bin.var)
```

    # [1] "The expected value of binomial distribution is 3"
    # [1] "The variance of binomial distribution is 2.1"

*Exercise:* Try different values of ![n](https://latex.codecogs.com/png.latex?n "n") and ![p](https://latex.codecogs.com/png.latex?p "p") and postulate a relationship with the expectation and variance.

Geometric
---------

The reandom variable is the first "success" in a string of independent binary trials and the distribution describes the probability of any nonnegative value. It may be pretty intuitive that since all the trials have the same probability of success, the distribution with have a geometric (exponential) form---try to figure out the exact formula for the probability density without looking it up!

``` r
p <- 0.3 # the probability of success
low <- 0 # minimum value
high <- 20 # maximum value
values <- low:high # vector of discrete values of the RV
probs <- dgeom(values, p)
barplot(probs, names.arg = values, xlab = 'values', ylab = 'probability',
        main = paste("geometric distribution with p=", p))
```

<img src="distributions_files/figure-markdown_github/unnamed-chunk-5-1.png" style="display: block; margin: auto;" />

``` r
geom.exp <- sum(values*probs)
paste("The expected value of geometric distribution is", geom.exp)
geom.var <- sum((geom.exp - values)^2*probs)
paste("The variance of geometric distribution is", geom.var)
```

    # [1] "The expected value of geometric distribution is 2.32030059650472"
    # [1] "The variance of geometric distribution is 7.52697882945385"

*Exercise:* Calculate the expectations and variances for different values of ![p](https://latex.codecogs.com/png.latex?p "p") and report how they are related.

Poisson
-------

Suppose that there is a discrete process that occurs with some average rate ![\\lambda](https://latex.codecogs.com/png.latex?%5Clambda "\lambda"), which describes the expected number of occurrences of these events in a unit of time. The Poisson random variable is the number of such occurrences, and the distribution describes the probability of any nonnegative value.

``` r
low <- 0 # minimum value
high <- 20 # maximum value
lambda <- 2 # Poisson rate
values <- low:high # vector of discrete values of the RV
probs <- dpois(values, lambda)
barplot(probs, names.arg = values, xlab = 'values', ylab = 'probability',
        main = paste("Poisson distribution with lambda=", lambda))
```

<img src="distributions_files/figure-markdown_github/unnamed-chunk-7-1.png" style="display: block; margin: auto;" />

``` r
pois.exp <- sum(values*probs)
paste("The expected value of Poisson distribution is", pois.exp)
pois.var <- sum((pois.exp - values)^2*probs)
paste("The variance of Poisson distribution is", pois.var)
```

    # [1] "The expected value of Poisson distribution is 1.99999999999987"
    # [1] "The variance of Poisson distribution is 1.99999999999777"

*Exercise:* Calculate the expectations and variances for different values of ![\\lambda](https://latex.codecogs.com/png.latex?%5Clambda "\lambda") and report how they are related.

Probability distribution examples (continuous)
==============================================

In the following examples with continuous variables we cannot calculate the means and variances directly from the density function. One way to do it is to produce a sample using the random number generator and calculate the mean and variance of that sample.

Uniform
-------

The continuous equivalent of the discrete uniform distribution.

``` r
low <- 0 # minimum value
high <- 10 # maximum values
number <- 100
values <- seq(low, high, length.out = number) # vector of discrete values of the RV
probs <- dunif(values, min=low, max = high)
plot(values, probs, t='l', xlab = 'values', ylab = 'density',
        main = paste("Uniform distribution on interval from ", low, "to ", high))
```

<img src="distributions_files/figure-markdown_github/unnamed-chunk-9-1.png" style="display: block; margin: auto;" />

``` r
n <- 1000 # sample size
unif.sample <- runif(n, low, high) # generate sample
unif.exp <- mean(unif.sample)
paste("The expected value of uniform distribution is", unif.exp)
unif.var <- var(unif.sample)
paste("The variance of uniform distribution is", unif.var)
```

    # [1] "The expected value of uniform distribution is 5.1684341000393"
    # [1] "The variance of uniform distribution is 8.26305092782467"

*Exercise:* experiment with the width of the interval to see how it affects the expectation and variance.

Exponential
-----------

The random variable describes the length of time between independent discrete events occurring with a certain rate, like we saw in the Poisson distribution.

``` r
low <- 0 # minimum value
high <- 20 # maximum values
number <- 100
r <- 0.5
values <- seq(low,high,length.out = number) # vactor of discrete values of the RV
probs <- dexp(values, r)
plot(values, probs, t='l', xlab = 'values', ylab = 'density',
        main = paste("Exponential distribution with rate=", r))
```

<img src="distributions_files/figure-markdown_github/unnamed-chunk-11-1.png" style="display: block; margin: auto;" />

``` r
n <- 1000 # sample size
exp.sample <- rexp(n, r) # generate sample
exp.exp <- mean(exp.sample)
paste("The expected value of exponential distribution is", exp.exp)
exp.var <- var(exp.sample)
paste("The variance of exponential distribution is", exp.var)
```

    # [1] "The expected value of exponential distribution is 2.112338991434"
    # [1] "The variance of exponential distribution is 5.09652107598591"

*Exercise:* What is the relationship between the rate and the expectation and variance?

Normal distribution
-------------------

The normal distribution, sometimes written ![N(\\mu, \\sigma)](https://latex.codecogs.com/png.latex?N%28%5Cmu%2C%20%5Csigma%29 "N(\mu, \sigma)") comes up everywhere (e.g., in the limit of the Poisson distribution for large ![n](https://latex.codecogs.com/png.latex?n "n")). The two parameters are simply the mean and the standard deviation. The reason for its ubiquity is that it is that any sum of a large number of independent random variables converges to the normal, formalized by the Central Limit Theorem:

> For a set of ![n](https://latex.codecogs.com/png.latex?n "n") IID random variables ![\\{X\_i\\}](https://latex.codecogs.com/png.latex?%5C%7BX_i%5C%7D "\{X_i\}") with mean ![\\mu](https://latex.codecogs.com/png.latex?%5Cmu "\mu") and standard deviation ![\\sigma](https://latex.codecogs.com/png.latex?%5Csigma "\sigma"), the sample mean ![\\bar X\_n](https://latex.codecogs.com/png.latex?%5Cbar%20X_n "\bar X_n") has the property:
>
> ![\\lim\_{n \\to \\infty} = \\frac{\\bar X\_n - \\mu}{\\sigma} = N(0,1)  ](https://latex.codecogs.com/png.latex?%5Clim_%7Bn%20%5Cto%20%5Cinfty%7D%20%3D%20%5Cfrac%7B%5Cbar%20X_n%20-%20%5Cmu%7D%7B%5Csigma%7D%20%3D%20N%280%2C1%29%20%20 "\lim_{n \to \infty} = \frac{\bar X_n - \mu}{\sigma} = N(0,1)  ")
>
>  where ![N(0,1)](https://latex.codecogs.com/png.latex?N%280%2C1%29 "N(0,1)") stands for the normal distribution with mean 0 and standard deviation 1.

``` r
low <- 0 # minimum value
high <- 10 # maximum values
number <- 100
mu <- 5
sigma <- 0.5 
values <- seq(low,high,length.out = number) # vactor of discrete values of the RV
probs <- dnorm(values, mu, sigma)
plot(values, probs, t='l',xlab = 'values', ylab = 'density',
        main = paste("Normal distribution with mean=", mu, "and sigma=", sigma))
```

<img src="distributions_files/figure-markdown_github/unnamed-chunk-13-1.png" style="display: block; margin: auto;" />

``` r
n <- 1000 # sample size
norm.sample <- rnorm(n, mu, sigma) # generate sample
norm.exp <- mean(norm.sample)
paste("The expected value of normal distribution is", norm.exp)
norm.var <- var(norm.sample)
paste("The variance of normal distribution is", norm.var)
```

    # [1] "The expected value of normal distribution is 4.99244615349142"
    # [1] "The variance of normal distribution is 0.248669598990351"

Application of normal distribution: confidence intervals:
=========================================================

The most important use of the normal distribution has to do with estimation of means, because the normal distribution describes the *sampling distributions* of means of IID samples. The mean of that sampling distribution is the mean of the population distribution that is being sampled, and the standard deviation is called the *standard error* and is related to the standard deviation of the population ![\\sigma\_X](https://latex.codecogs.com/png.latex?%5Csigma_X "\sigma_X") as follows: ![\\sigma\_{SE} = \\sigma/n](https://latex.codecogs.com/png.latex?%5Csigma_%7BSE%7D%20%3D%20%5Csigma%2Fn "\sigma_{SE} = \sigma/n"), where ![n](https://latex.codecogs.com/png.latex?n "n") is the sample size.

``` r
numsamples <- 1000
size <- 100
# compute mean for different samples
samplemeans <- replicate(n = numsamples, mean(sample(0:10, size, replace = TRUE)))
break_points <- seq(min(samplemeans), max(samplemeans), 
                    (max(samplemeans) - min(samplemeans)) / 20)
hist(samplemeans, breaks = break_points, freq = FALSE, 
     cex.axis = 1.5, cex.lab = 1.5,
     main= '1000 means of samples of size 100')
sigma <- 10 / sqrt(12) / sqrt(size)
mu <- 5
range <- seq(min(samplemeans), max(samplemeans), sigma / 100)
lines(range, 
      dnorm(range, mu, sigma),
      t = 'l', lwd = 3, col = 2, lty = 1, cex.axis = 1.5, cex.lab = 1.5)
```

<img src="distributions_files/figure-markdown_github/unnamed-chunk-15-1.png" style="display: block; margin: auto;" />

*Exercise:* Try using different distributions from above and see if the sample means still converge to the normal distribution.

The following script calculates a confidence interval based on a sample.

``` r
# Computing confidence intervals
qnorm(0.5) # the value that divides the density function in two
qnorm(0.95) # the value such that 95% of density is to its left 
size <- 100 # sample size
alpha <- 0.95 # significance level
sample <- runif(size)
s <- sd(sample) / sqrt(size) # standard error
z <- qnorm((1 - alpha) / 2) # z-value
left <- mean(sample) + s * z
right <- mean(sample) - s * z
print(right)
print(left)
```

    # [1] 0
    # [1] 1.644854
    # [1] 0.5736382
    # [1] 0.4546334

*Exercise:* Modify that script to report whether the confidence interval captures the true mean. Use a loop structure (as in the script above) to generate 1000 sample means and report how many of them are within the theoretical confidence interval. Does this match the fraction you expect from the significance level? Try different significance levels and sample sizes and report what you discover.
