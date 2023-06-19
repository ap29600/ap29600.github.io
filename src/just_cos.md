---
geometry: margin=1in
title: Just `cos`
subtitle: lessons learned implementing trig functions without floating point numbers
---

## Notation:
- $\sin$ and $\cos$ are mathematical functions taking the argument in radians
- $\text{tsin}$ and $\text{tcos}$ are mathematical functions taking the
  argument in turns
- `sin` and `cos`   are the `C` functions that calculate $\sin$ and $\cos$ on
  `double`s
- `tsin` and `tcos` are the `C` functions that calculate $\text{tsin}$ and
  $\text{tcos}$ on `double`s 
- `isin` and `icos` are the `C` functions that calculate $\text{tsin}$ and
  $\text{tcos}$ on fixed point numbers

## Choosing turns over radians
This is a tough decision for a mathematician to do, but $\pi$ is irrational and
computer numbers are rationals. What's more, they are (almost always) multiples
of a power of 2. This means that, in floating point numbers of __any__
precision, `sin`$(\pi) \ne 0$, because the argument is rounded before being
passed to the function. If we use turns, though, more of the meaningful angles
are exactly representable: `tsin`$(\frac 1 2) = 0$; in fact, $\forall n \in
\mathbb{N},$ `tsin`$(\frac n 2) = 0$ and `tcos`$(\frac n 2) = \pm 1$, always,
with exactly zero error.

Moreover (and this is completely anecdotal) the most common operation that a
user does before feeding an exact value to `sin` or `cos` is likely to multiply
it by some integer multiple of $\pi$, and the first thing that (at least some)
mathematical library does is immediately divide by $2\pi$.

## Polynomial interpolation
A common strategy for evaluating transcendental functions on computer hardware
is to find a good enough interpolating polynomial for a small subset of the
function's domain, and then to refer the other values back to that subset by
some symmetry.

For our purposes, we will first focus on finding an interpolating polynomial
for the $\cos$ function on the range $[-\frac\pi 4, \frac\pi 4)$ (or rather for
the $\text{tcos}$ function on the corresponding range $[-\frac 1 8, \frac 1
8)$). Then the $\text{tcos}$ function can be calculated by multiplying the angle by 4:
then the integer part mod 4 tells us what function to apply:

- if $\left\lfloor4x\right\rceil \equiv 0 \pmod{4}$, calculate $+$`icos`$\left(2^{63} \cdot 8\left(x - \frac{\left\lfloor4x\right\rceil}{4}\right)\right)$
- if $\left\lfloor4x\right\rceil \equiv 1 \pmod{4}$, calculate $-$`isin`$\left(2^{63} \cdot 8\left(x - \frac{\left\lfloor4x\right\rceil}{4}\right)\right)$
- if $\left\lfloor4x\right\rceil \equiv 2 \pmod{4}$, calculate $-$`icos`$\left(2^{63} \cdot 8\left(x - \frac{\left\lfloor4x\right\rceil}{4}\right)\right)$
- if $\left\lfloor4x\right\rceil \equiv 3 \pmod{4}$, calculate $+$`isin`$\left(2^{63} \cdot 8\left(x - \frac{\left\lfloor4x\right\rceil}{4}\right)\right)$

We want a few good properties of this polynomial; first, we want it to be a
close approximation, which means that for all values of $x \in [-\frac 1 8,
\frac 1 8)$, the absolute error $|p(x) - \text{tcos}(x)|$ is smaller than some
defined threshold. Since the goal is implementing trigonometric functions for
posit arithmetic, we're looking for an error in the order of $2^{-61}$, which
corresponds to about 1ULP (unit in last place) for the posit value
$\frac{\sqrt{2}}{2}$ (which is the smallest output value we need to care about);
this also means that we need to do significantly better
than `double` precision arithmetic, which only has 53 bits of precision and
as such can get away with an error below $2^{-54}$.

Second, we want it to be cheap to evaluate: for polynomials this means to have
as small of a degree as possible, in our case a degree 14 polynomial is enough,
however the terms of odd degree will be zero due to $\text{tcos}$
being an even function, and only 6 terms (plus the constant term) will be used.

## I'm going to be real with you chief, I'm not doing that by hand.
Obviously, finding the best polynomial to approximate $\text{tcos}$ would take
forever to do by hand, considering that the coefficients need about the same
number of significant digits that we want in our result.
Fortunately interpolating polynomials are very well researched, and getting
a computer to do the heavy lifting for us will be pretty easy; the first order of
business is getting some accurate samples of the $\text{tcos}$ function.

Standing on the shoulders of giants, we know that for maximum numerical accuracy
it's a good idea __not__ to use equidistant samples, which would incur
[Runge's phenomenon](https://en.wikipedia.org/wiki/Runge%27s_phenomenon), but
instead to choose a set of
[Chebyshev nodes](https://en.wikipedia.org/wiki/Chebyshev_nodes), obtained through the formula: 
$$
x_k = a + (b-a) \cos\left(\frac k n \pi\right)
$$
where $[a, b]$ is the interval we want to approximate our function on.

Ok then, we ask Wolfram Alpha (which supposedly calculates $\cos$ accurately
enough for our purposes) with a series of queries like this:

`round(cos(k/14 * pi) * 2^63)`

`round(cos(cos(k/14 * pi) * pi/4) * 2^63)`

and get back our data:

 | $x \cdot 2^{63}$        | $y \cdot 2^{63}$      |
 | ---                     | ---                   |
 | `-9223372036854775808`  | `6521908912666391106` |
 | `-8992122843167040398`  | `6649062829911898522` |
 | `-8309971062287876938`  | `7008946939873054357` |
 | `-7211122632928341039`  | `7538462657309403407` |
 | `-5750678403727967667`  | `8139439240617060549` |
 | `-4001871146622878209`  | `8693001419304471778` |
 | `-2052393359867478633`  | `9082872288129848112` |
 | `+0000000000000000000`  | `9223372036854775808` |
 | `+2052393359867478633`  | `9082872288129848112` |
 | `+4001871146622878209`  | `8693001419304471778` |
 | `+5750678403727967667`  | `8139439240617060549` |
 | `+7211122632928341039`  | `7538462657309403407` |
 | `+8309971062287876938`  | `7008946939873054357` |
 | `+8992122843167040398`  | `6649062829911898522` |
 | `+9223372036854775808`  | `6521908912666391106` |

now it's time to do some actual math. I chose haskell as a language because it has
good support for arbitrary precision rational numbers: if any of our values so
much as touch a floating point representation, we will be losing precision.

```haskell

xs' = [-9223372036854775808
      ,-8992122843167040398
      ,-8309971062287876938
      ,-7211122632928341039
      ,-5750678403727967667
      ,-4001871146622878209
      ,-2052393359867478633
      , 0000000000000000000
      , 2052393359867478633
      , 4001871146622878209
      , 5750678403727967667
      , 7211122632928341039
      , 8309971062287876938
      , 8992122843167040398
      , 9223372036854775808
      ]

ys' = [6521908912666391106
      ,6649062829911898522
      ,7008946939873054357
      ,7538462657309403407
      ,8139439240617060549
      ,8693001419304471778
      ,9082872288129848112
      ,9223372036854775808
      ,9082872288129848112
      ,8693001419304471778
      ,8139439240617060549
      ,7538462657309403407
      ,7008946939873054357
      ,6649062829911898522
      ,6521908912666391106
      ]

xs::[Rational]
xs = map (/(2^63)) xs'

ys::[Rational]
ys = map (/(2^63)) ys'

```

One way to find the interpolating polynomial given a set of nodes is through
[Newton's polynomial](https://en.wikipedia.org/wiki/Newton_polynomial#Main_idea)
form: we construct the lower triangular matrix as a list of lists
(here I rounded the values for printing).

```haskell
newtonMatrix = map (takeWhile (/=0) . scanl (*) 1 . zipWith (flip (-)) xs . repeat) xs
```

```
[[1.0]
 [1.0,2.5072087818176395e-2]
 [1.0,9.903113209758087e-2,7.324247883844538e-3]
 [1.0,0.2181685175319702,4.212756181137467e-2,5.0189675689328046e-3]
 [1.0,0.37651019814126646,0.13232003255213892,3.671603905143758e-2,5.81367932 ...
 [1.0,0.5661162608824418,0.30629390422737474,0.1430653277020612,4.97792579256 ...
 [1.0,0.7774790660436856,0.5849806747155206,0.3968789301591433,0.221978572120 ...
 [1.0,1.0,0.9749279121818236,0.8783796973249267,0.6867449009293668,0.42817844 ...
...]
```

then the coefficients to Newton's form can be obtained by solving the triangular system:

```haskell
newtonCoefficients = zipWith fn newtonMatrix ys
  where fn row y = (y - ((init row) `dot` newtonCoefficients)) / (last row)
        as `dot` bs = sum $ zipWith (*) as bs
```

```
[0.7071067811865476,0.5498566944938549,-0.22502818791042026,-5.310916131098175 ...
```

finally we can convert the coefficients from Newton's form to power series form, by
distributing the products:

```haskell
powerSeries [c] _ = [c]
powerSeries (c:cs) (x:xs) = (c+constantTerm):rest
  where psc' = powerSeries cs xs
        (constantTerm:rest) = zipWith (-) (0:psc') $ map (x*) (psc'++[0])

coefficients = powerSeries newtonCoefficients xs
```

then we can print our result, converted to fixnum:
```haskell
main = mapM (print . round . (*(2^63))) coefficients
```

```
9223372036854775808
0
-2844719788994575527
0
146230515361076815
0
-3006744454122068
0
33119841825907
0
-226999764641
0
1060732338
0
-3557527
```

The first is the constant term, which we can verify is equal to $2^{63}$, and
represents `icos`$(0)$; after that every other term is zero, which matches up
with our prediction; the remaining terms alternate signs, similar to the Taylor
coefficients for $\cos$' power series. In fact, they are nearly identical
accounting for the rescaling.

## Evaluating the polynomial
The formula we want to evaluate now is:
\begin{align*}
    \text{tcos}(\alpha) \sim 1
    &- \frac{2844719788994575527}{2^{63}} \alpha^2
     + \frac{146230515361076815} {2^{63}} \alpha^4
     - \frac{3006744454122068}   {2^{63}} \alpha^6 \\
    &+ \frac{33119841825907}     {2^{63}} \alpha^8
     - \frac{226999764641}       {2^{63}} \alpha^{10}
     + \frac{1060732338}         {2^{63}} \alpha^{12}
     - \frac{3557527}            {2^{63}} \alpha^{14}
\end{align*}

since our input and output are fixed point numbers, what we will end up is:

\begin{align*}
    \text{icos}(t) = 2^{63}
    &- 2844719788994575527 \left(\frac{t}{2^{63}}\right)^2
     + 146230515361076815  \left(\frac{t}{2^{63}}\right)^4
     - 3006744454122068    \left(\frac{t}{2^{63}}\right)^6 \\
    &+ 33119841825907      \left(\frac{t}{2^{63}}\right)^8
     - 226999764641        \left(\frac{t}{2^{63}}\right)^{10}
     + 1060732338          \left(\frac{t}{2^{63}}\right)^{12}
     - 3557527             \left(\frac{t}{2^{63}}\right)^{14}
\end{align*}
Now, we obviously can't calculate $\frac t {2^{63}}$ right away, as that would
always be $0$ or $\pm 1$ when rounded to an integer and we would lose all our
precision; we also would like as many of our divisions to have a denominator of
exactly $2^{64}$, since that's equivalent to taking the higher word of a merged
128-bit register, which on x86_64 is very cheap (good ol' `_mulx_u64`).

we perform the following substitution:
\begin{align*}
    x_2 = \frac{t^2}{2^{62}} &= \left(\frac{t}{2^{63}}\right)^2 2^{64} \\
    x_{2n+2} = \frac{x_{2n} \cdot x_2}{2^{64}} &= \left(\frac{t}{2^{63}}\right)^{2n+2} 2^{64} \\
    \text{icos}(t) = 2^{63}
    &- \frac{x_2    \cdot 2844719788994575527}{2^{64}}
     + \frac{x_4    \cdot 146230515361076815 }{2^{64}}
     - \frac{x_6    \cdot 3006744454122068   }{2^{64}}\\
    &+ \frac{x_8    \cdot 33119841825907     }{2^{64}}
     - \frac{x_{10} \cdot 226999764641       }{2^{64}}
     + \frac{x_{12} \cdot 1060732338         }{2^{64}}
     - \frac{x_{14} \cdot 3557527            }{2^{64}}
\end{align*}
We implement that in `C` as:
```c
uint64_t icos(int64_t t) {
	static const uint64_t c[] = {
		2844719788994575527,
		146230515361076815,
		3006744454122068,
		33119841825907,
		226999764641,
		1060732338,
		3557527,
	};

	uint64_t x2  = ((int128_t)t    * (int128_t)t)   >> 62;
	uint64_t x4  = ((uint128_t)x2  * (uint128_t)x2) >> 64;
	uint64_t x6  = ((uint128_t)x4  * (uint128_t)x2) >> 64;
	uint64_t x8  = ((uint128_t)x6  * (uint128_t)x2) >> 64;
	uint64_t x10 = ((uint128_t)x8  * (uint128_t)x2) >> 64;
	uint64_t x12 = ((uint128_t)x10 * (uint128_t)x2) >> 64;
	uint64_t x14 = ((uint128_t)x12 * (uint128_t)x2) >> 64;

	uint64_t res = (uint64_t)INT64_MAX + 1;
	res -= (uint128_t)x2  * (uint128_t)c[0] >> 64;
	res += (uint128_t)x4  * (uint128_t)c[1] >> 64;
	res -= (uint128_t)x6  * (uint128_t)c[2] >> 64;
	res += (uint128_t)x8  * (uint128_t)c[3] >> 64;
	res -= (uint128_t)x10 * (uint128_t)c[4] >> 64;
	res += (uint128_t)x12 * (uint128_t)c[5] >> 64;
	res -= (uint128_t)x14 * (uint128_t)c[6] >> 64;

	return res;
}
```

Here we're exploiting the fact that modern compilers offer support for the
types `__int128` and `unsigned __int128`, which we don't really use fully: we
only care about the upper half of the result which on x86_64 gets stored in
`rdx` after a multiplication.

Now all our multiplications are between numbers in the range $[0,2^{64})$,
there's no integer division in sight and... wait, what's that?

```c 
icos(INT64_MIN);
9223372036854775808
```

Yeah uhm, it turns out that `((int128_t)INT64_MIN * (int128_t)INT64_MIN) >> 62`
is exactly `UINT64_MAX + 1`, and it makes us overflow. For reference,
`INT64_MIN` was one of our sample points and the correct value there was
`6521908912666391106`.

The fix is pretty easy, just subtract 1 from `x2` if `t==INT64_MIN` to take it
back in the representable range of `uint64_t`. We will verify that the error
introduced this way is negligible.
```c
uint64_t icos(int64_t t) {
	static const uint64_t c[] = {
		2844719788994575527,
		146230515361076815,
		3006744454122068,
		33119841825907,
		226999764641,
		1060732338,
		3557527,
	};

	uint64_t x2  = ((int128_t)t * (int128_t)t - (t==INT64_MIN)) >> 62;
	uint64_t x4  = ((uint128_t)x2  * (uint128_t)x2) >> 64;
	uint64_t x6  = ((uint128_t)x4  * (uint128_t)x2) >> 64;
	uint64_t x8  = ((uint128_t)x6  * (uint128_t)x2) >> 64;
	uint64_t x10 = ((uint128_t)x8  * (uint128_t)x2) >> 64;
	uint64_t x12 = ((uint128_t)x10 * (uint128_t)x2) >> 64;
	uint64_t x14 = ((uint128_t)x12 * (uint128_t)x2) >> 64;

	uint64_t res = (uint64_t)INT64_MAX + 1;
	res -= (uint128_t)x2  * (uint128_t)c[0] >> 64;
	res += (uint128_t)x4  * (uint128_t)c[1] >> 64;
	res -= (uint128_t)x6  * (uint128_t)c[2] >> 64;
	res += (uint128_t)x8  * (uint128_t)c[3] >> 64;
	res -= (uint128_t)x10 * (uint128_t)c[4] >> 64;
	res += (uint128_t)x12 * (uint128_t)c[5] >> 64;
	res -= (uint128_t)x14 * (uint128_t)c[6] >> 64;

	return res;
}
```
