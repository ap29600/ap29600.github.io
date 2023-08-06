---
title: APL Problem Solving Competition 2023
fontsize: 14pt
---

This year I participated in the APL Problem Solving Competition
organized by Dyalog Ltd, I will document my thoughts about the problems
I solved here.

# Part 1: Bioinformatics

## Task 1

the first problem starts off pretty tame: write a dfn that converts a DNA
sequence to RNA, where both input and outputs are character arrays.

The solution doesn't need much explanation, it's just replacing the `'T'`s
with `'U'`s.

```apl
rna  ← {'U'@('T'∘=)⍵}
```

## Task 2

I had a little more fun with the next one, as I avoided writing down two
separate look up tables by arranging my data properly:

the statement calls for a dfn that returns the reverse complement of a DNA
sequence.

- First we reverse the sequence: `{⌽⍵}`
- Then we find look up each element in the string `'ACGT'`: `{'ACGT'⍳⌽⍵}`
- The reverse of the same string dictates the replacements: `{(⊂'ACGT'⍳⌽⍵)⌷⌽'ACGT'}`

With a little massaging we end up with the final solution.

```apl
revc ← {'ACGT'(⊂⍤⍳⌷⌽⍤⊣)⌽⍵}
```

## Task 3

This task calls for a dfn that returns the protein coded by a RNA sequence,
terminating at the first stop codon. Setting up a lookup table is a logical
first step:

```apl
t ← ↑('UUU' 'F')('CUU' 'L')('AUU' 'I')('GUU' 'V')
t⍪← ↑('UUC' 'F')('CUC' 'L')('AUC' 'I')('GUC' 'V')
 ...
t⍪← ↑('UGG' 'W')('CGG' 'R')('AGG' 'R')('GGG' 'G')
r_cod amm ← ↓⍉t
```

The helper function `r_trs` does the heavy lifting, splitting the sequence
in codons of 3 (`x←{↓(⌊3÷⍨≢⍵)3⍴⍵}r`), looking up each one (`y←r_cod⍳x`) and
assinging the corresponding amino-acid (`amm⌷⍨∘⊂y`).

the whole function can be written tacitly, and the final solution just
applies it to it's input and truncates the result at the first `'$'`.

```apl
r_trs ← amm⌷⍨∘⊂r_cod⍳{↓(⌊3÷⍨≢⍵)3⍴⍵}
prot ← {'$'((¯1+⍳⍨)↑⊢)r_trs⍵}
```

## Task 4

This was a parsing problem, of which I've done entirely too many doing
Advent of Code in APL last year: it all boils down to partitioning and
enclosing and partition-enclosing.

We need to write a dfn that accepts the path to a fasta file and returns
a list of header-sequence pairs corresponding to the DNA sequences stored
in the file.

For reference, a fasta file looks something like this:

```fasta
>Rosalind_99
AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGA
TTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG
>Rosalind_2748
ATCAGGCTACCGTGTTTGCGGACGGGGGCTTAATCTAGCTTCTCATCTCAGCGACGTCTC
CTTGTTGGCACAGCGGTGGCAGGAGGTCCCCGCCGAGGAGCACATCGACCTTTCGGTGTA
...
```

The lines starting with `>` are headers, and the following lines contain the DNA sequences.

The solution is fairly pedestrian:
- Read the file's lines with `⊃⎕NGET⍵1`
- Identify the header lines with `('>'=⊃¨)`, and use them to partition-enclose the input.
- For each chunk:
    - `,/1↓⍵` concatenates all but the header
    - `⊂' '(1↓(¯1+⍳⍨)↑⊢)⊃⍵` gets the header up to the first space character.

```apl
readFASTA ← {{(⊂' '(1↓(¯1+⍳⍨)↑⊢)⊃⍵),,/1↓⍵}¨(('>'=⊃¨)⊂⊢)⊃⎕NGET⍵1}
```

## Task 5

The task is to write a function that takes a fasta file, extracts the first
(only) DNA sequence, and returns all the proteins coded by it, beginning with a
'M' amino-acid and ending with a stop codon. I approached the problem by
writing a helper function that would turn a sequence of DNA into the amino-acid
sequences of it's 6 reading frames.

Our final solution will look like this (for some definition of `crf` and `aas`):

```apl
orf ← {crf⍤aas⊃⌽⊃readFASTA⍵}
```

- We just defined `readFASTA`, which will give us the sequences.
- `⊃⌽⊃` takes the last element in the first header-sequence pair.
- `aas` will be defined later, it outputs the amino-acid sequences for the six reading
  frames of the DNA sequence.
- `crf` will be defined later, it accepts a list of amino-acid sequences and outputs
  their subsequences that begin with `M` and end right before a `$`.

### Defining `aas`

Here I had the choice between the easy and intuitive solution by composition,
or the "smart" solution that probably has no real benefit in this case.

I obviously chose the latter.

What I could have done was to extract the reading frames with the help of
`revc`, then convert each DNA sequence into RNA with `rna`, finally transcribe
into amino-acids with `r_trs`: It would have looked something like this:

```apl
aas ← r_trs¨rna¨(,0 1 2∘.↓,⍥⊂∘revc⍨)
```

however, RNA and DNA are dual representations of the same data, so I opted to
convert my lookup table into DNA representation, and write a function to
parallel `r_trs`: this saves us the RNA conversion, and for big inputs it may
perform better (though probably not by a significant amount).

```apl
d_cod ← 'T'@('U'∘=)¨r_cod
d_trs ← amm⌷⍨∘⊂d_cod⍳{↓(⌊3÷⍨≢⍵)3⍴⍵}
aas   ← d_trs¨(,0 1 2∘.↓,⍥⊂∘revc⍨)
```

### Retrospective

In hindsight, I could have been a little smarter still: observing that the
possible codons are just 3-digit numbers in the base 4 number system with
digits `'ACGT'`, we can leverage the decode primitive.

`{1+4⊥⍉¯1+'ACGT'⍳(⌊3÷⍨≢⍵)3⍴⍵}` can substitute the lookup into `d_cod`, and
`d_trs` can be rewritten as `d_trs ← amm⌷⍨∘⊂{1+4⊥⍉¯1+'ACGT'⍳(⌊3÷⍨≢⍵)3⍴⍵}`,
provided that `amm` contains the amino-acid sorted by the respective codon's
numeric value.

There is no more need for the `d_cod` or `r_cod` tables at all! we could even
parametrise the solution and pass either `'ACGT'` or `'ACGU'` to get
both `r_trs` and `d_trs`. The final definition would have looked something like
this:

```apl
trs ← {amm[1+4⊥⍉¯1+⍺⍳(⌊3÷⍨≢⍵)3⍴⍵]}
aas ← 'ACGT'∘trs¨(,0 1 2∘.↓,⍥⊂∘revc⍨)
```

But we would still be doing duplicate work in the lookups, as well as giving up some
opportunities for semplification: this is what it would look like after addressing
those issues.

```apl
amm ← 'KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV$Y$YSSSS$CWCLFLF'
aas ← {{amm[1+4⊥⍉(⌊3÷⍨≢⍵)3⍴⍵]}¨,0 1 2∘.↓,⍥⊂∘(5-⌽)¯1+'ACGT'⍳⍵}
```

A nice insight comes from this layout: the last base is ignored in most codons,
as can be observed in `4 4 4⍴amm` (many rows are constant).

### Defining `crf`

Speculation aside, to complete the task we would like - given a list of the
amino-acids coded by the reading frames - to extract any subsequence of one of
them, which starts with `'M'` and ends just before a `'$'`.

Now I find myself in a bit of an embarassing situation: I have two solutions to
this task, which diverge at this point. I definitely submitted the first one,
but I am not sure if I ended up overriding my submission with the second
solution; as a result, I am not sure which version ended up being graded. I
will just describe them both and worry about it some other time.

### Solution 1

We define a 1-modifier `_f` which takes a monadic function and applies it to
each value in `⍵`, then (expecting `⍺⍺` to return a list for each application)
concatenates the results together. This is heavily inspired by rust's
`flat_map`.

```apl
_f ← {⊃,/⍺⍺¨⍵}
```

To construct the solution:

- Partition each reading frame by `$`:
  - `'$'∘≠⊆⊢` is a tacit function that splits on `$`
  - `-'$'≠⊢/` returns `0` if the sequence ends in `$`, `¯1` otherwise.
    In the latter case, we need to drop the last chunk from the partitioning.
  - `(-'$'≠⊢/)↓'$'∘≠⊆⊢` puts them together to give the chunks that end just before `$`.
- Given each chunk:
  - `'M'∘(=⊂⊢)` partition-encloses at the positions of `M`.
  - `,⍨\⍤⌽` performs a reverse-catenate-scan in reverse (what a mouthful) which
    yields the subsequences starting with `M` and ending at the end of the chunk.

What's left is to do some plumbing with `_f` and make sure the results are unique:

```apl
crf ← {∪'M'∘(,⍨\⍤⌽(=⊂⊢))_f((-'$'≠⊢/)↓'$'∘≠⊆⊢)_f⍵}
```

### Solution 2

The second solution is arguably a little more interesting, as it doesn't
directly translate from something you might write in a scalar language.

- Concatenate all the reading frames together, separating with `'$'`: \
  `s d←(,/,(⊂(+\≢¨))),∘'$'¨⍵`
  - to avoid losing information, we also store the positions of the
    artificially inserted `$`s, which are available as the cumulative
    sum of the lengths.
- Find the positions of `M` and `$` in the merged input with `'M$'(⍸=)¨⊂s`.
- Match each of the `M`s to the corresponding `$`'s position with
  binary search: `i←{↑⍺(⍵[1+⍵⍸⍺])}/ ...`
  - note that there is always a `$` at the end of input, so each `M` finds
    a corresponding `$`.
- Remove any values where the stop codon was inserted in step 1: \
  `b e←¯1+↓(~i[2;]∊d)/i`
- Finally use the indices to select the relevant subsequences: `∪b↓¨e↑¨⊂s`.

Putting it all together, we get the following:

```apl
crf ← {
	s d←(,/,(⊂(+\≢¨))),∘'$'¨⍵
	b e←¯1+↓(~i[2;]∊d)/i←⊃{↑⍺(⍵[1+⍵⍸⍺])}/'M$'(⍸=)¨⊂s
	∪b↓¨e↑¨⊂s
}
```

### Retrospective

Looking back on it, the last line might have been better expressed as
`∪b{s[⍺+⍳⍵-⍺]}¨e` to avoid materializing too many big arrays, but the performance
implications of it depend on the implementation.

# Part 2: Potpourri

## Task 1

The task is to write a dfn that computes or validates a vehicle identification number;
there are quite a few details involved, but the long and short of it is that characters in
the VIN contribute to a total score according to their value and position, and the score
modulo 11 determines the correct value of the 9th digit.

The solution handles a few cases:

- if the VIN contains 16 legal characters, calculate what character can be placed
  at position 9 to make a valid VIN, and do so.
- if the VIN contains 17 legal characters, check that the 9th character matches the
  one that we calculate from the rest of the VIN.
- otherwise, return `¯1`.

```apl
vinc ← ⎕D,⎕A~'IOQ'
vin ← {
	~×/⍵∊vinc: ¯1
	16=≢⍵: (calc(⊣@9)⊢)8(↑,'0',↓)⍵
	17=≢⍵: (calc=9∘⌷)⍵
	¯1
}
```

We delegate the messy part to the `calc` helper function, that takes a 17 character vector
and returns the character that, placed at position 9, would make the VIN valid.

### Defining `calc`

I set up the tables:

- `num` contains the numeric values to substitute to the characters
- `mul` contains the values that digits at each position must be multiplied by.

```apl
num  ← ((65∘≤×1+9|83∘≤+-∘65)+(57∘≥×-∘48))⎕ucs vinc
mul  ← 10 0@8 9⊢17↑2↓⌽,⍨⍳10
calc ← {(⎕D,'X')⌷⍨1+11|mul+.×num⌷⍨⊂vinc⍳⍵}
```

The procedure is as follows:

- Look up each character's value with `num⌷⍨⊂vinc⍳⍵`
- Since the score is calculated via a sum-of-products,
  the choice falls on inner product: `mul +.× ...`
- take the modulo and pick a digit: `(⎕D,'X')⌷⍨1+11| ...`

### Retrospective

the tables defined above look very bad and are not at all insightful, what I should really have written is:

```apl
num  ← (⍳9),(⍳8),((⍳9)~6 8),1↓⍳9
mul  ← (⊢,10 0 9,⊢)⌽1↓⍳8
```

... or even just the raw numbers, it's not like we're starved for hard disk space.

## Task 2

This task asked for a function that sorts strings representing software
versions according to their release number. The structure of the input is
either:

- a string in the format `"[a-zA-Z]+-[a-zA-Z]+-[0-9]+\.[0-9]+\.[0-9]+"`
- a list of strings in the above format.

The solution leverages total array ordering, and the bulk of the problem
just boils down to parsing the numbers in the string.

```apl
parts ← (-∘~⍨/'.'⎕vfi⊃)@3('-'∘≠⊆⊢)
sortVersions ← {1=≡⍵: ,⊂⍵ ⋄ ((⊂∘⍋parts¨)⌷⊢)⍵}
```

If the input is a flat character array, then just return it wrapped in a
1-element list. If it's a list of strings, split each entry by `-` into the
fields (`org`, `package`, `version`). Then parse the three numbers in `version`
using `⎕vfi`: this suggests a nice extension of the problem statement: if one
of the fields in `version` is not a number, `⎕vfi` will return a `0`
verification bit and a `0` result: this means that the expression
`result-~verification` returns `result` if valid, `¯1` otherwise.

This allows us to sort versions like `foo-bar-10.x.0` in between
`foo-bar-9.9.9` and `foo-bar-10.0.0`, which is kind of nice to have.

With that out of the way, we have a list of entries in the form: \
`org package (major minor patch)`. \
Applying `⍋` gives us a sorting permutation, which by the rules of total
array ordering conforms to the problem statement.

## Task 3

The task is to write a dfn that, given a sorted list of coin denominations
(positive integers) as left argument and a total amount (a positive integer) as
right argument, returns a matrix, where each row contains the coefficients of a
linear combination of the elements in the left argument which sums to the right
argument. Of course, the matrix should be exhaustive and not contain any
duplicates.

Another way to look at this is that it returns the tallest matrix `res`
that satisfies the following constraints:

```apl
⍵∧.=res+.×⍺
res≡∪res
```

The solution I submitted is a classic recursive one, with just a little bit
of care put into the performance of it:

```apl
makeChange ← {
	0≠⍵|⍨∨/⍺: (0,≢⍺)⍴0
	1=≢⍺: ⍪⍵÷⊃⍺
	⊃⍪/(((⊂¯1↓⍺)∇¨⍵-(⊃⌽⍺)×⊢),¨⊢)0,⍳⌊⍵÷⊃⌽⍺
}
```

- It is a well known mathematical fact that for a finite set $A$ of integers,
  there exists a linear combination in $A$ with integer coefficients that sums to
  $n$ if and only if $\text{GCD}(A) | n$ (where the GCD of a finite set is given
  by the pairwise GCD reduction of the elements in some arbitrary ordering, and
  $|$ is the divisibility symbol). This allows us to immediately return
  `(0,≢⍺)⍴0` if `⍵` is not divisible by `∨/⍺`, which prunes some branches.

- The second line handles the general base case, and it is only reached if there
  is exactly one trivial solution. We can avoid checking divisibility because the
  first line catches the case where `0≠⍵|⊃⍺`.

- The third line does all the recursing:
    - `0,⍳⌊⍵÷⊃⌽⍺` determines how many coins of the largest denomination can be
      taken away from `⍵`.
    - `⍵-(⊃⌽⍺)×⊢` calculates the amount left over from each operation.
    - `(⊂¯1↓⍺)∇¨⍵-(⊃⌽⍺)×⊢` recurses on each case, removing the largest denomination.
    - a solution of the original problem can be formed by appending $n$ to a
      solution of the reduced problem that takes away $n$ coins of the largest
      denomination: therefore `(((⊂¯1↓⍺)∇¨⍵-(⊃⌽⍺)×⊢),¨⊢) ... ` gives a list of
      matrices of the solutions
    - `⊃⍪/` concatenates the solutions together.
    - It may look suspicious that we're not deduplicating the lines, but it can be
      proven that every solution gotten this way is already unique:
        - the first two cases always return zero or one result, and are trivially
          unique.
        - by inductive hypothesis, results returned by a recursive call are
          internally unique.
        - results obtained by two different recursive calls are distinct, since
          they are joined with distinct values of $n$.

### Retrospective

After some thinking, I came up with a more interesting solution for this
problem, which I didn't submit as the competition had already closed: the idea
is to reduce the amount of recursive calls as much as possible, by allowing the
function to take an array right argument and process values of `⍵` in parallel,
sharing the recursive subproblems between them.

```apl
makeChange ← {
	1≡≢⍺: (0=⍵|⍨⊃⍺)↑¨1 1∘⍴¨⍵
	d←¯1+⍳¨1+⌊⍵÷i←⊢/⍺
	u←∪∊r←⍵-i×d
	s←(¯1↓⍺)∇u
	⊃⍣(0≡≢⍴⍵)⊢r{⊃⍪/s[u⍳⍺],¨⍵}¨d
}
```
