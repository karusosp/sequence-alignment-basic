---
title: Introduction to Sequence Alignment
author: Karso Suryo Putro
date: today
format:
  hugo:
    code-fold: true
  pdf: default
---


One of the most important thing to do in biology is to compare biological information, which mostly exist in the form of sequences. This is called sequence alignment and is of great importance as the comparison between sequences can provide important biological interpretation, such as finding evolutionary distance, predicting sequence's structure and function based on annotated sequence database, or detecting anomalies in a biological system. It is not an exageration to say that sequence alignment is the backbone of bioinformatics, the study of complex biological data. Here, I would like to explore the common techniques of sequence alignment by implementing the algorithms in R from scratch with the hope of understanding these fundamental techniques even more.

In this blog, firstly I would like to describe the need of efficient algorithm by demonstrating that the naive approach to align two sequences, that is by iterating over all possible alignment and find the best one by some scoring system, is not at all feasible in the real world as the operating cost of such method is exceedingly high. And after showing that naive approach is not feasible, I will discuss about a more efficient approach in the form of **dynamic programming**, of which the difference between *local alignment* and *global alignment* is introduced.

# Introduction: The Basic of Sequence Alignment

Suppose we have two sequences, say both have the length of 6; the first sequence is "ACGTAG" and the second sequnce is "ACATAC". What is the most obvious way to compare those two sequences? One may think that the way to do it is by stacking both sequences and try to give some form of evaluation to it, that is by quantifying how similar are those sequences.

    ACGTAG
    || ||
    ACATAC

In this case we have four similar matched words (or nucleotides) and two mismatchs. One of the most obvious way to evaluate the alignment is by providing scoring system. For example, for each match a score of $+1$ will be added whereas mismatch will give a score of $-1$. In biological sequence, the mismatch usually could be interpreted as mutation event in the form of *substitution* where one nucleotide is substituted with another kind of nucleotide. Therefore the above alignment is having a score $4-2=2$. This is straightforward enough. But there's still problem. Suppose if we use the same scoring system, what if the sequences are best aligned in different way? For instance we have two sequences: "ACGTAG" and "TACGTC". If we do the same as the above, we get the following alignment where nothing is matched at all.

    ACGTAG
          
    TACGTC

In this case, the alignment give a score of $-12$. But if we take a look at both sequences, a different alignment where we shift one of the sequence to the left or right will give us better alignment score:

    -ACGTAG
     ||||
    TACGTC-

Here we have four matchs and four mismatchs, that give us a score of $4 + -(4)$ = 0, which is higher than the previous one. However, in doing so we introduce something new: a *gap* which is denoted by the dashed symbol. And this actually made a biological sense because sequences organism are subject to *insertion* and *deletion* and gap is the representation of it. Therefore, our previous scoring system can be updated to provide gap penalty which in our case is given arbitrarily, where each gap would give score of $-1$. Once we introduce gap to the alignment, there will be complexity because we can add gap to every possible position on the both sequences.

For example if we try to align two sequences, both with the length of just 2: "AT" and "GA"; we have these possible alignment:

    alignment: AT | -AT | -AT | A-T | A-T | AT- | AT- | 
               GA | GA- | G-A | GA- | -GA | G-A | -GA |
    score    : -2 | -1  | -3  | -3  | -3  | -3  | -3  |

    alignment: --AT | -A-T | -AT- | A--T | A-T- | AT-- |
               GA-- | G-A- | G--A | -GA- | -G-A | --GA |
    score    :  -4  |  -4  |  -4  |  -4  |  -4  |  -4  |

For each alignment, we use the previous scoring system and find the best alignment. In this case the best alignment is the one with $-1$ score. Thus, using a scoring system that capture match, mistmatch, and gap information, we could implement a form of algorithm in a computer machine. This is called *naive approach* to the alignment problem because we try to naively brute-force every possible alignment and find the best one

## Implementing Naive Approach Algorithm in R

Here, I try to implement the brute-force approach of sequence alignment in R. The basic algorithm for the naive approach algorithm is given as follow:

1.  CONSTRUCT matrix representation for the two sequences that will be aligned
2.  GENERATE all the possible alignment recursively
3.  SCORE each of the alignment
4.  RETURN the alignment with the best score

Those are easy enough to understand, but try to be careful with the recursive programming part as it is the one which somewhat harder to understand, at least for me. I have constructed the R script that follow that algorithm's logic faithfully which you can see at [my github repo](https://github.com/karusosp/sequence-alignment-basic).

For example, suppose we have two sequences, each has the length of just 5: "ATCAG" and "GATCA". If I run my `naive_algorithm()` function, I will get the following:

<details class="code-fold">
<summary>Code</summary>

``` r
source("scripts/naive_alignment.R")
seq1 = "ATCAG"
seq2 = "GATCA"
naive_alignment(seq1, seq2, 
                match_score = 1,
                mismatch_score = -1,
                gap_score = -1)
```

</details>

    Total number of possible alignments: 1683 

    BEST ALIGNMENT(S) with score: 2 
    -ATCAG 
    GATCA- 

As you can see from the result, the function iterate over all possible alignments and the total possibility is $1683$ when the sequence is just 5 characters length. What if the length were higher then? Well, I would not recommend you try to run my script beyond 6 character length. The operation cost is too high and you can imagine how impractical if it is used in the real world scenario where most biological sequences have the length of more than one hundred. In fact, I would demonstrate the amount of all possible alignment for two 100-character long sequences. The formula for doing so is given by Delannoy Number:

$$
D(m, n) = \sum_{k=0}^{\min(m,n)} {m \choose k} {n \choose k} 2^k
$$

Where $m$ and $n$ is the amount of first and second sequence, whereas $k$ denote the amount of non-gap in each alignment (for further understanding, read [wikipedia page for Delannoy Number](https://en.wikipedia.org/wiki/Delannoy_number)). And if we implement that formula in an R code, we get the following:

<details class="code-fold">
<summary>Code</summary>

``` r
delannoy <- function(m, n) {
  k <- 0:min(m, n)
  sum(choose(m, k) * choose(n, k) * 2^k)
}

delannoy(100, 100)  
```

</details>

    [1] 2.053717e+75

You see that just with two 100-character long sequences, we almost reach **the eddington number** (approximately $1.57 \times 10^{79}$), a number that represent **the total amount of all atom in the universe!**. We could also illustrate it by creating a graph where the number of possibility increase on logarithmic scale, where as you can see in the below graph, 20-long sequence already contain more than $10^{10}$ possible alignments, which is already computationally hopeless.

<img src="sequence-alignment-basic.markdown_strict_files/figure-markdown_strict/all%20possible%20allignment%20plot-1.png" style="width:60.0%" data-fig-align="center" />

# Smarter Approach: Introduction to Dynamic Programming

As we've discussed earlier, the naive approach to sequence alignment cannot work in real life precisely because it need to run for all possible alignment. A different method to find the best alignment therefore is needed to ensure optimal running time and efficient computational cost. The way to do this is by a programming technique called **dynamic programming**, which refer to a class of computing technique where a large and complex problem is broken down into simpler and overlapping sub-problems. But how exactly does this technique is being applied in the context of sequence alignment?

In order to understand, let's firslty clarify our main problem in alignment before. As we've discussed, the sequence alignment can be represented as a two-dimensional matrix. The main problem is how to find the best possible path from the START position, denoted as (1,1) in R, to the END position, denoted as (nrow, ncol). And here the best path is determined by a scoring system that is being applied to the all possible resulting alignment.

The key insight is instead of applying the scoring for each alignment, **the score can be immediately determined for each cell within the matrix based on a recurring simple question**: *what is the maximum score for this cell given three possible moves: diagonal (match = +1/mismatch = -1), downward and rightward (gap = -1)?*. And that question is being repetitively asked for every cell within the matrix. And this is much simpler problem than determining all the possible path from the start to end position. If that is a bit too abstract, let's try to imagine it.

::: {#fig-complex-label} <img src="figure/scoring_mat_step.png" style="width:65.0%" /> Step-by-step illustration for filling the alignment matrix with scores. The image is taken from the awesome [Kenko Wong's blog](https://www.kenkoonwong.com/blog/dynamic-programming/) :::

And as you can see in the image above, the recurring problem are highlighted with red boxes. And after all the cells within the matrix are filled with the best score, we can find the best alignment by tracing back the path it require to get to the bottom-right cell and this is a trivial matter because we can directly store the move that correspond with the best score during the matrix filling process.

Why this works? The justification for this approach is the fact that for an optimal sequence alignment, the subsequence alignment also needed to be optimal as well. Therefore we can built from the ground up, from the smallest possible subsequence alignment and reuse the solution to gradually built up the full sequence alignment. And this made sense because if the subsequence alignment is unoptimal, then the whole sequence alignment is not the most optimal.

By using this dynamic programming technique, we reduce the scale problem from having to calculate the whole possibility which increase exponentially (approximately $2^{m \times n}$), to a process where the scale only increase polynomially (approximately $m \times n$).

## Needleman-Wuncsh Algorithm (Global Alignment)

One of the gold standard in sequence alignment that uses dynamic programming is the Needleman-Wuncsh algorithm. The core logic of the algorithm can be described in few steps:

1.  GENERATE the 2D matrix representation of the sequence alignment
2.  FILL each interior cell with the maximum score from three possible moves: diagonal (match/mismatch), downward (gap), rightward (gap) and simultaneously record which move produced the best score in a pointer matrix.
3.  TRACEBACK \[nrow, ncol\] to \[1, 1\] by following the pointer matrix, and directly build the two aligned sequences.

For demosntration purpose, I have implemented the algorithm in R which you can directly access within this github repo: [sequence-alignment-basic](https://github.com/karusosp/sequence-alignment-basic). And below is the result of aligning two 50-character long sequences using the script.

<details class="code-fold">
<summary>Code</summary>

``` r
source("scripts/needleman-wunsch.R")
seq1 <- "TCTTCACCACCATGGAGAAGGCGATACTGGATACATACATAGCATACATA"
seq2 <- "TATACGGCCATGGCATAGATTCGATCATGTACACAATGACATAGACAGTG"
result <- needleman_wunsch(seq1, seq2)
cat(" Seq1  :", result$align1, "(Seq Length ", nchar(seq1), ")", "\n",
    "Seq2  :", result$align2, "(Seq Length ", nchar(seq2), ")", "\n",
    "Score :", result$score )
```

</details>

     Seq1  : TCTTCACCACCATGG-AGA-AGGCGATAC-TGGATAC-AT-ACATAGCATACA-TA (Seq Length  50 ) 
     Seq2  : T-AT-ACGGCCATGGCATAGATTCGAT-CATGTACACAATGACATAG---ACAGTG (Seq Length  50 ) 
     Score : 14

As you can see, the algorithm works wonder, even despite my bad code implementation. It can run for over 500bp sequence without problems. In fact, with much better code implementation it can handle up to thousands sequence. The fact that it is working very fast and with reliable result made this technique a gold standard in sequence alignment. However, it must be noted that this algorithm can only be applied as pairwise-alignment (two sequences only) and global-alignment (it align the whole sequence). Different cases of alignment, for instance multiple sequence alignment or local alignment between unequal-sized sequences, need different algorithms.

# Smith-Waterman Algorithm

Suppose we have two sequences with different length. What we would love to do is to find the region for the larger sequence that correspond with our short sequence with the highest similarity. This kind of case is similar to sequence database search where we want to find out which region from the database sequence where we want to find the kind of sequence we have in hand. As we've discussed previously, using global alignment for this kind of case is not appropriate. What we should do is to align our sequences using one form of local alignment.

<details class="code-fold">
<summary>Code</summary>

``` r
seq1 = "GACATAGACAGATACACAGATAGACAGATAGACAGAGATGACACACCGTCGTCACAGTCAGATCAGATGGGATAGCCCAGAGTTTGCACAGTAGC"
seq2 = "CACACAGTTTAC"
```

</details>

We want to align those two sequences using SWA. The step for doing so is actually somewhat similar to the NWA. What we want to do firstly is to create 2 dimensional matrix representation of our sequences.

<details class="code-fold">
<summary>Code</summary>

``` r
generate_2d_matrix = function(seq1, seq2){
  n_col <- nchar(seq1) + 1
  n_row <- nchar(seq2) + 1

  mat         <- matrix(NA, nrow = n_row, ncol = n_col)
  colnames(mat) <- c(".", strsplit(seq1, "")[[1]])
  rownames(mat) <- c(".", strsplit(seq2, "")[[1]])

  return(mat)
}
mat <- generate_2d_matrix(seq1,seq2)
```

</details>

After we create the matrix representation, what we want to do is to fill the matrix's cells with score. The scoring schema in this step is different from Needleman-Wunsch Algorithm. In that algorithm, we fill the first row and first column with the gap penalty. However, in SWA we just put 0 in both the first row and column.

<details class="code-fold">
<summary>Code</summary>

``` r
score_matrix = function(mat, 
                        seq1, seq2, 
                        match, 
                        mismatch, 
                        gap_penalty){
  n_row <- nrow(mat)
  n_col <- ncol(mat)
  
  seq1_chars <- strsplit(seq1, "")[[1]]
  seq2_chars <- strsplit(seq2, "")[[1]]
  
  # Boundary initialization — all zeros for local alignment
  mat[1, ] <- 0
  mat[, 1] <- 0
  
  # Fill interior cells
  for (i in 2:n_row){
    for (j in 2:n_col){
      
      # Diagonal: match or mismatch between seq2[i-1] and seq1[j-1]
      # (subtract 1 because row/col 1 is the boundary, not a character)
      if (seq2_chars[i - 1] == seq1_chars[j - 1]){
        score_diagonal <- mat[i-1, j-1] + match
      } else {
        score_diagonal <- mat[i-1, j-1] + mismatch
      }
      
      # Up: gap inserted in seq1
      score_above <- mat[i-1, j] + gap_penalty
      
      # Left: gap inserted in seq2
      score_left  <- mat[i, j-1] + gap_penalty
      
      # Zero-floor is part of the recurrence, not a post-processing step
      mat[i, j]   <- max(0, score_diagonal, score_above, score_left)
    }
  }
  
  return(mat)
}

mat <- score_matrix(mat, 
                    seq1, seq2, 
                    match = 1, 
                    mismatch = -1,
                    gap_penalty = -1)
```

</details>

After the matrix is filled with scores, now we worked backward to find the most optimal path that correspond to the optimal local alignment between two sequences. To do this, we firstly need to identify the cell with the highest score and walk from there following the best possible route until it encounter zero cell, in which the alignment stop.

<details class="code-fold">
<summary>Code</summary>

``` r
traceback = function(mat, seq1, seq2, match, mismatch, gap_penalty){
  
  seq1_chars <- strsplit(seq1, "")[[1]]
  seq2_chars <- strsplit(seq2, "")[[1]]
  
  # --- Find starting cell (highest score) ---
  best_score <- 0
  best_i     <- 1
  best_j     <- 1
  
  for (i in 1:nrow(mat)){
    for (j in 1:ncol(mat)){
      if (mat[i, j] > best_score){
        best_score <- mat[i, j]
        best_i     <- i
        best_j     <- j
      }
    }
  }
  
  # --- Walk backwards until hitting a zero cell ---
  aligned_seq1 <- ""
  aligned_seq2 <- ""
  
  i <- best_i
  j <- best_j
  
  while (mat[i, j] != 0){
    
    current_score <- mat[i, j]
    
    # Recompute each neighbor's contribution to current cell
    if (seq2_chars[i - 1] == seq1_chars[j - 1]){
      diag_score <- mat[i-1, j-1] + match
    } else {
      diag_score <- mat[i-1, j-1] + mismatch
    }
    up_score   <- mat[i-1, j] + gap_penalty
    left_score <- mat[i, j-1] + gap_penalty
    
    # Move to the neighbor that produced the current score
    if (diag_score == current_score){
      # Match or mismatch — both sequences contribute a character
      aligned_seq1 <- paste0(seq1_chars[j - 1], aligned_seq1)
      aligned_seq2 <- paste0(seq2_chars[i - 1], aligned_seq2)
      i <- i - 1
      j <- j - 1
      
    } else if (up_score == current_score){
      # Gap in seq1 — only seq2 contributes a character
      aligned_seq1 <- paste0("-", aligned_seq1)
      aligned_seq2 <- paste0(seq2_chars[i - 1], aligned_seq2)
      i <- i - 1
      
    } else {
      # Gap in seq2 — only seq1 contributes a character
      aligned_seq1 <- paste0(seq1_chars[j - 1], aligned_seq1)
      aligned_seq2 <- paste0("-", aligned_seq2)
      j <- j - 1
    }
  }
  
  return(list(
    score    = best_score,
    seq1     = aligned_seq1,
    seq2     = aligned_seq2
  ))
}

traceback(mat, seq1, seq2, 1, -1,-1)
```

</details>

    $score
    [1] 6

    $seq1
    [1] "ACACAG"

    $seq2
    [1] "ACACAG"

After we work our way to the end and get the result we want, we can summarize our code into:

# Real Case Example
