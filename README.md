# Sequence Alignment: From Naive to Dynamic Programming

> Companion repository for the blog series **[Introduction to Sequence Alignment](https://karso.xyz/blog/intro-to-seq-alignment-1/)** — a learning and demonstration project that implements fundamental sequence alignment algorithms in R from scratch.

---

## What Is This About?

Sequence alignment is the process of comparing two or more biological sequences (like DNA) to find similarities. It is one of the most foundational tasks in bioinformatics, with applications ranging from inferring evolutionary relationships to predicting gene function.

This repository exists to make the *how* and *why* behind the algorithms tangible. Instead of relying on black-box tools, the scripts here implement the algorithms step by step in plain R, so you can follow the logic directly.

Two approaches to pairwise sequence alignment are covered:

1. **Naive Approach (Brute-Force)** — enumerate every possible alignment, score each one, and return the best. Simple to understand, but computationally catastrophic for real-world sequences.
2. **Needleman-Wunsch Algorithm** — a dynamic programming solution that finds the globally optimal alignment in polynomial time, making it practical even for longer sequences.

---

## Repository Structure

```
sequence-alignment-basic/
├── scripts/
│   ├── naive_alignment.R       # Brute-force alignment via recursive enumeration
│   └── needleman-wunsch.R      # Global alignment via dynamic programming
├── figure/                     # Output figures used in the blogpost
└── sequence-alignment-basic.qmd  # Quarto source document for the blogpost
```

---

## How to Use the Scripts

Both scripts are self-contained and can be sourced directly in any R session. No additional packages are required.

### Naive Alignment

```r
source("scripts/naive_alignment.R")

seq1 <- "ATCAG"
seq2 <- "GATCA"

naive_alignment(seq1, seq2,
                match_score    = 1,
                mismatch_score = -1,
                gap_score      = -1)
```

**Note:** Keep sequences to 6 characters or fewer. The number of possible alignments grows according to the Delannoy number — two sequences of length 100 already yield roughly 2 × 10⁷⁵ possibilities, which is computationally hopeless.

### Needleman-Wunsch (Global Alignment)

```r
source("scripts/needleman-wunsch.R")

seq1 <- "TCTTCACCACCATGGAGAAGGCGATACTGGATACATACATAGCATACATA"
seq2 <- "TATACGGCCATGGCATAGATTCGATCATGTACACAATGACATAGACAGTG"

result <- needleman_wunsch(seq1, seq2)

cat(" Seq1  :", result$align1, "\n",
    "Seq2  :", result$align2, "\n",
    "Score :", result$score)
```

This handles sequences of 500+ characters without issue.

---

## Read the Full Blogpost

The blogpost walks through the motivation, the math, and the R code in an accessible way — no prior bioinformatics background required.

👉 **[Introduction to Sequence Alignment, Part 1](https://karso.xyz/blog/intro-to-seq-alignment-1/)**

Topics covered in Part 1:
- Why sequence alignment matters in biology
- Why the naive approach is infeasible (with the Delannoy number as proof)
- How dynamic programming solves the problem efficiently
- Step-by-step walkthrough of the Needleman-Wunsch algorithm

---

## Reproducing the Blogpost

The `.qmd` file is the Quarto source for the blogpost. To render it locally, make sure you have [Quarto](https://quarto.org/) installed, then run:

```bash
quarto render sequence-alignment-basic.qmd
```

---

## What's Next

This is Part 1 of an ongoing series. Planned topics for future parts include:

- Local Alignment: Smith-Waterman Algorithm
- Introduction to Bioconductor packages in R
- Extracting biological insight from alignment results

---

*Written by [Karso](https://karso.xyz) · Contact: [X](https://x.com/karso_sp) | [email](mailto:karsosuryo.p@gmail.com)*
