#-----------------------------------------------------------------------------#
#                     Needleman-Wunsch Algorithm                              #
#-----------------------------------------------------------------------------#
# The needleman_wunsch() is a function that takes two sequences and
# finds the best global alignment using the Needleman-Wunsch algorithm.
#
# The general algorithm is described in the following steps:
# 1. Generate a 2D (n+1) x (m+1) matrix where n and m are the lengths of
#    seq1 and seq2. The +1 accounts for the gap row and gap column.
#    Both sequences become the dimension labels of the matrix.
#    Initialize base cases: first row and first column with cumulative
#    gap penalties.
# 2. Fill each interior cell with the maximum score from three possible
#    moves: diagonal (match/mismatch), above (gap in seq1), left (gap in seq2).
#    Simultaneously record which move produced the best score in a pointer
#    matrix.
# 3. Traceback from [n_row, n_col] to [1, 1] by following the pointer matrix,
#    and directly build the two aligned sequences during traversal.



# MAIN FUNCTION ---------------------------------------------------------------
needleman_wunsch <- function(seq1, seq2,
                             match_score    =  1,
                             mismatch_score = -1,
                             gap_penalty    = -1) {
  # Step 1. Generate and initialize 2D matrix with base cases
  mat <- generate_2d_matrix(seq1, seq2, gap_penalty)
  
  # Step 2. Fill score matrix and record pointers simultaneously
  result      <- fill_matrix(mat, seq1, seq2,
                             match_score, mismatch_score, gap_penalty)
  score_mat   <- result$score
  pointer_mat <- result$pointer
  
  # Step 3. Traceback: follow pointers and build alignment strings directly
  best_alignment <- traceback(pointer_mat, seq1, seq2)
  best_alignment$score <- score_mat[nrow(score_mat), ncol(score_mat)]

  # Step 4. Print result in consistent format
  return(best_alignment)
}


# SUBFUNCTIONS ----------------------------------------------------------------

# score_pair()
# Return the score for aligning two characters.
score_pair <- function(char1, char2, match_score, mismatch_score) {
  if (char1 == char2) {
    return(match_score)
  } else {
    return(mismatch_score)
  }
}


# generate_2d_matrix()
# Create and initialize the needed matrix.
#
# matrix structure:
#   dimensions: (n_row) x (n_col) where n_row = nchar(seq2) + 1,
#               n_col = nchar(seq1) + 1. The +1 accounts for the gap
#               row and gap column (labeled ".").
#
# BASE CASES:
#   mat[1, 1] = 0           : aligning two empty sequences costs nothing
#   mat[1, j] = (j-1) * gap : aligning nothing with first j-1 chars of seq1
#                             requires j-1 gaps
#   mat[i, 1] = (i-1) * gap : aligning first i-1 chars of seq2 with nothing
#                             requires i-1 gaps

generate_2d_matrix <- function(seq1, seq2, gap_penalty) {
  n_row <- nchar(seq2) + 1
  n_col <- nchar(seq1) + 1
  
  mat <- matrix(NA, nrow = n_row, ncol = n_col)
  
  colnames(mat) <- c(".", strsplit(seq1, "")[[1]])
  rownames(mat) <- c(".", strsplit(seq2, "")[[1]])
  
  # Initialize base cases
  mat[1, ] <- seq(0, n_col - 1) * gap_penalty  # first row
  mat[, 1] <- seq(0, n_row - 1) * gap_penalty  # first column
  
  return(mat)
}


# fill_matrix()
# Fill all interior cells using the recurrence relation and simultaneously 
# record which direction produced the best score.
#
# Reccuring relation:
#   mat[i, j] = max(
#     mat[i-1, j-1] + score_pair(seq2[i-1], seq1[j-1]),  <- diagonal
#     mat[i-1, j]   + gap_penalty,                        <- above
#     mat[i,   j-1] + gap_penalty                         <- left
#   )
#
# RETURNS: list(score = ..., pointer = ...)
#   $score   : completed score matrix
#   $pointer : matrix recording "diagonal", "above", or "left" per cell
#
# TIE-BREAKING PRIORITY: diagonal > above > left
#   because when scores are equal, prefer aligning characters over
#   introducing gaps. This made biological sense because gaps are 
#   insertion/deletion mutation and therefore contain higher evolutionary cost
fill_matrix <- function(mat, seq1, seq2,
                        match_score, mismatch_score, gap_penalty) {
  
  seq1_chars <- strsplit(seq1, "")[[1]]
  seq2_chars <- strsplit(seq2, "")[[1]]
  
  n_row <- nrow(mat)
  n_col <- ncol(mat)
  
  # Initialize pointer matrix
  # NA for base case cells: traceback never visits them
  pointer_mat <- matrix(NA, nrow = n_row, ncol = n_col)
  rownames(pointer_mat) <- rownames(mat)
  colnames(pointer_mat) <- colnames(mat)
  
  # Fill interior cells only 
  for (i in 2:n_row) {
    for (j in 2:n_col) {
      
      # Characters this cell is responsible for aligning
      char_seq2 <- seq2_chars[i - 1]
      char_seq1 <- seq1_chars[j - 1]
      
      # Compute all three candidate scores
      score_diagonal <- mat[i - 1, j - 1] + score_pair(char_seq2, char_seq1,
                                                       match_score, mismatch_score)
      score_above    <- mat[i - 1, j] + gap_penalty
      score_left     <- mat[i, j - 1] + gap_penalty
      
      # Record best score and its corresponding pointer simultaneously
      best_score <- max(score_diagonal, score_above, score_left)
      mat[i, j]  <- best_score
      
      if (best_score == score_diagonal) {
        pointer_mat[i, j] <- "diagonal"
      } else if (best_score == score_above) {
        pointer_mat[i, j] <- "above"
      } else {
        pointer_mat[i, j] <- "left"
      }
    }
  }
  
  return(list(
    score   = mat,
    pointer = pointer_mat
  ))
}


# traceback()
# PURPOSE: Follow the pointer matrix from [n_row, n_col] to [1, 1] and
#          directly build the two aligned sequences during traversal.
#
# RETURNS: list(align1 = ..., align2 = ...)
#   $align1 : aligned version of seq1 (may contain "-" gaps)
#   $align2 : aligned version of seq2 (may contain "-" gaps)
#
# MOVE INTERPRETATION:
#   "diagonal" -> both sequences contribute a character (align them together)
#   "above"    -> seq2 contributes a character, seq1 contributes a gap
#   "left"     -> seq1 contributes a character, seq2 contributes a gap
#
# NOTE ON EDGE CASES:
#   When i == 1 (top row) or j == 1 (left column), pointer_mat is NA.
#   The only valid move in each case is forced: left or above respectively.

traceback <- function(pointer_mat, seq1, seq2) {
  
  seq1_chars <- strsplit(seq1, "")[[1]]
  seq2_chars <- strsplit(seq2, "")[[1]]
  
  i <- nrow(pointer_mat)
  j <- ncol(pointer_mat)
  
  # Build alignment as character vectors, prepending at each step
  align1 <- c()
  align2 <- c()
  
  while (i > 1 || j > 1) {
    
    if (i == 1) {
      # Top row: only valid move is left
      align1 <- c(seq1_chars[j - 1], align1)
      align2 <- c("-", align2)
      j <- j - 1
      
    } else if (j == 1) {
      # Left column: only valid move is above
      align1 <- c("-", align1)
      align2 <- c(seq2_chars[i - 1], align2)
      i <- i - 1
      
    } else {
      move <- pointer_mat[i, j]
      
      if (move == "diagonal") {
        align1 <- c(seq1_chars[j - 1], align1)
        align2 <- c(seq2_chars[i - 1], align2)
        i <- i - 1
        j <- j - 1
        
      } else if (move == "above") {
        align1 <- c("-", align1)
        align2 <- c(seq2_chars[i - 1], align2)
        i <- i - 1
        
      } else {
        align1 <- c(seq1_chars[j - 1], align1)
        align2 <- c("-", align2)
        j <- j - 1
      }
    }
  }
  
  return(list(
    align1 = paste(align1, collapse = ""),
    align2 = paste(align2, collapse = "")
  ))
}



