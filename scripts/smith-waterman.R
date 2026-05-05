#------------------------- SMITH-WATERMAN ALGORITHM ----------=---------------#
# Algorithm Steps:
# 1. CREATE THE MATRIX
# 2. SCORE THE MATRIX 
# 3. TRACEBACK
# 4. RETURN SEQUENCE ALIGNMENT AND SCORE
#-----------------------------------------------------------------------------#

# MAIN FUNCTION ---------------------------------------------------------------
smith_waterman = function(seq1, seq2, match = 1, mismatch = -1, gap_penalty = -1){
  mat     <- generate_2d_matrix(seq1, seq2)
  mat     <- score_matrix(mat, seq1, seq2, match, mismatch, gap_penalty)
  result  <- traceback(mat, seq1, seq2, match, mismatch, gap_penalty)
  return(result)
}

# SUBFUNCTION -----------------------------------------------------------------

# STEP 1: CREATE THE MATRIX
# Produces an (len(seq2)+1) x (len(seq1)+1) matrix filled
# with NA, with sequence characters as row/col names.
# The extra row and column are for the gap-boundary (index 0).

generate_2d_matrix = function(seq1, seq2){
  n_col <- nchar(seq1) + 1
  n_row <- nchar(seq2) + 1
  
  mat         <- matrix(NA, nrow = n_row, ncol = n_col)
  colnames(mat) <- c(".", strsplit(seq1, "")[[1]])
  rownames(mat) <- c(".", strsplit(seq2, "")[[1]])
  
  return(mat)
}

score_matrix = function(mat, seq1, seq2, match, mismatch, gap_penalty){
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



# STEP 3: TRACEBACK
# Logic:
#   (a) Start at the cell with the highest score — this is
#       where the best local alignment ENDS.
#   (b) At each cell, recompute what each of the three
#       neighbors would have contributed. Move to whichever
#       neighbor matches the current score.
#   (c) Stop when the current cell is 0 — this marks where
#       the local alignment BEGINS.
#   (d) Because we walk backwards, we prepend characters
#       rather than append, so the final strings are already
#       in the correct order.
#
# On ties: diagonal is checked first (prefer match over gap),
# then up, then left. Ties mean equally optimal alignments;
# this just picks one deterministically.


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

