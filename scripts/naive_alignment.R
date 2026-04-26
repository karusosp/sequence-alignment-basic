#-----------------------------------------------------------------------------#
#                       Naive Sequence Alignment                              #
#-----------------------------------------------------------------------------#
# The naive sequence alignment is a function that take two sequences and 
# find the best alignment by way of brute-force: finding each possible alignment
# and determine the best one based on a scoring system.
# 
# The general algorithm can be described in several steps: 
# 1. Generate two dimensional m x n matrix where m and n is the length of 
#    the first and second sequence sequentially, and both sequences
#    become the dimension name of the matrix
# 2. Find all the possible paths from the starting position (0,0) within the 
#    matrix to the end position (m,n)
# 3. Determine from each completed path all the possible sequence alignment
# 4. Give a score for each sequence alignment and find the best one
# 5. Return the best alignment, its score, and the total possible alignment

# MAIN FUNCTION: naive_algorithm() --------------------------------------------
naive_alignment <- function(seq1, seq2, 
                            match_score = 1, 
                            mismatch_score = -1, 
                            gap_score = -1) {
  
  # Step 1. Generate matrix and its Dimension Name
  mat <- generate_2d_matrix(seq1, seq2)
  nrow <- nrow(mat)
  ncol <- ncol(mat)
  
  # Step 2. Generate all possible paths
  all_paths <- generate_paths(1, 1, nrow, ncol)
  total_paths <- length(all_paths)
  best_score <- -Inf
  best_alignments <- list()
  # Step 3. Determine all possible path from all paths
  for (path in all_paths) {
    # Convert path to alignment strings
    alignment <- path_to_alignment(path, seq1, seq2)
    
 # Step 4. Score each alignment and find the best one

    
    # SCORING PARAMETER 
    match_score = match_score
    mismatch_score = mismatch_score
    gap_score = gap_score
    
    score <- score_alignment(alignment$align1, alignment$align2, 
                             match_score, mismatch_score, gap_score)
    
    # Track the best scoring alignment
    if (score > best_score) {
      best_score <- score
      best_alignments <- list(alignment)
    } else if (score == best_score) {
      # Keep track of ties too
      best_alignments[[length(best_alignments) + 1]] <- alignment
    }
  }
  
  # Return the results: total possible alignment and the best alignment 
  cat("Total number of possible alignments:", total_paths, "\n\n")
  cat("BEST ALIGNMENT(S) with score:", best_score, "\n")
  for (aln in best_alignments) {
    cat(aln$align1, "\n")
    cat(aln$align2, "\n\n")
  }
}
# SUBFUNCTION ------------------------------------------------------------------

# generate_2d_matrix()
# function: generate 2 dimensional matrix with the size of m x n, where 
#           m and n are the length of the first and second sequence.
# input: two sequence (seq1, and seq2)
# output: two dimensional matrix with the sequence as the dimension name

generate_2d_matrix <- function(seq1, seq2) {
  nrow <- nchar(seq2) + 1  # +1 for the gap row
  ncol <- nchar(seq1) + 1  # +1 for the gap column
  
  mat <- matrix(NA, nrow = nrow, ncol = ncol)
  
  # Label columns with seq1 characters, rows with seq2 characters
  colnames(mat) <- c(".", strsplit(seq1, "")[[1]])
  rownames(mat) <- c(".", strsplit(seq2, "")[[1]])
  
  return(mat)
}


# generate_paths
# function: determine all the possible path from the starting position (1,1)
#           to end position (nrow, ncol, where each path is a combination of
#           three legal moves: "D" (diagonal), "R" (right), "Down" (down).
#           The function is invoked inside the main naive_algorithm() function
#
# desc: i = current row, j = current column
#           We start at (1,1) and want to reach (nrow, ncol) by recursion
generate_paths <- function(i, j, nrow, ncol, 
                           current_path = c(), 
                           all_paths = list()) {
  
  # BASE CASE: reached bottom-right corner
  # Save the completed path and return
  if (i == nrow && j == ncol) {
    all_paths[[length(all_paths) + 1]] <- current_path
    return(all_paths)
  }
  
  # BOUNDARY CASE: at last row, can only move right
  if (i == nrow) {
    all_paths <- generate_paths(i, j + 1, nrow, ncol,
                                c(current_path, "R"),
                                all_paths)
    return(all_paths)
  }
  
  # BOUNDARY CASE: at last column, can only move down
  if (j == ncol) {
    all_paths <- generate_paths(i + 1, j, nrow, ncol,
                                c(current_path, "Down"),
                                all_paths)
    return(all_paths)
  }
  
  # GENERAL CASE: try all three moves
  # Move diagonally (consume from both sequences)
  all_paths <- generate_paths(i + 1, j + 1, nrow, ncol,
                              c(current_path, "D"),
                              all_paths)
  # Move right (consume from seq1, gap in seq2)
  all_paths <- generate_paths(i, j + 1, nrow, ncol,
                              c(current_path, "R"),
                              all_paths)
  # Move down (gap in seq1, consume from seq2)
  all_paths <- generate_paths(i + 1, j, nrow, ncol,
                              c(current_path, "Down"),
                              all_paths)
  
  return(all_paths)
}


# path_to_alignment()
# function: from each possible path generated, determine the alignment based
#           on the vector moves inside the path. The logic is that diagonal
#           move represent a match/mismatch and therefore exhaust both sequence
#           length (i + 1, j + 1), whereas downward and rightward move 
#           represent gap and thus exhaust only one of the sequence length
path_to_alignment <- function(path, seq1, seq2) {
  seq1_chars <- strsplit(seq1, "")[[1]]  # e.g. c("A", "T")
  seq2_chars <- strsplit(seq2, "")[[1]]  # e.g. c("G", "A")
  
  align1 <- c()  # will build up seq1's alignment string
  align2 <- c()  # will build up seq2's alignment string
  
  i <- 0  # tracks position in seq2
  j <- 0  # tracks position in seq1
  
  for (move in path) {
    if (move == "D") {
      # Diagonal: consume one char from each sequence
      i <- i + 1
      j <- j + 1
      align1 <- c(align1, seq1_chars[j])
      align2 <- c(align2, seq2_chars[i])
      
    } else if (move == "R") {
      # Right: consume from seq1, insert gap in seq2
      j <- j + 1
      align1 <- c(align1, seq1_chars[j])
      align2 <- c(align2, "-")
      
    } else if (move == "Down") {
      # Down: insert gap in seq1, consume from seq2
      i <- i + 1
      align1 <- c(align1, "-")
      align2 <- c(align2, seq2_chars[i])
    }
  }
  
  # Collapse character vectors into strings
  return(list(
    align1 = paste(align1, collapse = ""),
    align2 = paste(align2, collapse = "")
  ))
}

# score_alignment()
# function: determine the score for each alignment by a scoring system
# Simple scoring: match = +1, mismatch = -1, gap = -1

score_alignment <- function(align1, align2, 
                            match_score, 
                            mismatch_score, 
                            gap_score) {
  
  chars1 <- strsplit(align1, "")[[1]]
  chars2 <- strsplit(align2, "")[[1]]
  
  total_score <- 0
  
  for (k in seq_along(chars1)) {
    if (chars1[k] == "-" || chars2[k] == "-") {
      # Gap in either sequence
      total_score <- total_score + gap_score
    } else if (chars1[k] == chars2[k]) {
      # Match
      total_score <- total_score + match_score
    } else {
      # Mismatch
      total_score <- total_score + mismatch_score
    }
  }
  
  return(total_score)
}