from enum import Enum, auto
import numpy as np
from dataclasses import dataclass
import copy
import argparse

# Construct grid
def initialize_matrix(reference, query):
    return np.full((len(reference) + 1, len(query) + 1), np.nan)

# Scoring system options
class Option(Enum):
    BASIC = auto() # Most basic scoring system
    LIBRE = auto() # Scoring system from LIBREtexts

# Select scoring system
def retrieve_scoring_system(option: Option):
    match option:
        case option.BASIC:
            return {"match" : 1, "mismatch" : -1, "indel" : -1}
        case option.LIBRE:
            return {"match" : 1, "mismatch" : -1, "indel" : -2}
        case _:
            raise ValueError("Invalid scoring system selected")

# Fill in the matrix
def fill_matrix(matrix, scoring_system, reference, query):
    # Fill matrix according to scoring system
    # reference = rows, query = columns (stay consistent)
    rows, columns = matrix.shape

    # Initialize row and column gap penalties
    match_s, mismatch_s, indel_s = scoring_system["match"], scoring_system["mismatch"], scoring_system["indel"]

    # Initialize first row / column
    matrix[0,:] = np.arange(0, indel_s * columns, indel_s)
    matrix[:,0] = np.arange(0, indel_s * rows, indel_s)    

    # fill matrix
    for row in range(1, rows):
        for column in range(1, columns):
            # Solve for 3 canidate scores
            # Select highest score, update cell accordingly
            top_s = matrix[row-1, column] + indel_s
            left_s = matrix[row, column - 1] + indel_s
            if query[column - 1] == reference[row - 1]:
                topleft_s = matrix[row - 1, column - 1] + match_s
            else:
                topleft_s = matrix[row - 1, column - 1] + mismatch_s

            matrix[row, column] = max(top_s, left_s, topleft_s)

    return matrix

# Trace back to the origin ( Retrieve optimal global alignment - singular)
def retrieve_alignment(matrix, reference, query, scoring_system):
    
    q_aln = ""
    r_aln = ""

    row, column = matrix.shape
    row -= 1
    column -= 1
    score = matrix[row, column]

    while not (row == 0 and column == 0):
        # Select which cell to pivot to, update ptrs / sequences accordingly

        # If direction from vertical
        if matrix[row - 1, column] == score - scoring_system["indel"]:
            q_aln = "_" + q_aln
            r_aln = reference[row - 1] + r_aln
            row = row - 1
            score = matrix[row, column]
        # If from horizontal
        elif matrix [row, column - 1] == score - scoring_system["indel"]:
            q_aln = query[column - 1] + q_aln
            r_aln = "_" + r_aln
            column = column - 1
            score = matrix[row, column]

        # If from match / mismatch
        else:
            q_aln = query[column - 1] + q_aln
            r_aln = reference[row - 1] + r_aln
            row -= 1
            column -= 1
            score = matrix[row, column]

    return [r_aln], [q_aln]

@dataclass
class Alignment:
    q_aln: str
    r_aln: str
    row: int
    column: int
    score: int

# Update the alignment based on the matrix and the path
def update_alignment(alignment, matrix, reference, query, path):
    col = alignment.column
    row = alignment.row

    match path:
        case 'vertical':
            alignment.q_aln = '_' + alignment.q_aln
            alignment.r_aln = reference[ row - 1 ] + alignment.r_aln
            alignment.row -= 1
             
        case 'diagonal':
            alignment.q_aln = query[ col - 1 ] + alignment.q_aln
            alignment.r_aln = reference[ row - 1 ] + alignment.r_aln
            alignment.column -= 1
            alignment.row -= 1

        case 'left':
            alignment.q_aln = query[ col - 1 ] + alignment.q_aln
            alignment.r_aln = '_' + alignment.r_aln
            alignment.column -= 1
    

    alignment.score = matrix[alignment.row, alignment.column]

def retrieve_all_alignments(matrix, reference, query, scoring_system):
    
    q_aln = "" # List
    r_aln = "" # List

    row, column = matrix.shape
    row -= 1 
    column -= 1
    score = matrix[row, column]

    active_alignments = [Alignment("", "", row, column, score)]
    completed_alignments = []

    # Increment over alignments - while not alignments complete
    all_alignments_retrieved = False
    while not all_alignments_retrieved:
        update_q = []
        for alignment in active_alignments:
            row = alignment.row
            column = alignment.column
            score = alignment.score

            # Locate branching alignments / store
            valid_paths = []
            if matrix[row - 1, column] == score - scoring_system["indel"]:
                valid_paths.append("vertical")
            if matrix[row, column - 1] == score - scoring_system["indel"]:
                valid_paths.append("left")
            # If match
            if reference[ row - 1 ] == query[ column - 1 ]:
                if matrix[row - 1, column - 1] == score - scoring_system["match"]:
                    valid_paths.append("diagonal")
            # If mismatch
            else:
                if matrix[row - 1, column - 1] == score - scoring_system["mismatch"]:
                    valid_paths.append("diagonal")

            # Copy alignment as a template
            alignment_template = copy.deepcopy(alignment)

            # Update first alignment
            update_alignment(alignment, matrix, reference, query, valid_paths[0])
            del valid_paths[0]

            # For next alignments (if next alignments), copy from copy and update each w/ corresponding path - append to q
            for path in valid_paths:
                branching_alignment = copy.deepcopy(alignment_template)
                update_alignment(branching_alignment, matrix, reference, query, path)
                update_q.append(branching_alignment)


        # Iterate over all alignments and remove completed alignments
        del_q = []
        for i, alignment in enumerate(active_alignments):
            if alignment.row == 0 and alignment.column == 0:
                completed_alignments.append(alignment)
                del_q.append(i)
        for i in del_q[::-1]: # Do in reverse order so I don't need to worry about reordering
            del active_alignments[i]

        # Add queue alignments to active alignments
        for alignment in update_q:
            active_alignments.append(alignment)
        update_q = []
    
        if len(active_alignments) == 0:
            all_alignments_retrieved = True

    # Extract the alignments
    references = [alignment.r_aln for alignment in completed_alignments]
    queries = [alignment.q_aln for alignment in completed_alignments]

    return references, queries

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", help = "Reference sequence (rows)")
    parser.add_argument("--query", help = "Query sequence (columns)")
    parser.add_argument("--scoring_system", choices = ["basic", "libre"], help = "Scoring system to use")
    parser.add_argument("--retrieve_all", action = 'store_true', default = False, help = "Retrieve all alignments")
    args = parser.parse_args()

    #reference = "GATTACA"
    #query = "GTCGACGCA"

    matrix = initialize_matrix(args.reference, args.query)

    #ss = retrieve_scoring_system(Option.BASIC)
    option_dict = {"basic" : Option.BASIC, "libre" : Option.LIBRE}
    option = option_dict[args.scoring_system]
    ss = retrieve_scoring_system(option)

    filled_matrix = fill_matrix(matrix, ss, args.reference, args.query)

    if args.retrieve_all:
        references, queries = retrieve_all_alignments(filled_matrix, args.reference, args.query, ss)
    else:
        references, queries = retrieve_alignment(filled_matrix, args.reference, args.query, ss)

    print("Alignment(s)")
    for r, q in zip(references, queries):
        print("--")
        print(f"{q}")
        print(f"{r}")
