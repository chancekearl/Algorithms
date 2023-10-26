#!/usr/bin/python3

from which_pyqt import PYQT_VER

if PYQT_VER == 'PYQT6':
    from PyQt6.QtCore import QLineF, QPointF
else:
    raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import math
import time
import random
import numpy as np

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1


# O(n) time. The code loops through the back pointers, puts them into a queue and then pops them off and makes the
# string to have up to 100 characters
# O(n) space complexity. The back pointers grow as the sequences grow up to align_length
def writeString(seq1, seq2, dp_matrix, num_rows, num_cols):
    operations = []
    current = dp_matrix[num_rows - 1][num_cols - 1]
    while True:
        if current.back is None:
            break
        if current.back.col < current.col:
            if current.back.row < current.row:
                operations.append("d")
            else:
                operations.append("l")
        else:
            operations.append("u")
        current = current.back
    i = 0
    while i <= 100:
        if len(operations) == 0:
            break
        op = operations.pop(-1)
        if op == 'l':
            seq1 = seq1[:i] + '-' + seq1[i:]
        elif op == 'u':
            seq2 = seq2[:i - 1] + '-' + seq2[i + 1:]
        i += 1

    return seq1, seq2


# O(n) time. The code loops through the back pointers, puts them into a queue and then pops them off and makes the
# string to have up to 100 characters. Same as above but for the banded array
# O(n) space complexity. The back pointers grow as the sequences grow up to align_length
def bandedWriteString(seq1, seq2, dp_matrix, prev_matrix, num_rows, num_cols):
    operations = []
    i = num_rows - 1
    j = MAXINDELS + len(seq2) - len(seq1)
    direction = prev_matrix[num_rows - 1][MAXINDELS + len(seq2) - len(seq1)]
    while True:
        if direction == "None":
            break
        if direction == "left":
            operations.append('l')
            direction = prev_matrix[i][j - 1]
            j -= 1
        if direction == "up":
            operations.append('u')
            if i > MAXINDELS:
                direction = prev_matrix[i - 1][j + 1]
                i -= 1
                j += 1
            else:
                direction = prev_matrix[i - 1][j]
                i -= 1
        if direction == "diag":
            operations.append('d')
            if i > MAXINDELS:
                direction = prev_matrix[i - 1][j]
                i -= 1
            else:
                direction = prev_matrix[i - 1][j - 1]
                i -= 1
                j -= 1
    q = 0
    while q <= 100:
        if len(operations) == 0:
            break
        op = operations.pop(-1)
        if op == 'l':
            seq1 = seq1[:q] + '-' + seq1[q:]
        elif op == 'u':
            seq2 = seq2[:q - 1] + '-' + seq2[q + 1:]
        q += 1

    return seq1, seq2


# Called by unbaded matrix. O(1) time complexity. O(1) space complexity. Compares the smallest cost in order of
# precedence set in the spec.
def findBackPointer(left, up, diag, a, b):
    if left is None or left == float("inf"):
        node = makeNode(100, 100)
        node.score = float("inf")
        left = node
    if up is None or up == float("inf"):
        node = makeNode(100, 100)
        node.score = float("inf")
        up = node
    if diag is None or diag == float("inf"):
        node = makeNode(100, 100)
        node.score = float("inf")
        diag = node
    back = left
    minScore = left.score
    toAdd = 5
    if up.score + 5 < minScore + 5:
        back = up
        minScore = up.score
    if a == b:
        diagVal = -3
        isMatch = True
    else:
        diagVal = 1
        isMatch = False
    if diag.score + diagVal < minScore + 5:
        back = diag
        toAdd = diagVal
    return back, isMatch, toAdd


class GeneSequencing:

    def __init__(self):
        pass

    # This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean
    # that tells you whether you should compute a banded alignment or full alignment, and _align_length_ tells you
    # how many base pairs to use in computing the alignment

    def align(self, seq1, seq2, banded, align_length):
        self.banded = banded
        self.MaxCharactersToAlign = align_length
        # calls banded matrix, O(kn) time and space. As explained above banded helper functions
        if banded:
            seq1 = seq1[:align_length]
            seq2 = seq2[:align_length]
            score, seq1, seq2 = self.isBanded(seq1, seq2)
        # calls not banded matrix, O(mn) time and space. As explained above not banded helper functions
        else:
            score, seq1, seq2 = self.notBanded(seq1, seq2)

        alignment1 = seq1[:100]
        alignment2 = seq2[:100]

        return {'align_cost': score, 'seqi_first100': alignment1, 'seqj_first100': alignment2}

    # Sets up M by N matrix filled with inf costs.
    def notBanded(self, seq1, seq2):
        num_rows = min(self.MaxCharactersToAlign, len(seq1)) + 1
        num_cols = min(self.MaxCharactersToAlign, len(seq2)) + 1
        dp_matrix = np.full([num_rows, num_cols], fill_value=float("inf"), dtype=Node)
        node = makeNode(0, 0)
        node.score = 0
        dp_matrix[0][0] = node
        # Set up base cases (first col and row)
        for i in range(1, num_rows):
            node = makeNode(i, 0)
            node.back = dp_matrix[i - 1, 0]
            node.score = INDEL * i
            dp_matrix[i][0] = node

        for j in range(1, num_cols):
            node = makeNode(0, j)
            node.back = dp_matrix[0][j - 1]
            node.score = INDEL * j
            dp_matrix[0][j] = node
        # fills in the rest of the matrix, does not revisit first row or col.
        for i in range(1, num_rows):
            for j in range(1, num_cols):
                node = makeNode(i, j)
                node.back, node.isMatch, toAdd = findBackPointer(dp_matrix[i][j - 1], dp_matrix[i - 1][j],
                                                                 dp_matrix[i - 1][j - 1], seq1[i - 1], seq2[j - 1])
                toAdd = toAdd + node.back.score
                node.score = toAdd
                dp_matrix[i][j] = node
        # Calls O(n) writeString func. (top of file)
        seq1, seq2 = writeString(seq1, seq2, dp_matrix, num_rows, num_cols)
        return dp_matrix[num_rows - 1][num_cols - 1].score, seq1[0:100], seq2[0:100]

    # Sets up N by K matrix where K = 7 filled with inf.
    def isBanded(self, seq1, seq2):
        # If the string lengths are too far apart banded will not work, returns as such
        if len(seq2) - len(seq1) > MAXINDELS:
            return float("inf"), "No Alignment Possible", "No Alignment Possible"
        num_rows = min(self.MaxCharactersToAlign, len(seq1)) + 1
        num_cols = 7
        dp_matrix = np.full([num_rows, num_cols], float("inf"))
        prev_matrix = np.full([num_rows, num_cols], "None")

        dp_matrix[0][0] = 0
        # Set up base cases according to spec.
        for i in range(1, MAXINDELS + 1):
            prev_matrix[i][0] = "up"
            dp_matrix[i][0] = i * INDEL

        for j in range(1, MAXINDELS + 1):
            prev_matrix[0][j] = "left"
            dp_matrix[0][j] = j * INDEL
        # Fills the N * K matrix with appropriate costs, does not revisit base case nodes
        for i in range(1, num_rows):
            for j in range(num_cols):
                # Fills within basecase rows where i <= 3
                if i <= MAXINDELS:
                    if j == 0 or j > i + MAXINDELS:
                        continue
                    else:
                        dp_matrix[i][j] = self.basecaseMax(i, j, dp_matrix, prev_matrix, seq1, seq2)
                # Adapts to banded size matrix with appropriate shift for back pointers.
                else:
                    if (i - MAXINDELS) + j > len(seq2):
                        continue
                    else:
                        dp_matrix[i][j] = self.bandedMax(i, j, dp_matrix, prev_matrix, seq1, seq2)

        score = dp_matrix[num_rows - 1][MAXINDELS + len(seq2) - len(seq1)]
        # Calls banded write string with O(n) complexity
        seq1, seq2 = bandedWriteString(seq1, seq2, dp_matrix, prev_matrix, num_rows, num_cols)
        return score, seq1, seq2

    # O(1) time and space complexity. Just compares two to three back pointers, depending on location of back pointers
    def bandedMax(self, i, j, dp_matrix, prev_matrix, seq1, seq2):
        offset = i - MAXINDELS
        minScore = float("inf")
        toAdd = 5
        if j > 0:
            minScore = dp_matrix[i][j - 1]
            direction = "left"
        if j < i + MAXINDELS - offset:
            if dp_matrix[i - 1][j + 1] + 5 < minScore + 5:
                minScore = dp_matrix[i - 1][j + 1]
                direction = "up"
        if seq1[i - 1] == seq2[(i - MAXINDELS) + j - 1]:
            diagVal = -3
        else:
            diagVal = 1
        if dp_matrix[i - 1][j] + diagVal < minScore + 5:
            minScore = dp_matrix[i - 1][j]
            direction = "diag"
            toAdd = diagVal
        prev_matrix[i][j] = direction
        score = minScore + toAdd
        return score

    # O(1) time and space complexity. Comparing three back values.
    def basecaseMax(self, i, j, dp_matrix, prev_matrix, seq1, seq2):
        minScore = dp_matrix[i][j - 1]
        direction = "left"
        toAdd = 5
        if dp_matrix[i - 1][j] + 5 < minScore + 5:
            minScore = dp_matrix[i - 1][j]
            direction = "up"
        if seq1[i - 1] == seq2[j - 1]:
            diagVal = -3
        else:
            diagVal = 1
        if dp_matrix[i - 1][j - 1] + diagVal < minScore + 5:
            minScore = dp_matrix[i - 1][j - 1]
            direction = "diag"
            toAdd = diagVal
        prev_matrix[i][j] = direction
        score = minScore + toAdd
        return score


# Used for unbanded. I switched for banded because it was getting too complex.
class Node:
    def __init__(self):
        self.score = float("inf")
        self.back = None
        self.row = 0
        self.col = 0
        self.isMatch = False


def makeNode(row, col):
    node = Node()
    node.row = row
    node.col = col
    return node
