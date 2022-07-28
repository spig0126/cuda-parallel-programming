import scipy.sparse as sparse
import scipy.stats as stats
import numpy as np
from numpy.random import default_rng
from typing import List
from io import TextIOWrapper
import sys

np.set_printoptions(threshold=sys.maxsize)


def generateSparseMatrixArray(w: int, h: int, d: int) -> List[int]:
    np.random.seed(100)
    rng = default_rng()
    rvs = stats.poisson(15, loc=10).rvs  # generate random numbers from Poisson distribution
    sparseMatrix = sparse.random(w, h, density=d, random_state=rng, data_rvs=rvs)  # create sparse matrix

    a = sparseMatrix.toarray()  # change sparse matrix to array
    a = a.astype(int)
    printMatrix(a.tolist())  # change array to list for further easier use
    return a.tolist()


def printMatrix(l: List[List[int]]):
    print("[", end="")
    for row, row_i in zip(l, range(len(l))):
        print(row, end="]\n\n" if row_i == (len(l) - 1) else "\n")


def printVector(l: List[int]):
    print(l)

def printSparseMatrixList(l):
    print("[")
    for row in l:
        print(row)
    print("]")


def saveData(f: TextIOWrapper, data: list, flag: bool, varname: str) -> None:
    if flag:  # if 1d vector
        f.write("int {}[{}] = ".format(varname, len(data)))
        f.write("{")
        f.write(', '.join(str(e) for e in data))
        f.write("};\n\n")

    else:  # 2d matrix
        f.write("int {}[{}][{}] = ".format(varname, len(data), len(data[0])))
        f.write("{")
        for row, row_i in zip(data, range(len(data))):
            f.write("{")
            f.write(', '.join(str(e) for e in row))
            f.write("}," if row_i != (len(data) - 1) else "}")
        f.write("};\n\n")


def COO(l):
    data = []
    col_index = []
    row_index = []
    for row, r in zip(l, range(len(l))):
        for e, c in zip(row, range(len(row))):
            if e != 0:
                data.append(e)
                col_index.append(c)
                row_index.append(r)
    print("\n----------COO-------------")
    print("data")
    printVector(data)
    print("col_index")
    printVector(col_index)
    print("row_index")
    printVector(row_index)
    return data, col_index, row_index


def CSR(l):
    data = []
    col_index = []
    row_ptr = []

    row_ptr.append(0)
    for row in l:
        for e in row:
            if e != 0:
                data.append(e)
                col_index.append(row.index(e))
        row_ptr.append(len(data))

    print()
    print("-----------------CSR----------------")
    print("# of nonzero elements", len(data))
    print("data")
    printVector(data)
    print("col_index")
    printVector(col_index)
    print("row_ptr")
    printVector(row_ptr)

    return data, col_index, row_ptr


def ELL(l):
    max = 0
    maxRow = 0
    for row in l:
        cnt = 0
        for e in row:
            if e != 0:
                cnt += 1
        if cnt > max:
            max = cnt
            maxRow = l.index(row)

    l_padding = [[] for _ in range(len(l))]
    l_padding_colIndex = [[] for _ in range(len(l))]
    for row, i in zip(l, range(len(l))):
        cnt = 0
        for e, c in zip(row, range(len(row))):
            if e != 0:
                l_padding[i].append(e)
                l_padding_colIndex[i].append(c)
                cnt += 1
        if cnt < max:
            for _ in range(max - cnt):
                l_padding[i].append(0)
                l_padding_colIndex[i].append(0)

    l_transposed = np.transpose(np.asarray(l_padding))
    l_colIndex_transposed = np.transpose(np.asarray(l_padding_colIndex))
    print()
    print("-----------------ELL----------------")
    print("padding")
    printMatrix(np.asarray(l_padding).tolist())
    print("padding col_index")
    printMatrix(np.asarray(l_padding_colIndex).tolist())

    print("transposed")
    printMatrix(l_transposed.tolist())
    print("transposed col_index")
    printMatrix(l_colIndex_transposed.tolist())
    return l_transposed, l_colIndex_transposed


if __name__ == '__main__':
    f = open("sparseMatrix.h", "a")

    smList = generateSparseMatrixArray(32, 32, 1 / 32)
    saveData(f, smList, False, "A")

    data, col_index, row_index = COO(smList)
    saveData(f, data, True, "data")
    saveData(f, col_index, True, "col_index")
    saveData(f, row_index, True, "row_index")

    data, col_index, row_ptr = CSR(smList)
    saveData(f, row_ptr, True, "row_ptr")

    l_transposed, l_colIndex_transposed = ELL(smList)
    saveData(f, l_transposed, False, "transposed")
    saveData(f, l_colIndex_transposed, False, "transposed_col_index")