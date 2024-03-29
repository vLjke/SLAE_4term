import numpy as np
from random import randint

# DOK cell sort test
# number of cells to sort
N = 10
with open('DOK_cell_sort.txt', 'w') as f:
    # write number of cells to get
    f.write(str(N) + '\n')
    for _ in range(N):
        i, j, val = randint(0, 100), randint(0, 100), randint(0, 100)
        f.write(str(i) + ' ' + str(j) + ' ' + str(val) + '\n')
f.close()


# CSR matrix multiply vector test 1
# 20x20 random matrix
m1 = np.random.random((20, 20)) * 10
# 20x1 random vector
v1 = np.random.random((20, 1)) * 3
# add some zeros
for _ in range(50):
    i, j = randint(0, 19), randint(0, 19)
    m1[i][j] = 0
# number of non-zero elements
N = np.count_nonzero(m1)
# result of multiplying
r1 = m1.dot(v1)
with open('CSR_matrix_mult_vector_1.txt', 'w') as f:
    # write matrix dimensions
    f.write(str(np.shape(m1)[0]) + ' ' + str(np.shape(m1)[1]) + '\n')
    # write number of non-zero elements
    f.write(str(N) + '\n')
    # write matrix in the next lines
    for i in range(np.shape(m1)[0]):
        for j in range(np.shape(m1)[1]):
            if m1[i][j] != 0:
                f.write(str(i) + ' ' + str(j) + ' ' + str(m1[i][j]) + '\n')
    # write vector in the next line
    for v in v1:
        for j in v:
            f.write(str(j) + ' ')
    f.write('\n')
    # write result of multiplying in the last line
    for r in r1:
        for j in r:
            f.write(str(j) + ' ')
    f.write('\n')
f.close()


# CSR matrix multiply vector test 2
# 200x500 random matrix
m2 = np.random.random((100, 200)) * 3
# 200x1 random vector
v2 = np.random.random((200, 1)) * 3
# add some zeros
for _ in range(5000):
    i, j = randint(0, 99), randint(0, 199)
    m2[i][j] = 0
# number of non-zero elements
N = np.count_nonzero(m2)
# result of multiplying
r2 = m2.dot(v2)
with open('CSR_matrix_mult_vector_2.txt', 'w') as f:
    # write matrix dimensions
    f.write(str(np.shape(m2)[0]) + ' ' + str(np.shape(m2)[1]) + '\n')
    # write number of non-zero elements
    f.write(str(N) + '\n')
    # write matrix in the next lines
    for i in range(np.shape(m2)[0]):
        for j in range(np.shape(m2)[1]):
            if m2[i][j] != 0:
                f.write(str(i) + ' ' + str(j) + ' ' + str(m2[i][j]) + '\n')
    # write vector in the next line
    for v in v2:
        for j in v:
            f.write(str(j) + ' ')
    f.write('\n')
    # write result of multiplying in the last line
    for r in r2:
        for j in r:
            f.write(str(j) + ' ')
    f.write('\n')
f.close()


# CSR matrix get n elements test
# 20x20 random matrix
m3 = np.random.random((20, 20)) * 4
# 20x1 random vector
v3 = np.random.random((20, 1)) * 5
# add some zeros
for _ in range(100):
    i, j = randint(0, 19), randint(0, 19)
    m3[i][j] = 0
# number of non-zero elements
N = np.count_nonzero(m3)
# number of elements to get
n = 10
with open('CSR_matrix_get_elems.txt', 'w') as f:
    # write matrix dimensions
    f.write(str(np.shape(m3)[0]) + ' ' + str(np.shape(m3)[1]) + '\n')
    # write number of non-zero elements
    f.write(str(N) + '\n')
    # write matrix in the next lines
    for i in range(np.shape(m3)[0]):
        for j in range(np.shape(m3)[1]):
            if m3[i][j] != 0:
                f.write(str(i) + ' ' + str(j) + ' ' + str(m3[i][j]) + '\n')
    # write number of elements to get
    f.write(str(n) + '\n')
    # write indexes and value of elements' to get
    for _ in range(n):
        i, j = randint(0, 19), randint(0, 19)
        f.write(str(i) + ' ' + str(j) + ' ' + str(m3[i][j]) + '\n')
f.close()


# CSR matrix GMRES method test
# 20x20 random matrix
m4 = np.random.random((20, 20)) * 4 + np.eye(20) * 8
# 20x1 random vector
v4 = np.random.random((20, 1)) * 5
# add some zeros
for _ in range(300):
    i, j = randint(0, 19), randint(0, 19)
    m4[i][j] = 0
# number of non-zero elements
N = np.count_nonzero(m4)
# SLAE solution
sol = np.linalg.solve(m4, v4)
with open('CSR_matrix_GMRES(n)_method_test.txt', 'w') as f:
    # write matrix dimensions
    f.write(str(np.shape(m4)[0]) + ' ' + str(np.shape(m4)[1]) + '\n')   
    # write number of non-zero elements
    f.write(str(N) + '\n')
    # write matrix in the next lines
    for i in range(np.shape(m4)[0]):
        for j in range(np.shape(m4)[1]):
            if m4[i][j] != 0:
                f.write(str(i) + ' ' + str(j) + ' ' + str(m4[i][j]) + '\n')
    for v in v4:
        for j in v:
            f.write(str(j) + ' ')
    f.write('\n')
    for s in sol:
        for j in s:
            f.write(str(j) + ' ')
    f.write('\n')
f.close()


# CSR transposed matrix test
# 50x20 random matrix
m6 = np.random.random((50, 20)) * 4
# add some zeros
for _ in range(500):
    i, j = randint(0, 49), randint(0, 19)
    m6[i][j] = 0
# transposed matrix
m6_T = m6.T
# number of non-zero elements
N = np.count_nonzero(m6)
with open('CSR_transposed_matrix_test.txt', 'w') as f:
    # write matrix dimensions
    f.write(str(np.shape(m6)[0]) + ' ' + str(np.shape(m6)[1]) + '\n')
    # write number of non-zero elements
    f.write(str(N) + '\n')
    # write matrix in the next lines
    for i in range(np.shape(m6)[0]):
        for j in range(np.shape(m6)[1]):
            if m6[i][j] != 0:
                f.write(str(i) + ' ' + str(j) + ' ' + str(m6[i][j]) + '\n')
    # write transposed matrix in the next lines
    # write matrix dimensions
    f.write(str(np.shape(m6_T)[0]) + ' ' + str(np.shape(m6_T)[1]) + '\n')
    # write number of non-zero elements
    f.write(str(N) + '\n')
    # write matrix in the next lines
    for i in range(np.shape(m6_T)[0]):
        for j in range(np.shape(m6_T)[1]):
            if m6_T[i][j] != 0:
                f.write(str(i) + ' ' + str(j) + ' ' + str(m6_T[i][j]) + '\n')

f.close()


# CSR matrix BiCG method test
# 20x20 random matrix
m5 = np.random.random((20, 20)) * 4 + np.eye(20) * 4
# 20x1 random vector
v5 = np.random.random((20, 1)) * 5
# add some zeros
for _ in range(300):
    i, j = randint(0, 19), randint(0, 19)
    m5[i][j] = 0
# number of non-zero elements
N = np.count_nonzero(m5)
# SLAE solution
sol = np.linalg.solve(m5, v5)
with open('CSR_matrix_BiCG_method_test.txt', 'w') as f:
    # write matrix dimensions
    f.write(str(np.shape(m5)[0]) + ' ' + str(np.shape(m5)[1]) + '\n')
    # write number of non-zero elements
    f.write(str(N) + '\n')
    # write matrix in the next lines
    for i in range(np.shape(m5)[0]):
        for j in range(np.shape(m5)[1]):
            if m5[i][j] != 0:
                f.write(str(i) + ' ' + str(j) + ' ' + str(m5[i][j]) + '\n')
    for v in v5:
        for j in v:
            f.write(str(j) + ' ')
    f.write('\n')
    for s in sol:
        for j in s:
            f.write(str(j) + ' ')
    f.write('\n')
f.close()

