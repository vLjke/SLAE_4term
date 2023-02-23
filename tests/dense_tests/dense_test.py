import numpy as np
from random import randint

# matrix multiply vector test 1
# 100x100 random matrix
m1 = np.random.random((100, 100)) * 3
# 100x1 random vector
v1 = np.random.random((100, 1)) * 3
# result of multiplying
r1 = m1.dot(v1)
with open('Dense_matrix_mult_vector_1.txt', 'w') as f:
    # write matrix dimensions in the first line
    f.write(str(np.shape(m1)[0]) + ' ' + str(np.shape(m1)[1]) + '\n')
    # write matrix in the next lines
    for i in m1:
        for j in i:
            f.write(str(j) + ' ')
        f.write('\n')
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

# matrix multiply vector test 2
# 30x50 random matrix
m2 = np.random.random((30, 50)) * 5
# 50x1 random vector
v2 = np.random.random((50, 1)) * 5
# result of multiplying
r2 = m2.dot(v2)
with open('Dense_matrix_mult_vector_2.txt', 'w') as f:
    # write matrix dimensions in the first line
    f.write(str(np.shape(m2)[0]) + ' ' + str(np.shape(m2)[1]) + '\n')
    # write matrix in the next lines
    for i in m2:
        for j in i:
            f.write(str(j) + ' ')
        f.write('\n')
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

# matrix get element test
# 10x10 random matrix
m3 = np.random.random((10, 10)) * 2
# number of elements to get
N = 10
with open('Dense_matrix_get_elem.txt', 'w') as f:
    # write matrix dimensions in the first line
    f.write(str(np.shape(m3)[0]) + ' ' + str(np.shape(m3)[1]) + '\n')
    # write matrix in the next lines
    for i in m3:
        for j in i:
            f.write(str(j) + ' ')
        f.write('\n')
    # write number of elements to get
    f.write(str(N) + '\n')
    # write indexes and value of elements' to get
    for k in range(N):
        i, j = randint(0, 9), randint(0, 9)
        f.write(str(i) + ' ' + str(j) + ' ' + str(m3[i][j]) + '\n')
f.close()

