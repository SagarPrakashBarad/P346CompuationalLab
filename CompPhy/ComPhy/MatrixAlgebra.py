global eps
import numpy as np
eps = 1/10**8 # Predefined epsilon, will be used for tolerance testing.

class Matrix:

    def __init__(self, dims, fill):

        self.rows = dims[1]
        self.cols = dims[0]
        self.is_symmetric = True
        self.is_dd = True

        if isinstance(fill,list):
            k = 0
            lst = []
            for i in range(dims[0]):
                col = []
                for j in range(dims[1]):
                    col.append(fill[k])
                    k += 1
                #col.append(x[i])
                lst.append(col)
            self.M = lst
        else:
            self.M = [[fill] * self.cols for i in range(self.rows)]

        for i in range(self.rows):
            for j in range(self.cols):
                if self.M[i][j] != self.M[j][i]:
                    self.is_symmetric = True



    def __str__(self):
        m = len(self.M)
        mtxStr = ''

        for i in range(self.rows):
            mtxStr += '-------'

        mtxStr += '------\n'
        for i in range(m):
            mtxStr += ('|' + ', '.join( map(lambda x: '{0:8.3f}'.format(x), self.M[i])) + '| \n')

        for i in range(self.rows):
            mtxStr += '-------'
        mtxStr += '------\n'

        return mtxStr

    def jacobi_iter(self,B,init_val):

        self.var = init_val
        iter = 0
        tol = init_val.copy()
        t = max(tol)

        print("System of Linear Equations:\n")
        for i in range(len(self.M[0])):
            row = ["{}x{}".format(self.M[i][j], j) for j in range(len(self.M[0]))]
            print(" + ".join(row), "=", B[i])


        while t > eps:
            for i in range(len(self.M[0])):
                sum = 0
                for j in range(len(self.M[1])):
                    if i != j:
                        sum += self.M[i][j]*self.var[j]
                vb = self.var[i]
                self.var[i] = (1/self.M[i][i])*(B[i] - sum)
                tol[i] = abs((self.var[i] - vb)*100/vb)

            print("Iteration: {} | Current Solution: {}".format(iter, self.var))
            iter = iter + 1
            t = max(tol)

        print("\nSolution converged, after {} iterations.".format(iter))

    def gauss_siedel(self,B,init_val):

        self.var = init_val
        self.var_new = self.var.copy()
        iter = 0
        tol = init_val.copy()
        t = max(tol)

        print("System of Linear Equations:\n")
        for i in range(len(self.M[0])):
                row = ["{}x{}".format(self.M[i][j], j) for j in range(len(self.M[0]))]
                print(" + ".join(row), "=", B[i])

        while t > eps:
            for i in range(len(self.M[0])):
                l = 0
                u = 0
                for j in range(len(self.M[1])):
                    if j < i:
                        l += self.M[i][j]*self.var_new[j]
                    elif j > i:
                        u += self.M[i][j]*self.var[j]
                vb = self.var_new[i]
                self.var_new[i] = (1/self.M[i][i])*(B[i] - l - u)
                tol[i] = abs((self.var[i] - vb)*100/vb)

            self.var = self.var_new.copy()

            print("Iteration: {} | Current Solution: {}".format(iter, self.var_new))
            iter = iter + 1
            t = max(tol)

        print("\nSolution converged, after {} iterations.".format(iter))


    def forward_substitution(self,B):

        self.y = []

        for i in range(len(self.M[0])):
            sum = 0
            if self.M[i][i] == 0: continue # skips singular pivot entries to avoid division by zero
            for j in range(i): # for all off daigonal elements
                sum += self.M[i][j]*self.y[j] # take sum of product non-pivot elements with corresponding y value
            self.y.append((B[i] - sum)/self.M[i][i])

    def backward_substitution(self):

        self.x = self.y.copy()

        self.x[self.rows -1] = [self.rows -1]/self.M[self.rows -1][self.cols -1]
        for i in range(self.rows - 1,-1,-1):
            sum = 0
            if self.M[i][i] == 0: continue # skips singular pivot entries to avoid division by zero
            for j in range(i+1,self.rows):
                sum += self.M[i][j]*self.x[j] # take sum of product non-pivot elements with corresponding y value
            self.x[i] = (1/self.M[i][i])*(self.y[i] - sum)

    def LUdecomp(self):

        (n,m) = (self.rows, self.cols)
        self.L = [[0] * m for i in range(n)]
        self.U = [[0] * m for i in range(n)]

        for i in range(n):
            max = i
            if abs(self.M[i][i]) <= eps:
                for j in range(i+1,n):
                    if abs(self.M[j][i]) > abs(self.M[max][i]):
                            max = j
                    (self.M[i], self.M[max]) = (self.M[max], self.M[i])

            for k in range(i,n):
                sum = 0
                for j in range(i):
                    sum += self.L[i][j]*self.U[j][k]
                self.U[i][k] = self.M[i][k] - sum

            for l in range(i,n):
                if i == l:
                    self.L[i][i] = 1
                elif self.U[i][i] != 0:
                    sum = 0
                    for j in range(i):
                        sum += self.L[l][j]*self.U[j][i]

                    self.L[l][i] = (self.M[l][i] - sum)/self.M[i][i]

    def Cholesky(self):

        (n,m) = (self.rows, self.cols)
        self.L = [[0] * m for i in range(n)]
        self.U = [[0] * m for i in range(n)]

def DDcheck(A,i): # checks whether a pivot is daigonally dominant
        sum = 0
        for j in range(len(A.M[0])):
            if i != j:
                sum += A.M[i][j]
            sum = abs(sum)
        if abs(A.M[i][i]) < sum:
             return 1


def makeDD(mat,b):
        (n, m) = (len(mat.M), len(mat.M[0]))

        for i in range(0,n):
         maxelement = i
         if DDcheck(mat,i) == 1: # if pivot is not daigonally dominant
                for j in range(n):
                    if abs(mat.M[j][i]) > abs(mat.M[maxelement][i]):
                     maxelement = j
                (mat.M[i], mat.M[maxelement]) = (mat.M[maxelement], mat.M[i]) # swap places with row having largest column element
                (b[i], b[maxelement]) = (b[maxelement], b[i]) # swaping to simualte an augemntated matrix

def gaussjordan(mat,b):
    (n, m) = (len(mat.M), len(mat.M[0]))
    """gaussjordan takes a matrix (2D array) as input and reduce it to Reduced Row Echelon Form.
        Returns True if successful, False if 'm' is singular."""

    # Main loop
    for i in range(0,n):
        # partial pivoting
        maxelement = i
        if abs(mat.M[i][i]) <= eps:
            for j in range(i+1,n):
                if abs(mat.M[j][i]) > abs(mat.M[maxelement][i]):
                    maxelement = j
            (mat.M[i], mat.M[maxelement]) = (mat.M[maxelement], mat.M[i])
            (b[i], b[maxelement]) = (b[maxelement], b[i])

        # Division of pivot row
        pivot = mat.M[i][i]
        for k in range(i,n):
            mat.M[i][k] /= pivot
        b[i] /= pivot

        # Elimination loop
        for j in range(n):
            if i == j or mat.M[j][i] == 0: continue

            ratio = mat.M[j][i]
            for k in range(i,n):
                mat.M[j][k] -= ratio*mat.M[i][k]
            b[j] -= ratio*b[i]

    return b, mat



def lu_decomposition(mat):
    """Perform LU decomposition of a square matrix."""
    n = len(mat)
    L = np.zeros((n, n))
    U = np.zeros((n, n))

    for i in range(n):
        for j in range(i, n):
            U[i][j] = mat[i][j] - sum(L[i][k] * U[k][j] for k in range(i))
            if i == j:
                L[i][i] = 1
            else:
                L[j][i] = (mat[j][i] - sum(L[j][k] * U[k][i] for k in range(i))) / U[i][i]

    return L, U

def forward_substitution(L, b):
    """Perform forward substitution to solve Ly = b."""
    n = len(b)
    y = np.zeros(n)

    for i in range(n):
        y[i] = b[i] - sum(L[i][j] * y[j] for j in range(i))

    return y

def backward_substitution(U, y):
    """Perform backward substitution to solve Ux = y."""
    n = len(y)
    x = np.zeros(n)

    for i in range(n - 1, -1, -1):
        x[i] = (y[i] - sum(U[i][j] * x[j] for j in range(i + 1, n))) / U[i][i]

    return x

def solve_system_lu_decomposition(A, b):
    """Solve a system of linear equations Ax = b using LU decomposition."""
    L, U = lu_decomposition(A)
    y = forward_substitution(L, b)
    x = backward_substitution(U, y)
    return x

def inverse_matrix_lu_decomposition(A):
    """Find the inverse of a matrix using LU decomposition."""
    n = len(A)
    I = np.eye(n)
    A_inv = []

    for i in range(n):
        b = I[:, i]
        x = solve_system_lu_decomposition(A, b)
        A_inv.append(x)

    return np.array(A_inv).T

