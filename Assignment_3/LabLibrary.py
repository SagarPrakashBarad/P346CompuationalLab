global eps
eps = 1/10**6

class Matrix:
    """Matrix Class
       Functions: Make a 2D array (Matrix) from a sting of numbers or make them from
       a single element
       INPUT:
       1) shape = (2,2), a = string (1 2 2 1)
       2) shape = (2,2), a = 3
       OUTPUT:
       1) |1 2|
          |2 1|
       2) |3 3|
          |3 3|
    """

    def __init__(self, dims, fill): # initialize the paramters of matrix like length
                                    # of rows and columns
        self.rows = dims[1]
        self.cols = dims[0]

        if isinstance(fill,list): # checks for string or not
            k = 0
            lst = []
            for i in range(dims[0]): # if string fills the array numbers of string
                col = []             # separated by space
                for j in range(dims[1]):
                    col.append(fill[k])
                    k += 1
                #col.append(x[i])
                lst.append(col)
            self.ele = lst
        else: # make every element the same number
            self.ele = [[fill] * self.cols for i in range(self.rows)]

        """2D array of class is accessed by using .ele attribute say mat.ele[l][k]"""

    def __str__(self): # print the 2D array in matrix form
        m = len(self.ele)
        mtxStr = ''

        mtxStr += '--------------------- output ---------------------\n'
        for i in range(m):
            mtxStr += ('|' + ', '.join( map(lambda x: '{0:8.3f}'.format(x), self.ele[i])) + '| \n')
        mtxStr += '--------------------------------------------------'

        return mtxStr

def forward_substitution(mat,b):
    """solves y against a lower triangular matrix"""
    y = []
    n = len(mat.ele)
    for i in range(n): # main loop
        sum = 0
        if mat.ele[i][i] == 0: continue # skips singular pivot entries to avoid division by zero
        for j in range(i): # for all off daigonal elements
            sum += mat.ele[i][j]*y[j] # take sum of product non-pivot elements with corresponding y value
        y.append((b[i] - sum)/mat.ele[i][i]) # append the value to new list

    return y # return the solution as list

def backward_substitution(mat,y):
    x = []
    n = len(mat.ele)
    for i in range(n):
        x.append(0)

    x[n-1] = y[n-1]/mat.ele[n-1][n-1]
    for i in range(n-1,-1,-1): # loop interates in reverse (main loop)
        sum = 0
        if mat.ele[i][i] == 0: continue # skips singular pivot entries to avoid division by zero
        for j in range(i+1,n):
            sum += mat.ele[i][j]*x[j] # take sum of product non-pivot elements with corresponding y value
        x[i] = (1/mat.ele[i][i])*(y[i] - sum)
        
    return x # return the solution as list

def is_symmetric(mat): # checks if matrix is is_symmetric or not?
    n = len(mat.ele)
    for i in range(n):
        for j in range(n):
            if mat.ele[i][j] != mat.ele[j][i]: # it checks that all (i,j) and (j,i) are equal, if not then they are not symmetric
                return False
    return True

def LUdecom(mat):
    (n, m) = (len(mat.ele), len(mat.ele[0]))
    L = Matrix(dims=(n,m),fill=0)
    U = Matrix(dims=(n,m),fill=0)

    for i in range(0,n):
            # partial pivoting
        maxelement = i
        if abs(mat.ele[i][i]) <= eps:
            for j in range(i+1,n):
                if abs(mat.ele[j][i]) > abs(mat.ele[maxelement][i]):
                    maxelement = j
            (mat.ele[i], mat.ele[maxelement]) = (mat.ele[maxelement], mat.ele[i])

    for i in range(0,n): # main loop
        for k in range(i,n):

            sum = 0
            for j in range(i): # secondary loop
                sum += (L.ele[i][j]*U.ele[j][k])
            U.ele[i][k] = mat.ele[i][k] - sum

        for k in range(i,n):
            if i == k:
                L.ele[i][i] = 1
            elif (U.ele[i][i] != 0):
                sum = 0
                for j in range(i):
                    sum += (U.ele[j][i]* L.ele[k][j])
                L.ele[k][i] = (mat.ele[k][i] - sum)/U.ele[i][i]

    return L,U


def solveLU(mat,b): # takes a matrix and b vector
    L,U = LUdecom(mat) # renders an L and U matrix from original matrix
    y = forward_substitution(L,b) # y the backward_substitution solution of L and b
    x = backward_substitution(U,y) # x the forward_substitution solution of U and y

    return x

def solvecholemat(mat,b): # works same as solveLU
    L,LT = cholesky(mat)
    y = forward_substitution(L,b)
    x = backward_substitution(LT,y)

    return x


def gaussjordan(mat,b):
    (n, m) = (len(mat.ele), len(mat.ele[0]))
    """gaussjordan takes a matrix (2D array) as input and reduce it to Reduced Row Echelon Form.
        Returns True if successful, False if 'm' is singular."""

    # Main loop
    for i in range(0,n):
        # partial pivoting
        maxelement = i
        if abs(mat.ele[i][i]) <= eps:
            for j in range(i+1,n):
                if abs(mat.ele[j][i]) > abs(mat.ele[maxelement][i]):
                    maxelement = j
            (mat.ele[i], mat.ele[maxelement]) = (mat.ele[maxelement], mat.ele[i])
            (b[i], b[maxelement]) = (b[maxelement], b[i])

        # Division of pivot row
        pivot = mat.ele[i][i]
        for k in range(i,n):
            mat.ele[i][k] /= pivot
        b[i] /= pivot

        # Elimination loop
        for j in range(n):
            if i == j or mat.ele[j][i] == 0: continue

            ratio = mat.ele[j][i]
            for k in range(i,n):
                mat.ele[j][k] -= ratio*mat.ele[i][k]
            b[j] -= ratio*b[i]

    return b, mat


def cholesky(mat):
    if is_symmetric(mat):
        (n, m) = (len(mat.ele), len(mat.ele[0]))
        L1 = Matrix(dims=(n,m),fill=0) # makes two null matrix
        L2 = Matrix(dims=(n,m),fill=0)

        for i in range(0,n):
            for j in range(0,i+1): # filling the matrix to get L1
                sum = 0
                for k in range(0,j):
                    sum += L1.ele[i][k]*L1.ele[j][k]
                if (i == j):
                    L1.ele[j][j] = (mat.ele[j][j] - sum)**(0.5)
                else:
                    L1.ele[i][j] = (1/L1.ele[j][j])*(mat.ele[i][j] - sum)
                L2.ele[j][i] = L1.ele[i][j] # taking transpose of L1 to get L2

        return L1, L2




def GaussSiedel(A,b,g):
     x = g.copy()
     t = 1
     iter = 0

     while t > eps: # termination when tolerance becomes less than eps
         for i in range(len(A.ele)):
             l = 0
             u = 0
             for k in range(len(A.ele)):
                        if k < i:
                            l += A.ele[i][k]*x[k]
                        elif k > i:
                            u += A.ele[i][k]*g[k]
                            
             x[i] = (1/A.ele[i][i])*(b[i] - l - u)
         iter += 1
         t = max([abs(x[i]-g[i]) for i in range(0,len(A.ele))]) # find difference of all elements of all elemnts of vectors
         g = x.copy()
         print("Iteration: {} | Current Solution: {}".format(iter, x))

     return x,iter

def DDcheck(A,i): # checks whether a pivot is daigonally dominant
        sum = 0
        for j in range(len(A.ele[0])):
            if i != j:
                sum += A.ele[i][j]
            sum = abs(sum)
        if abs(A.ele[i][i]) < sum:
             return 1

def makeDD(mat,b):
        (n, m) = (len(mat.ele), len(mat.ele[0]))

        for i in range(0,n):
         maxelement = i
         if DDcheck(mat,i): # if pivot is not daigonally dominant
                for j in range(i+1,n):
                    if abs(mat.ele[j][i]) > abs(mat.ele[maxelement][i]):
                     maxelement = j
                (mat.ele[i], mat.ele[maxelement]) = (mat.ele[maxelement], mat.ele[i]) # swap places with row having largest column element
                (b[i], b[maxelement]) = (b[maxelement], b[i]) # swaping to simualte an augemntated matrix

def jacobi(mat,b,g):
     x = g.copy()
     t = 1
     iter = 0

     while t > eps:
         iter = iter + 1
         for i in range(0,len(mat.ele)):
             sum = 0
             for j in range(0,len(mat.ele)):
                 if i != j:
                     sum += mat.ele[i][j]*g[j]
             x[i] = (1/mat.ele[i][i])*(b[i]- sum)
             t = abs((x[i] - x[i])*100/x[i])
             g[i] = x[i]
            
         print("Iteration: {} | Current Solution: {}".format(iter, x))      
    
     return x,iter
