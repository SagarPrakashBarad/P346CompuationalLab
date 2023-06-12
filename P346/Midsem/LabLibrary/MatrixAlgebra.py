global eps
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
                sum += mat.ele[i][j]*y[j] # take sum of product non-pivot elements with corresponding y value
            self.y.append((B[i] - sum)/self.M[i][i])

    def backward_substitution(self,y):

        self.x = y.copy()

        self.x[self.rows -1] = y[self.rows -1]/self.M[self.rows -1][self.cols -1]
        for i in range(self.rows - 1,-1,-1):
            sum = 0
            if self.M[i][i] == 0: continue # skips singular pivot entries to avoid division by zero
            for j in range(i+1,self.rows):
                sum += self.M[i][j]*self.x[j] # take sum of product non-pivot elements with corresponding y value
            self.x[i] = (1/mat.ele[i][i])*(y[i] - sum)

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
                        sum += self.L[K][j]*self.U[j][i]

                    self.L[k][i] = (self.M[k][i] - sum)/self.M[i][i]

    def Cholesky(self):

        (n,m) = (self.rows, self.cols)
        self.L = [[0] * m for i in range(n)]
        self.U = [[0] * m for i in range(n)]
