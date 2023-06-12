global max_iter, eps
max_iter = 10
eps = 10**(-3)

class Eigenvalue:

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


        def eigenvalue(self,x0):

            # matmul subroutine
            self.V1 = []
            self.V2 = []
            for i in range(self.rows):
                self.V1.append(0)
                self.V2.append(0)
            self.ev1 = 1
            ev0 = 0
            k = 1
            rn = 0
            rd = 0
            while abs(self.ev1 - ev0) > eps or k < max_iter:
                ev0 = self.ev1
                for i in range(k):
                    for l in range(self.rows):
                        for m in range(self.cols):
                            self.V1[l] += self.M[l][m]*x0[m]

                    x0 = self.V1.copy()

                for i in range(self.rows):
                    for l in range(self.cols):
                        self.V2[i] += self.M[i][l]*self.V1[l]
                    rn += self.V1[i]*self.V2[i]
                    rd += self.V1[i]*self.V1[i]

                self.ev1 = rn/rd
                k += 1

            return k

        def eigenvectors(self):
            self.AM = self.M.copy()
            self.ev = []
            mod = 0
            for i in range(self.rows):
                self.AM[i][i] -= self.ev1

            self.ev.append(-1*(self.AM[1][1]*self.AM[0][2] - self.AM[1][2]*self.AM[0][1]))
            self.ev.append(self.AM[1][0]*self.AM[0][2] - self.AM[1][2]*self.AM[0][0])
            self.ev.append(self.AM[1][1]*self.AM[0][0] - self.AM[1][0]*self.AM[0][1])

            for i in range(self.rows):
                mod += self.ev[i]**2

            mod = (mod)**(1/2)

            for i in range(self.cols):
                self.ev[i] /= mod
