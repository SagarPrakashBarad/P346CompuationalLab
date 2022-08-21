""" A matrix class capable of generating a matrix (2d array) from a sequence of from input.txt"""
class Matrix:

    def __init__(self, dims, fill):

        self.rows = dims[0]
        self.cols = dims[1]

        if isinstance(fill,list):
            k = 0
            lst = []
            for i in range(dims[0]):
                col = []
                for j in range(dims[1]):
                    col.append(fill[k])
                    k += 1
                lst.append(col)
            self.ele = lst
        else:
            self.ele = [[fill] * self.cols for i in range(self.rows)]

    def __str__(self): # printing a matrix
        m = len(self.ele)
        mtxStr = ''

        mtxStr += '--------------------- output ---------------------\n'
        for i in range(m):
            mtxStr += ('|' + ', '.join( map(lambda x: '{0:8.3f}'.format(x), self.ele[i])) + '| \n')
        mtxStr += '--------------------------------------------------'

        return mtxStr

def matmul(A, B): #matrix-matrix multiplication

		if isinstance(B, Matrix):
			C = Matrix( dims = (B.cols, A.rows), fill = 0)

			#Multiply the elements in the same row of the first matrix
			#to the elements in the same col of the second matrix
			for i in range(A.rows):
				for j in range(B.cols):
					acc = 0

					for k in range(A.rows):
						acc += A.ele[i][k] * B.ele[k][j]

					C.ele[j][j] = acc

		return C

def dotprod(C, D): # dot product

        dp = 0
        for i in range(C.cols):
            dp += C.ele[i][0] * D.ele[i][0]

        return dp

"""complex class for creating and complex object"""
class mycomplex:

    def __init__(self,x ,y):
        self.real = x
        self.imag = y
        self.mod = (self.real**2 + self.imag**2)**(1/2)

    def __str__(self):
        cnstr = ('{} + {}i'.format(self.real, self.imag))
        return cnstr

"""Addition and subtraction of complex numbers functions"""
def add(A,B):

    C = mycomplex(x=0,y=0)
    C.real = A.real + B.real
    C.imag = A.imag + B.imag
    C.mod = (C.real**2 + C.imag**2)**(1/2)

    return C

def mul(A,B):
    D = mycomplex(x=0,y=0)
    D.real = A.real*B.real - A.imag*B.imag
    D.imag = A.imag*B.real + B.imag*A.real
    D.mod = (D.real**2 + D.imag**2)**(1/2)

    return D

"""Python implementation of a factorial function and a sum of n odd numbers function"""
def factorial(n):

    res = 1
    for i in range(2, n+1):
        res *= i
    return res

def oddsum(n):

    sum = 0
    for i in range(1, n):
        sum += (2*i + 1)
    return sum
"""Python implementation of AP, GP and HP series for common difference 1.5 and common ratio 0.5."""
def progressionsum(n,r,d):

    sum = []
    a = 1.0
    for i in range(0,3):
        sum.append(a)

    for i in range(2,n):
        sum[0] += a + (i-1)*d
        sum[1] += a*(r)**(i-1)
        sum[2] += 1/(a + (i-1)*d)

    return sum
"""Python implementation of AP, GP and HP series for common difference 1.5 and common ratio 0.5."""
def summnationseries(n):
    sum = []
    a = 1
    for i in range(1,n+1):
        a += ((-1)**(i+1))/2**i
        b = float("{0:.4f}".format(a))
        sum.append(b)

    return sum
