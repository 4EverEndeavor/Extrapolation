from math import factorial 
from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np
from numpy import sin,cos
import decimal
import time
import random

# this is to calculate runtime
startTime = time.clock()

# set context for decimal module
decimal.Context(prec = 50, rounding = 'ROUND_DOWN')
decimal.getcontext().prec = 50

# several different functions to define data points
def f(x):
    # This adds a Gaussian distribution to the data points 
    noise = random.gauss(1,0)
    # test against 4 differenet kinds of functions
    #f = sin(x)
    f = cos(x)
    #f = 3*sin(x) + np.log(x) + np.sqrt(x)
    #f = np.exp(x)
    #f = 3*x
    #f = np.log(x)
    # with noise
    return decimal.Decimal(noise * f)
    # or without
    #return decimal.Decimal(f)
    #return f

# data points
# distance between each point
d = 1
# number of points
n = 30
# x and y arrays to hold data points
y = []
x = []
# there will be n data points
for i in range(1,n+1):
    # each x[i] = i*d, f(x[i]) = f(i*d)
    y.append(f(i*d))
    #x.append(i)
    x.append(d*i)

'''This is to compare the Talor series with the Lagrange polynomial.
The assumption is that a Taylor series centered around the middle data
point is roughly equal to the Lagrange polynomial. For a Taylor series
of cos not centered at zero, the derivatives are required for each
term. Instead of calculating derivatives for each term, two summations
were created incorperating sin and cos. '''
def cosExpansion(x,a):
    total = 0
    terms = n
    # cos terms, or even terms
    for k in range(int(terms/2)):
        numerator = (-1)**k * cos(a) * (x-a)**(2*k)
        denominator = factorial(2*k)
        total += (numerator / denominator)
    # sin terms, or odd terms
    for k in range(1,int(terms/2)+1):
        numerator = (-1)**k * sin(a) * (x-a)**(2*k-1)
        denominator = factorial(2*k-1)
        total += (numerator / denominator)
    return total


# This is used in several functions. It comes up frequently 
# because of how the matrix system was solved. 
def nChoosek(n,k):
    # factorials of negative numbers do not exist
    if n < k:
        # for our purposes, it is convenient to define them as zero
        return 0
    else:
        # binomial coefficient formula 
        r = factorial(n)/(factorial(n-k)*factorial(k))
        return decimal.Decimal(r)

# Elements of the reduced upper triangular matrix
def Mij(i,j):
    summation = decimal.Decimal('0')
    for k in range(0,i):
        r = decimal.Decimal(pow((k+1)*d,j-1)*pow(-1,i-k-1))
        #summation += nChoosek(i-1,k)*pow((k+1)*d,j-1)*pow(-1,i-k-1)
        summation += nChoosek(i-1,k)*r
    return int(summation)

# Creates the upper triangular matrix
def upperMatrix(n,d):
    row = range(1,n+1)
    col = range(1,n+1)
    tempRow = []
    tempArr = []
    for i in col:
        for j in row:
            tempRow.append(Mij(i,j))
        tempArr.append(tempRow)
        tempRow = []
    return np.matrix(tempArr)

# Elements of the y column after row reduction is performed
'''def Bi(i):
    summation = 0
    for k in range(0,i):
        summation += nChoosek(i-1,k)*y[k]*pow(-1,i-k-1)   
    return summation'''

# Create a matrix system, so that each matrix element is 
# only computed once. 
M = upperMatrix(n,d)
'''B = []
for i in range(1,n+1):
    B.append(Bi(i))'''

'''This creates a matrix containing Pascal's triangle with negative 
diagonals. It appears in the formula because of the binomial 
expansion coefficients created from row reduction.  '''
'''def pascal(numberOfPoints):
    # create array of zeros
    P = np.zeros([numberOfPoints,numberOfPoints])
    # first row is all ones 
    P[0][:] = 1
    # diagonal is all ones
    for i in range(numberOfPoints):
        P[i][i] = 1
    # This makes the triangle 
    for row in range(numberOfPoints-1):
        for col in range(numberOfPoints-1):
            # by taking the two elements above it
            r = row + 1
            c = col + 1
            # and adding them together. 
            P[r][c] = P[r][c-1] + P[r-1][c-1]
    # Now make every other diagonal negative
    for row in range(numberOfPoints):
        for col in range(numberOfPoints):
            r = row
            c = col
            P[r][c] *= (-1)**(col+row)
    # return a matrix
    P = np.matrix(P)
    return P'''

'''This creates a matrix containing Pascal's triangle with negative 
diagonals. It appears in the formula because of the binomial 
expansion coefficients created from row reduction.  '''
def pascal(numberOfPoints):
    # create array of zeros
    P = np.zeros([numberOfPoints,numberOfPoints],dtype = object)
    # first row is all ones 
    P[0][:] = decimal.Decimal(1)
    # diagonal is all ones
    for i in range(numberOfPoints):
        P[i][i] = decimal.Decimal(1)
    # This makes the triangle 
    for row in range(numberOfPoints-1):
        for col in range(numberOfPoints-1):
            # by taking the two elements above it
            r = row + 1
            c = col + 1
            # and adding them together. 
            P[r][c] = P[r][c-1] + P[r-1][c-1]
    # Now make every other diagonal negative
    for row in range(numberOfPoints):
        for col in range(numberOfPoints):
            r = row
            c = col
            #P[r][c] = decimal.Decimal(P[r][c]* (-1)**(col+row))
            P[r][c] *= (-1)**(col+row)
    # return a matrix
    P = np.matrix(P)
    return P

'''This function returns an array of all the possible 
permutations of decreasing, non-equal, numbers in a given 
range. This list will be used to select certain matrix elements 
later on. The length of the decreasing sequence is given by 
the number of matrix indexes that need to be swapped.''' 
def swapIndexes(numberOfSwaps, start, stop):
    # This function should return a predetermined amount of 
    # permutations. 
    possiblePermutations = nChoosek(stop-start-1,numberOfSwaps)
    
    # create a list of indexes
    indexes = []
    for i in range(numberOfSwaps):
        indexes.append(stop - i - 1)
    
    # current index = last  
    currentIndex = (numberOfSwaps-1)
    
    # this is what the function will return 
    allPermutations = []
    # add the initial permutation to the list
    '''Deep copy was necessary to make a brand new list, that was 
    not effected by operations on the previous list, and vice versa.'''
    i = deepcopy(indexes)
    allPermutations.append(i)
    # create all the correct subscript permutations
    while possiblePermutations > 0:
        # if at the last element
        if currentIndex == numberOfSwaps - 1:
            # and it's greater than two
            if indexes[currentIndex] > start+1:
                # decrement the value
                indexes[currentIndex] -= 1
                # and tally off one permutation 
                possiblePermutations -= 1
                # and add it to the list
                i = deepcopy(indexes)
                allPermutations.append(i)
            # otherwise, if it is two already  
            elif indexes[currentIndex] == start+1:
                # and there's only one swap
                if currentIndex == 0:
                    # we are done!
                    break
                # otherwise
                else:
                    # move to previous number 
                    currentIndex -= 1
            # just in case something squirly happens...
            else:
                print('value: ', indexes[currentIndex])
                print('current: ', currentIndex)
                print('start: ', start)
                print('swapIndexes ERROR')
                print(indexes[currentIndex])
                break
        # if back to the 1st element
        elif currentIndex == 0:
            # and it's greater than two
            if indexes[currentIndex] > indexes[currentIndex+1]+1:
                # decrement the value
                indexes[currentIndex] -= 1
                # and tally off one permutation
                possiblePermutations -= 1
                # make all the trailing numbers 
                # n, n-1, n-2, ...
                for i in range(1,numberOfSwaps):
                    # each element is one less than previous 
                    indexes[i] = indexes[i-1] - 1
                # and return to last element to start over
                currentIndex = numberOfSwaps - 1 
                # add the permutation to the list
                i = deepcopy(indexes)
                allPermutations.append(i)
            # and it is at its lowest value
            else:
                # we're done! 
                possiblePermutations = 0
        # or if somewhere in the middle, 
        else:
            # and greater than the next element + 1
            if indexes[currentIndex] > indexes[currentIndex+1]+1:
                # decrement the value
                indexes[currentIndex] -= 1
                # and tally off a permutation
                possiblePermutations -= 1
                # fix trailing numbers
                for i in range(currentIndex + 1,numberOfSwaps):
                    # each element is one less than previous 
                    indexes[i] = indexes[i-1] - 1
                # and return to last element
                currentIndex = numberOfSwaps - 1
                # add the permutation to the list
                i = deepcopy(indexes)
                allPermutations.append(i)
            # and value cannot be decremented
            else:
                # move to previous element
                currentIndex -= 1
    # return the entire list of subscripts
    return allPermutations
    
''' Now the array of indexes from swapIndexes is used to select 
matrix elements from Mij. Each element in the array represents a 
sequence of swaps performed on matrix elements. After the swaps 
are performed, the sequence represents a set of matrix elements 
that are to be multiplied together. This next function, returns 
the sum of all multiplied matrix elements for each given number
of swaps. '''
def permutationSwapper(numberOfSwaps, start, stop):
    # this situation does come up, and will always equal 1
    if start == stop:
        return 1
    # create the matrix subscripts
    matrixSubscripts = []
    for j in range(start,stop):
        matrixSubscripts.append([j,j+1])
    # in each term, we'll need one product with zero swaps
    if numberOfSwaps == 0:
        # define a product
        product = 1
        # duplicate list of subscripts
        copy = deepcopy(matrixSubscripts)
        for element in range(len(copy)):
            product *= M.item(copy[element][0]-1,copy[element][1]-1)
        return product
    # get the swap positions 
    allPermutations = swapIndexes(numberOfSwaps, start, stop)
    # each product of matrix elements will be added to this
    summation = 0
    # use the permutations to perform swaps on the subscripts
    for eachPermutationSet in range(len(allPermutations)):
        # create a copy of the subscripts
        copy = deepcopy(matrixSubscripts)
        for eachSwap in range(numberOfSwaps):
            # designate the postion of the swap
            pos = allPermutations[eachPermutationSet][eachSwap]-start
            # perform the swap
            copy[pos-1][1],copy[pos][1] = copy[pos][1],copy[pos-1][1]
        # reset product to 1
        product = 1
        # Take all the subscripts
        for element in range(len(copy)):
            # and map them to matrix elements,
            # multiply them together
            product *= M.item(copy[element][0]-1,copy[element][1]-1)
        # and add the product to the sum
        summation += product 
    # finally, return the sum
    return summation


'''Now we sum together all the possible swaps and divide. The 
following function was listed in the write up, and Power Point as
Rij, or ratios. These are the stored elements of the matrix file. '''
def ratio(start, stop):
    # the first term is always one
    if start == stop:
        #return (1/M.item(start-1,start-1))
        element = M.item(start-1,start-1)
        return 1/decimal.Decimal(element)
    # this will add alternating terms in the numerator
    numerator = 0
    for swaps in range(stop-start):
        numerator += ((-1)**(swaps+stop-start)) * permutationSwapper(swaps,start,stop)
        
    # create a denominator
    denominator = 1
    for element in range(start-1, stop):
        denominator *= M.item(element,element)
    
    #return numerator
    '''As the exponents get larger, or statistical error is 
    introduced, higher accuracy is necessary to create a suitable fit.
    In an effort to accomplish that, the decimal module was called. 
    Although the decimals can now be calculated to an arbitrary 
    accuracy, the matrix file will only hold numbers of a maximum size.
    more work will need to be done to alleviate this. '''
    return decimal.Decimal(numerator)/decimal.Decimal(denominator)
    #return (numerator / denominator)
    #return denominator


# This function creates a matrix to store all ratios for future use.
'''def createRatioMatrix(size):
    # create an n x n array of zeros 
    matrix = np.zeros((size,size))
    # page through each row..
    for r in range(size):
        # ... and column
        for c in range(size):
            # to create and upper triangular matrix
            if r <= c:
                # poplulate elements with the ratio function
                matrix[r][c] = ratio(r+1,c+1)
                print(matrix[r][c],ratio(r+1,c+1))
    # retrun the matrix
    return matrix'''

# This function creates a matrix to store all ratios for future use.
def createRatioMatrix(size):
    # create an n x n array of zeros 
    matrix = []
    # page through each row..
    for r in range(size):
        row = []
        # ... and column
        for c in range(size):
            # to create and upper triangular matrix
            if r <= c:
                # poplulate elements with the ratio function
                row.append(ratio(r+1,c+1))
            else:
                row.append(0)
        matrix.append(row)
    # retrun the matrix
    return np.array(matrix)

''' This creates the term matrix, and only needs doing once.
Creating a ratio matrix of 25 x 25 took a runtime of 
5 hours, 38 minutes and 17.7 seconds. '''
# new xps macine created 30 x 30 24 hours, 23 minutes
##################   ratioMatrix = createRatioMatrix(n)

''' Some of these calculations are very lengthy, so 
it pays to keep track of how long they take. After 20 data 
points, the run times are longer than one hour, and rapidly 
go up from there. '''
def runtime():
    # startTime was initiated at the beginning of the program
    t = time.clock() - startTime
    hours = t // (60*60)
    minutes = t // 60
    seconds = t % 60
    # dispaly runtime
    print('\truntime: ', minutes, 'minutes', seconds,' seconds')


'''The following two functions read and write matrixes to 
a file. Calculating the matrix elements requires hefty runtime, 
and they need only be calculated once, as these figures do not
depend on the data points at all. Once calculated, the matrixes
can be called very quickly from a text file. '''
def writeMatrix(matrix, matrixFile):
    # open the file 
    f = open(matrixFile, 'w')
    # page through rows
    for row in range(len(matrix)):
        # create an empty string to add onto
        line = ''
        # and columns 
        for col in range(len(matrix[row])):
            # add each matrix element to the string
            line += str(matrix[row][col])
            # separate by a comma, so the matrix can be read easily
            line += ','
        # write the string to the file
        f.write(line)
        # separate the rows by newline escapes sequences 
        f.write('\n')
    # close the file. 
    f.close()
    
# Write the matrix to a file. This only needs to be done once
############### writeMatrix(ratioMatrix, 'RatioMatrix.txt')
############### writeMatrix(ratioMatrix, 'xpsMatrix.txt')
   
# Retrieve the stored matrix from the file 
def readMatrix(rows, matrixFile):
    # open the file
    f = open(matrixFile, 'r')
    # create an emtpy list
    matrix = []
    # page through lines in the file to create 
    for row in range(rows):
    #for row in f:
        # rows in the matrix
        r = f.readline()
        # turn the file string into a list by splitting at the ,'s
        r = r.split(',')
        # remove the '\n' from the end of the list
        r.pop()
        # convert each element in the list to float
        for i in range(len(r)):
            #r[i] = float(r[i])
            r[i] = decimal.Decimal(r[i])
        # add each row list to the matrix list
        matrix.append(r)
    # turn the array into a numpy matrix
    matrix = np.matrix(matrix)
    # close the file 
    f.close()
    # and return the stored matrix 
    return matrix

# Computes each Lagrange coefficient my multiplying truncated matrixes
def pCoefficient(i,p,matrix):
    # for indexing reasons,
    i -= 1
    # row vector of y value
    yVector = np.matrix(y)[:,0:len(matrix)]
    ''' This part is the essence of the algorithm. 
    Each coefficient can be computed by simply multiplying 
    a predifined, one-size-fits-all matrix by the y values.'''
    matrix = np.matrix(matrix)
    # formula can be found in write up as well
    r = yVector * p[:,i:n] * matrix.getT()[i:,i]
    return r.item(0,0)
    #return float(yVector * p[:,i:n] * matrix.getT()[i:,i])
    
    
# store each coefficient in an array
def createCoefficientArray(numberOfPoints, ratioMatrix):
    # negate every other row, so that the 1st term is always positive
    for i in range(numberOfPoints):
        ratioMatrix[i] *= (-1)**(2*numberOfPoints)
    # create a Pascal matrix for computing coefficients
    p = pascal(numberOfPoints)
    # list for storing calculated coefficients
    coefficients = []
    # append each coefficient to the list
    for i in range(1,numberOfPoints+1):
        coefficients.append(pCoefficient(i,p,ratioMatrix))
    # return the list
    return coefficients


''' This is a cascading algorithm for computing p values.
It multiplies, then adds, to prevent the exponential terms
from getting too large. The end result of any fit is going to 
be on the same order as the function value. However, higher order
terms in series solutions tend to grow very large, and are prone
to computational error. 

p = a1 + a2*x + a3*x^2 + a4*x^3 + a5*x^4

We already know that these huge terms will cancel eachother out 
to a great degree, so splitting up the exponent can keep the numbers 
within a managable range. 

  = a1 + {a2 + [a3 + (a4 + a5*x)*x]*x}*x '''
def cascade(x, coefficients):
    # start with last coefficient
    m = len(coefficients)
    p = coefficients[m-1]
    # loop m-1 times
    for i in range(m-1):
        # multiply
        p *= decimal.Decimal(x)
        # then add
        p += coefficients[m-i-2]
        # and repeat
    # return p value
    return p

'''
# Power rule derivative of power series
def firstDerivative(x, coefficients):
    p = 0
    # the first coefficient drops off
    for i in range(1, n):
        p += i*coefficients[i] * x**(i-1)
        p = float(p)
    return p'''

'''Because this is a power series, any derivative can easily 
be computed by applying a simple power rule. '''
def derivative(x, coefficients,j):
    # initial power series 
    p = 0
    # convert x to decimal
    x = decimal.Decimal(x)
    for i in range(j,n):
        # generalized power rule 
        r = factorial(i)/factorial(i-j)
        # conver r to decimal
        r = decimal.Decimal(r)
        # add to the power series 
        p += r * coefficients[i] * x**(i-j)
    # If derivative number exceeds polynomial order, 
    if j >= n:
        # all terms drop off.
        return 0
    return p

# Once again, apply power rule, but this time for 1st integral
def integral(x, coefficients):
    # initial power series
    p = decimal.Decimal(0)
    # convert x to decimal
    x = decimal.Decimal(x)
    # the first coefficient drops off
    for i in range(n):
        # copy i so it may be converted to decimal
        c = decimal.Decimal(1/(i+1))
        e = decimal.Decimal(i)
        # power series summation
        p += c * coefficients[i] * x**(e+1)
        #p = float(p)
    return p

# Creates, fits, and plots all functions. 
def display(numberOfPoints):
    # midpoint
    m = (n*d)/2
    # arrays for the Lagrange, Taylor, integral and derivatives
    pfit,tfit,intP,dpdx,dpdx4,dpdx8,dpdx12,dpdx16 = [],[],[],[],[],[],[],[]
    # read in the matrix
    # ratios = readMatrix(numberOfPoints,'RatioMatrix.txt')[:numberOfPoints,:numberOfPoints]
    ratios = readMatrix(numberOfPoints,'xpsMatrix.txt')[:numberOfPoints,:numberOfPoints]
    # now create the coefficients
    coefficients = createCoefficientArray(numberOfPoints,ratios)
    # this is the x axis of the fit plot
    xRange = np.arange(-3,n*d+3,0.01)
    #xRange = np.arange(-1,n+20,0.01)
    # populate y values for pfit
    for i in range(len(xRange)):
        #pfit.append(cascade(xRange[i],coefficients))
        pfit.append(cascade((1/d)*xRange[i],coefficients))
     
    # populate y values for derivatives 
    '''for i in range(len(xRange)):
        #pfit.append(cascade(xRange[i],coefficients))
        #dpdx.append(derivative((1/d)*xRange[i],coefficients,1))
        dpdx.append(derivative((1/d)*xRange[i],coefficients,1))'''
    '''for i in range(len(xRange)):
        #pfit.append(cascade(xRange[i],coefficients))
        dpdx8.append(derivative((1/d)*xRange[i],coefficients,8))
    for i in range(len(xRange)):
        #pfit.append(cascade(xRange[i],coefficients))
        dpdx12.append(derivative((1/d)*xRange[i],coefficients,12))
    for i in range(len(xRange)):
        #pfit.append(cascade(xRange[i],coefficients))
        dpdx16.append(derivative((1/d)*xRange[i],coefficients,16))'''
       
    # populate y values for integral
    '''for i in range(len(xRange)):
        #pfit.append(cascade(xRange[i],coefficients))
        adjustment = integral(m, coefficients)
        intP.append(integral((1/d)*xRange[i],coefficients)-adjustment)'''
        
    # populate y values for tfit
    '''for i in range(len(xRange)):
        tfit.append(cosExpansion(xRange[i],(n*d)/2))'''
        
    # plot of fit
    plt.plot(xRange,pfit,'b')
    # taylor series
    #plt.plot(xRange,tfit,'m--')
    # derivatives
    #plt.plot(xRange,dpdx,'g--')
    '''plt.plot(xRange,dpdx8,'--')
    plt.plot(xRange,dpdx12,'-.')
    plt.plot(xRange,dpdx16, ':')'''
    # integral 
    #plt.plot(xRange,intP,'c:')
    
    # Plot the points
    plt.plot(x,y,'r*')
    # add pretty arrows to explain stuff 
    plt.title('int(p), p prime, Taylor Series, and Points')
    #plt.annotate('p(x)', xy=(2*np.pi, 1), xytext=(2*np.pi, 2), arrowprops=dict(facecolor='black', shrink=0.05,width = 0.1, headwidth = 4))
    #plt.annotate('p\'(x)', xy=(np.pi/2, -1), xytext=(np.pi/2, -2), arrowprops=dict(facecolor='black', shrink=0.05,width = 0.1, headwidth = 4))
    #plt.annotate('Taylor', xy=(18, 2), xytext=(13, 2), arrowprops=dict(facecolor='black', shrink=0.05,width = 0.1, headwidth = 4))
    #plt.annotate('int(p)', xy=(7*np.pi/2, -1), xytext=(7*np.pi/2, -2), arrowprops=dict(facecolor='black', shrink=0.05,width = 0.1, headwidth = 4))
    plt.xlim(min(x)-1,max(x)+2)
    
    # plot y limits
    top = float(max(y))
    bottom = float(min(y))
    plt.ylim(bottom - 1.5,top + 1.5)

    plt.show()


'''This function will input a ratio matrix, and calculate
the difference between |p(x[i])-y[i]| for each known data point.'''
def minimizationFunction(matrix,testPoint,numberOfPoints):
    # make expanded ratio matrix, with new test point
    for component in range(1,len(testPoint)+1):
        matrix[component][-1] = testPoint[component-1]
    # compute coefficients
    coefficients = createCoefficientArray(numberOfPoints, matrix)
    # loop through each data point
    error = 0
    for dataPoint in range(numberOfPoints):
        # calculate difference |p(x[i])-y[i]|
        a = y[dataPoint]
        b = cascade(d*(dataPoint+1), coefficients)
        error += abs(a - b)
    return error
    
'''This function is a proof of concept, that a monte carlo algorithm
can in fact be used to expand the ratio matrix by another column. 
The runtime however, is actually much longer than calculating the 
ratios directly. Perhaps, a new, smarter optimization algorithm could
be employed here. Ultimately, the real problem here, is that each 
point adds another dimension to the optimization. That is, an n 
dimensional matrix will require a n-2 dimensional optimization, because
the top row and diagonal are already known. '''
def monteCarlo(existing):
    print('This might take a while. ')
    # expand the matrix by one column and row
    newSize = existing + 1
    
    # create a copy of the term matrix, but with a new row and column
    monteCarloMatrix = np.zeros((newSize,newSize))
    monteCarloMatrix = np.array(monteCarloMatrix)
    
    # read in the partially computed ratio matrix
    partialMatrix = readMatrix(existing,'RatioMatrix.txt')[:existing,:existing]
    # go through each element of the ratio matrix and copy it
    for row in range(existing):
        for col in range(existing):
            monteCarloMatrix[row][col] = partialMatrix.item(row,col)
     
    # the top row is always +1 or -1
    monteCarloMatrix[0][-1] = (-1)**(2*newSize)
    
    # the bottom right element will always be 1/M(n+1,n+1)
    monteCarloMatrix[-1][-1] = 1/Mij(newSize,newSize)
    
    # create an n-tuple point to center the initial search on 
    center = []
    for row in range(1,newSize-1):
        center.append(monteCarloMatrix[row][-2])
        
    #center = [-2.08,1.458,-0.4166]
    # initial search radius
    radius = 100
    # to make sure the loop does not execute forever
    counter = 0
    # set an initial tolerance
    error = 100
    # how accurate the result will be
    tolerance = 0.1
    # each time the loop will test a predefined number of points
    numberOfPoints = 20
    # this loop will execute until the point converges 
    while(abs(error) > tolerance):
        # create an array to store the test points
        testPointArray = []
        # make numberOfPoints random points
        for points in range(numberOfPoints):
            # each point will be a tuple of size n
            point = []
            for components in range(len(center)):
                # choose a random value within the radius of the center
                component = random.uniform(center[components]-radius,center[components]+radius)
                # add the value to the point
                point.append(component)
            # add the point to the array
            testPointArray.append(point)
            
        # test each point 
        for point in range(numberOfPoints):
            #print(error, radius)
            newError = minimizationFunction(monteCarloMatrix,testPointArray[point],newSize)
            #oldError = minimizationFunction(monteCarloMatrix,center,newSize)
            #print('error: ', error, 'radius: ', radius)
            # compare error 
            if newError < error:
                #print(newError,radius)
                error = newError
                center = testPointArray[point]
            # if the radius gets too small without converging, 
            if radius < tolerance:
                #radius = 2*newSize*error
                radius = error
        # reduce the search radius 
        #radius = radius / 2
        #radius = newSize*error
        radius *= 0.9
                
        # failsafe
        counter += 1
        if counter > 100000:
            break
    # display results
    print('finding the right column with Monte Carlo method')
    print('radius: ' , radius)
    print('error: ' , error)
    print('number of loops: ', counter)
    print(monteCarloMatrix)
    
'''The last function computes an estimate for an unknown y value at 
the last interval. The explanation behind this algorithm is somewhat
lengthy, and is included in the presentation write up.  '''
def extrapolation():
    # create the ratio vector
    ratioVector = readMatrix(n+1,'RatioMatrix.txt')[0,0:n+1]
    ratioVector = np.matrix(ratioVector)
    ratioVector = ratioVector.getT()
    # create truncated Pascal
    P = pascal(n+1)[0:-1,:]
    # set up y vector
    yVector = np.matrix(y)
    # pfit of known values
    matrix = readMatrix(n, 'RatioMatrix.txt')[0:n,0:n]
    coefficients = createCoefficientArray(n,matrix)
    # midpoint
    m = n*d / 2
    # compute yn
    truncated = yVector*P*ratioVector
    pfit = cascade(m,coefficients)
    yNext = (-1)**(n+1)*(pfit-truncated)
    xNext = (n+1)*d
    # plotting the prediciton
    top = float(max([max(y),yNext]))
    bottom = float(min([min(y),yNext]))
    plt.ylim(bottom - 1.5,top + 1.5)
    plt.plot(xNext,yNext,'c*')
    print(yNext)
    # arrow annotation
    plt.annotate('extrapolation',
                 xy=(xNext,float(yNext)), 
                 xytext=(xNext-10,float(yNext)+0.5), 
                 arrowprops=dict(facecolor='black', shrink=0.05,width = 0.1, headwidth = 4))
    

#extrapolation()
#monteCarlo(4)
display(n)
plt.show()

runtime()




