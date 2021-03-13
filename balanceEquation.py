# chem equations solver 
import numpy as np
import numpy.linalg as lin 
import re
import sympy
import itertools

#------------------------------------------------------------------------------ Helper Functions
def join(l, sep):
  """from user fabian789 on Code Review Stack Exchange: separates a list into a string with
  as passed in separator""" 
  
    out_str = ''
    for i, el in enumerate(l):
        out_str += '{}{}'.format(el, sep)
    return out_str[:-len(sep)]

def removeNums(l):
   """used to remove all numeric values from list"""
   return [value for value in l if not value.isdigit()]

def removeBlanks(l):
   """used to remove all blank '' characters from a list"""
   return [value for value in l if value != '']

# make re equation for finding all elements in equation: 
def findalphs(string, removenums=False, removeblanks=False):
    """passes in a string, separates into a list of homogenous type substrings. 

    e.g "fooptown524happy" -> ["fooptown", "524", "happy"]

        string:         a string you want segmented.
        removenums:     optionally remove all nums from output list
                        (defaults to FALSE)
        removeblanks:   optionally remove all blanks ('') from output list
                        (defaults to FALSE)
    
    """
    # if it is a single character (one instance of one element) just return it
    if len(string) == 1:
        return list(string)

    separates = re.split('(\d+)', string)
    
    if removenums == True:
        separates = removeNums(separates)
    if removeblanks == True:
        separates = removeBlanks(separates)
    
    return separates

#-----------------------------------------------------------------Main Function (see README for input specs)

def solveEquation(equation):

    # split equation in to lhs and rhs
    splitter = equation.split(" -> ")
    lhs = splitter[0]
    rhs = splitter[1]
    #print(lhs, rhs)

    # split into terms
    lsplits = lhs.split(" + ")
    ltermnum = len(lsplits)

    rsplits = rhs.split(" + ")
    rtermnum = len(rsplits)

    vars = np.zeros((1, (rtermnum + ltermnum)))
    #print(vars)

    allelements = []

    #find how many elements there are
    for term in lsplits:
        alphs = findalphs(term, removeblanks=True, removenums=True)
        #print(alphs)
        allelements = allelements + alphs
    
    uniqueelements = list(set(allelements))
    elementdict = {key: val for val, key in enumerate(uniqueelements)}
    #print(elementdict)

    elementamt = len(uniqueelements)
    
    # if no polyatomics
    if True:
        lhsvectors = []
        rhsvectors = [] 

        # lhs first, get vectors
        for term in lsplits:

            vec = np.zeros((elementamt, 1))

            combos = findalphs(term, removeblanks=True, removenums=False)

            # since each term in lsplits might have several occurences of diff. elems, (e.g Au3Pb4)
            # iterate through each element in each term and log each term's elements and occurences in vector
            termelems = [] 
            for i in range(0, len(combos), 2):
                currentelem = combos[i]
                currentamt = int(combos[i + 1])
                termelems.append((currentelem, currentamt))

            # now, loop over element / occurence pairs in the term (might be more than 1)
            # combo in format ('Au', 3), ('Pb', 4)
            for combo in termelems:
                vec[elementdict[combo[0]]] = combo[1]

            lhsvectors.append(vec)

        # now, rhs
        for term in rsplits:

            vec = np.zeros((elementamt, 1))
            combos = findalphs(term, removeblanks=True, removenums=False)
            termelems = [] 

            for i in range(0, len(combos), 2):
                currentelem = combos[i]
                currentamt = int(combos[i + 1])
                termelems.append((currentelem, currentamt))
            #print(termelems)

            # now, loop over element / occurence pairs in the term (might be more than 1)
            # combo in format ('Au', 3), ('Pb', 4)
            for combo in termelems:
                vec[elementdict[combo[0]]] = combo[1]

            rhsvectors.append((-1*vec))
            

        # since vectors on rhs will be put over the equal sign to create system and solve for 0:
        
        from sympy import Matrix
        A = Matrix( (np.concatenate((lhsvectors + rhsvectors), axis=1)) )

        b = Matrix( np.zeros((np.shape(A)[0], 1)) )
        
        # calculate nullspace vectors
        x = A.nullspace()
        N  = ((np.array(x)).transpose().astype('float'))
        print(N)

        varnum = np.shape(N)[1]

        # create list of integers to multiply each vector by, and take the one with smallest norm, int and nonneg
        trynums = list(range(1, 20))
        alltries = trynums*varnum

        # taking the combinations of this list makes it so that each variable gets to be one value at one time
        combinations = list(itertools.combinations(alltries, varnum))
        
        candidateVectors = []
        for vars in combinations:
            isDec = None
            # do Nx = 0 where N is a matrix of nullspace vectors and x is each of our combinations of vars
            x = np.asarray(vars).reshape((varnum, 1))
            vector = np.dot(N, x)

            # get rid of anything with decimals
            for item in vector:
                item = item[0]

                if not (item).is_integer():
                    isDec = True
                    break
            
            if isDec == True:
                continue 
            else:
                candidateVectors.append(vector)

        
        bestvector = [lin.norm(candidateVectors[0]), candidateVectors[0]]
        for vec in candidateVectors:
            norm = lin.norm(vec)

            if norm < bestvector[0]:
                bestvector[0] = norm
                bestvector[1] = vec
    

    coefficients = (bestvector[1].flatten().astype('int32')).tolist()


    valinbest = 0 # which coefficient to put in front of term

    lsplitsstatic = lsplits.copy()

    for term in lsplitsstatic:
        insertterm = lsplits.index(term)
        coefficient =  int(str(coefficients[valinbest]))
        lsplits.insert(insertterm, coefficient)
        valinbest += 1
    
    rsplitsstatic = rsplits.copy()
    for term in rsplitsstatic:
        insertterm = rsplits.index(term)
        coefficient =  int(str(coefficients[valinbest]))
        rsplits.insert(insertterm, coefficient)
        valinbest += 1
    

    # add coefficients
    lcompounds = []
    for compound in range(0, len(lsplits), 2):
        lcompounds.append(str(lsplits[compound]) + " " + str(lsplits[compound + 1]))

    rcompounds = []
    for compound in range(0, len(rsplits), 2):
        rcompounds.append(str(rsplits[compound]) + " " + str(rsplits[compound + 1]))

    # you can turn this into a void function if you'd like and print output instead of returning it
    return (join(lcompounds, ' + '), " -> ", join(rcompounds, ' + '))
