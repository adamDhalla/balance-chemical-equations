# balance-chemical-equations by adam dhalla (adamdhalla.com)

A rough but functional way to balance chemical equations that do NOT undergo polyatomic decomposition
(polyatomics must be in same groups on left hand side and right hand side). 

#---------------------------------------------------------------------------------------------------------------INPUT SYNTAX
There is a standard input syntax for all chemical equations. We take an input to our single function, solveEquation(string).

For regular: AXBY + CZ -> CZAX + BY 
For polyatomics: (AB)X + CY -> (AB)Y + C
______________________________________________
Where A, B, C are elements (e.g Au, C, Pb)
Where X, Y, Z are subscripts (how much corresponding element occurs)

IMPORTANT I: Please, if you are dealing with an element that only occurs once, e.g Na or Cl in NaCl 
             Instead input (Na1Cl1). Always explicitly refer to the subscript in the input.

IMPORTANT II: For Polyatomics, don't include the subscript inside the polyatomic - Don't Do (NO3), just (NO)

Sample inputs: solveEquation("Au3 + Br2 -> Au2Br5") 
               solveEquation("(CrO)4 + Pb5N2 -> Pb2(CrO)4 + N3")

#---------------------------------------------------------------------------------------------------------------UNDER THE HOOD
For a complete guide on how this algorithm works, I wrote a detailed medium post + diagrams, picture and bckg information: 
https://adamdhalla.medium.com/balancing-chemical-equations-using-python-4b9086a92a7c

This system exploits the fact that a chemical equation can be made into a linear system of equations. For each term in the chemical 
equation, I assign a vector that records which elements (and how many) are present in each term. Each of these vectors multiplies
by the coefficient (x1, x2, ... xn) that we are solving for.

Then, we combine these vectors into a matrix A (move all vectors to left-hand side) and solve for the nullspace of A. 

The problem is that using a nullspace basis of A as an answer will probably not work. It will most likely include zero as an element 
(we cannot have 0 coefficients when balancing), negatives (can't have negative molecules) or decimals (can't have fractions of mole-
cules). So, we need a natural-numbered, integer-valued vector, that is also as small as possible (we don't want multiples of our answer).

Since there is no effective way to find the smallest positive, integer nullspace vectors, we resort to a brute force method. We need to 
systematically 'try out' a range of values for each coefficient to the nullspace vectors, take the combinatinos, keep the viable ones 
(natural-numbered, integer-valued) in a list of vectors. 

Then, we take the vector with the smallest L2 norm in the list, which is our correct coefficient vector - the smallest, natural-numbered, 
integer-valued set of coefficients for our problem. 

Then, I display these back in the format of the chemical equations.



Then, I try to solve for the nullspace of a systems of equations matrix A. 
