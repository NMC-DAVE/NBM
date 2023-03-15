import numpy as np

np = 2000
q = [ np//20, np//10, (3*np)//20, np//5, np//4, (3*np)//10, \
    (7*np)//20, (2*np)//5, (9*np)//20, np//2, (11*np)//20, (3*np)//5, (13*np)//20, \
    (7*np)//10, (3*np)//4, (4*np)//5, (17*np)//20, (9*np)//10, (19*np)//20]
            
e = np.array(q,dtype=float) / float(np) 
print (e)