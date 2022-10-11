from scipy import optimize as opt
import numpy as np

def test(x):
    
    jibe = [2*x[0] - 4*x[1],9 - x[0] - x[1]]
    return jibe

x0 = np.array([1,2])


roots = opt.root(test,x0)
print(roots)