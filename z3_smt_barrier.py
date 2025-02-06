# Computes IBC and plot them for 1d system using SMT solver
from z3 import * # SMT solver
import numpy as np
from matplotlib import pyplot as plt
import time

def f(x):
    return(x/2)

# monomial for linear barrier certificate template
def ret_monoms(x):
    return(1, x)

# barrier certificate template - linear
def Bar(i,x,coeffs_smt): #B_i, x - states
    return(np.dot(ret_monoms(x),coeffs_smt[i]))

def realToFloat(var): #z3 real to python float
    var = var.as_fraction()
    return(float(var.numerator) / float(var.denominator))

leng = len(ret_monoms(0)) # length of monomial terms
dim = 1 # system state dimension

# IBC computation. kmax = 0 -> standard BC
def IBC(kmax, Xo, Xu, X, length = leng):
    # kmax - max k value for B_k in IBC, Xo - initial, Xu - unsafe, X - state set, length - hyperparameter pass
    coeffs_smt = [] # store coefficients for all B_i
    for i in range(kmax+1): # for each B_i upto max k
        t = [] # store coefficients for each B_i
        for j in range(length): # for each term in monomial
            c = Real("coeff"+str(i)+","+str(j)) # construct a coefficient
            t.append(c)
        coeffs_smt.append(t)
    
    xvar = [] # state variables
    for i in range(dim):
        xvar.append(Real("x"+str(i)))

    # initial state
    initial_cond = ForAll(xvar, 
                    Implies(And(Xo[0] <= xvar[0], xvar[0] <= Xo[1]), # if initial state
                            Bar(0,xvar[0],coeffs_smt) <= 0)) # B_0 condition

    # unsafe state
    quer = And(True)
    for i in range(kmax+1): # for each B_i
        quer = And(quer, Bar(i,xvar[0],coeffs_smt) > 0) # condition on B_i

    unsafe_cond = ForAll(xvar, 
                         Implies(And(Xu[0] <= xvar[0], xvar[0] <= Xu[1]), #if unsafe
                                quer)) # then check conditions on all B_i

    # 3rd condition
    quer = And(True)
    for j in range(kmax): # 3rd condition B_i(x) => B_(i+1)(f(x))
        last_cond = ForAll(xvar,
                    Implies(And(X[0] <= xvar[0], xvar[0] <= X[1], Bar(j, xvar[0], coeffs_smt) <= 0),
                        Bar(j+1, f(xvar[0]), coeffs_smt) <= 0))
        quer = And(quer, last_cond)
    
    # 4th (last) condition B(k) => B(k)
    last_cond = ForAll(xvar, 
                    Implies(And(X[0] <= xvar[0], xvar[0] <= X[1], Bar(kmax, xvar[0], coeffs_smt) <=0),
                        Bar(kmax, f(xvar[0]), coeffs_smt) <=0) )
    last_cond = And(quer, last_cond)
    
    # gather all conditions
    Bar_cond = And(initial_cond, unsafe_cond, last_cond)
    
    # compile everything and solve
    s = Solver()
    s.add(Bar_cond)

    print(f"k = {kmax}:")
    start_time = time.time()
    status = s.check() # z3 expecting to run s.check() only once so save it in variable
    print(f"Status: {status}")
    print(f"Execution time: {time.time() - start_time}") # ~ 643s with only first two conditions
    if (status == sat):
        return(s.model())
    elif (status == unsat):
        # c = s.unsat_core()
        return(None)

k_max = 3 # max bound
Xo = [2, 2.3] # initial set
Xu = [1.6,1.9] # unsafe set
X = [0,3] # state set
colors = ['blue','m','g','r','b','k',"pink","olive","orange","purple"]
for k in range(k_max+1): # goes from 0 to k_max
    model = IBC(k,Xo,Xu,X)
    if model != None: # if model successfully computed
        state = np.linspace(X[0],X[1],100)
        B = np.zeros((k+1,len(state))) # need +1 because of 0 indexing
        # sort model coefficients alphabetically by name (ascendingly for each i in B_i)
        sorted_model = sorted([(var, model[var]) for var in model], key = lambda x: str(x[0]))
        start = 0 # starting index for each B_i
        for i in range(k+1): # for i from 0 upto and including k
            # collect relevant coefficients for each B_i (start to (not including) start+leng gives the right number of monomial terms)
            # and reverse it to get highest to lowest degree ([::-1])
            coeffs = np.array([realToFloat(val) for (_, val) in sorted_model[start:start+leng]])[::-1]
            start += leng # shift starting index
            # coeffs needs to be highest to lowest degree i.e. coeffs[0]*x + coeffs[1]*1.
            print(f"i={i}: {coeffs}")
            B[i,:] = np.polyval(coeffs, state) # evaluate barrier certificate
            plt.plot(state, B[i,:], label="B"+str(i), color=colors[i]) # plot barrier certificate
            # plot 0-sublevel set of barriers (B(x) <= 0). y1, y2 decided based on plot without them
            plt.fill_between(state, y1 = -1.5625, y2 = 2, where= (B[i,:]<=0), facecolor=colors[-i-1], alpha=.2)
        
        # plot initial set region (only x axis)
        plt.fill_between(state, np.min(B), np.max(B), where= np.logical_and(Xo[0] <= state, state <= Xo[1]), facecolor='green', alpha=.3,label="Initial")
        # plot unsafe set region (only x axis)
        plt.fill_between(state, np.min(B), np.max(B), where= np.logical_and(Xu[0] <= state, state <= Xu[1]), facecolor='red', alpha=.3,label="Unsafe")
        plt.axhline(0,color='k') # horizontal line at B = 0
        plt.grid(True) 
        plt.legend(loc="upper left")
        plt.xlabel("x(t)")
        plt.savefig(f"./media/ibc_k{k}.svg", bbox_inches='tight', format='svg', dpi=1200)
        plt.show()
        break # once IBC found, no need to check or plot for larger k
