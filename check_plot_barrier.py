# check if barrier from Julia satisfies IBC conditions and if so, plot
# Julia might give wrong answer occassionally due to numerical approximation
from z3 import * # SMT solver
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import time

# Barrier for lotka-volterra model from paper
def Bar(j,x): 
    # j - index for B_j, x - states
    v, p = x[0], x[1]
    if j == 0: # B_j
        return (0.13305*v**3 + 0.16826*v**2*p - 0.11002*v*p**2 + 0.30322*p**3 - 1.52807*v**2 - 3.55885*v*p + 2.1532*p**2 + 5.77075*v + 3.22429*p + 5.34744)
    else:
        return (0.24199*v**3 + 0.22602*v**2*p + 0.14704*v*p**2 + 0.1285*p**3 + 2.99671*v**2 - 0.13358*v*p - 0.5105*p**2 - 6.71922*v + 1.58361*p - 12.69347) 

# system dynamics
def f(x):
    # x - states
    v, p = x[0], x[1]

    # constants
    alpha = 1.1; #growth rate -> prey
    beta  = 0.4; #death rate -> prey
    gamma = 0.4; #death rate -> predator
    delta = 0.1; #growth rate -> predator
    T = 0.1; # sampling time

    return([v + T*(alpha*v*(1-v) - beta*v*p),
            p + T*(gamma*p*(p*0-1) + delta*v*p)])

xvar = []
for i in range(2): #2 states
    xvar.append(Real("x"+str(i)))

# initial, unsafe sets and state space
# First 2 for v, last 2 for p
Xo = [6,7,2,3] # initial set (v,p).
Xu = [3,5,0,3] # unsafe set
X = [0,10,0,5] # state space

# initial states
init_states = And(Xo[0] <= xvar[0], xvar[0] <= Xo[1], Xo[2] <=xvar[1], xvar[1] <= Xo[3]) # variable range
init_cond = (Bar(0,xvar)<= 0) # condition on barrier B_0
init_check = ForAll(xvar, Implies(init_states, init_cond)) # condition to check

# unsafe_states
unsafe_states = And(Xu[0] <=xvar[0], xvar[0] <= Xu[1], Xu[2] <=xvar[1], xvar[1] <= Xu[3])
unsafe_cond = And(Bar(0,xvar) > 0, Bar(1,xvar) > 0) # condition on all barriers for unsafe
unsafe_check = ForAll(xvar, Implies(unsafe_states, unsafe_cond))

# 3rd condition B_0(x) <=0 => B_1(f(x))<=0
all_states = And(X[0] <= xvar[0], xvar[0] <= X[1], X[2] <=xvar[1], xvar[1] <= X[3])
third_cond = And(all_states, Bar(0, xvar) <= 0)
third_check = ForAll(xvar, Implies(third_cond, Bar(1, f(xvar)) <= 0))

# 4th (last) condition B_1(x) => B_1(f(x))
last_cond = And(all_states, Bar(1, xvar) <= 0)
last_check = ForAll(xvar, Implies(last_cond, Bar(1, f(xvar)) <= 0))

Bar_cond = And(init_check,unsafe_check,third_check, last_check) # IBC schedule

s = Solver()
# s.add(init_check)
# s.add(unsafe_check)
# s.add(third_check) # sometimes takes too long with this and next line (i.e. all 4)
# s.add(last_check)
s.add(Bar_cond) # sometimes takes too long

print("Solving...")
start_time = time.time()
status = s.check() # check if conditions are satisfied
print(f"Status: {status}")
print(f"Execution time: {time.time() - start_time}") # ~ 250s with all conditions
if (status == sat): # plot if barrier works
    n = 100
    v,p = np.linspace(X[0],X[1],n+1), np.linspace(X[2],X[3],n+1) # state space
    x,y = np.meshgrid(v,p)
    fig, ax = plt.subplots()
    
    # filled contour for 0-sublevel set of barriers (B(x) <= 0)
    ax.contourf(v, p, Bar(0,[x,y]),levels=0,colors=['blue','w'],alpha=0.9)
    ax.add_patch(Rectangle((6,2), 0, 0,color = 'blue',alpha=0.45,label="B0")) # added for legend, no rectangle plotted
    ax.contourf(v, p, Bar(1,[x,y]),levels=0,colors=['orange','w'],alpha=0.5)
    ax.add_patch(Rectangle((0,0), 0, 0,color = 'orange',alpha=0.5,label="B1")) # added for legend, no rectangle plotted
    
    # plot initial and unsafe regions
    ax.add_patch(Rectangle((Xo[0], Xo[2]), Xo[1]-Xo[0], Xo[3]-Xo[2],color = '#90ee90',alpha=0.8,label="Initial"))
    ax.add_patch(Rectangle((Xu[0], Xu[2]), Xu[1]-Xu[0], Xu[3]-Xu[2],color = 'red',alpha=0.5,label="Unsafe"))
    ax.axis('scaled')
    ax.grid()
    ax.legend(loc="lower right")
    plt.xlabel("$v$")
    plt.ylabel("$p$")
    # transparency issue with direct saving eps
    plt.savefig("./media/ibc_2d_k1.svg", bbox_inches='tight', format='svg', dpi=1200) # save as svg, convert to eps online
    plt.show()