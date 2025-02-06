# vector field
@polyvar x
vars = [x]

# simple dynamics. same as the z3_smt_barrier example
f = [x/2]

#initial, unsafe sets and state space
Xo = [2, 2.3] # initial set
Xu = [1.6,1.9] # unsafe set
X = [0,3] # state space

# semi-algebraic set description polynomial of the initial/unsafe set and state space, >=0 by default
go = [x-Xo[1], Xo[2]-x] #initial
gu = [x-Xu[1], Xu[2]-x] #unsafe
g = [x-X[1], X[2]-x] #whole
