# vector field - system variables
@polyvar v p
vars = [v,p]

# constants
alpha = 1.1; #growth rate -> prey
beta  = 0.4; #death rate -> prey
gamma = 0.4; #death rate -> predator
delta = 0.1; #growth rate -> predator
T = 0.1; # sampling time

# system dynamics - discrete-time Lotka-Volterra-type model
f = [v + T*(alpha*v*(1-v) - beta*v*p),
     p + T*(gamma*p*(p*0-1) + delta*v*p)]

#initial, unsafe sets and state space
Xo = [6,7,2,3] # initial set (v,p)
Xu = [3,5,0,3] # unsafe set
X = [0,10,0,5] # state space

# semi-algebraic set description polynomial of the initial/unsafe set and state space, >=0 by default
go = [v-Xo[1], Xo[2]-v, p-Xo[3], Xo[4]-p] #initial
gu = [v-Xu[1], Xu[2]-v, p-Xu[3], Xu[4]-p] #unsafe
g  = [v-X[1],  X[2]-v,  p-X[3],  X[4]-p] #whole
