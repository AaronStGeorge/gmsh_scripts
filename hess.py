from dolfin import *
import pylab as p

mesh = '../data/mesh.xml'
ub   = '../data/ubmag.xml'

mesh = Mesh(mesh)

V    = FunctionSpace(mesh, "CG", 1)
U    = Function(V)
File(ub) >> U

# compute the Hessian 
Hxx  = TrialFunction(V)
Hxy  = TrialFunction(V)
Hyy  = TrialFunction(V)
phi  = TestFunction(V)

a_xx = Hxx * phi * dx
L_xx = - U.dx(0) * phi.dx(0) * dx

a_xy = Hxy * phi * dx
L_xy = - U.dx(0) * phi.dx(1) * dx

a_yy = Hyy * phi * dx
L_yy = - U.dx(1) * phi.dx(1) * dx

Hxx  = Function(V)
Hxy  = Function(V)
Hyy  = Function(V)

solve(a_xx == L_xx, Hxx)
solve(a_xy == L_xy, Hxy)
solve(a_yy == L_yy, Hyy)

# Get the value of Hessian at each node
Hxx_a = project(Hxx,V).vector().array()
Hxy_a = project(Hxy,V).vector().array()
Hyy_a = project(Hyy,V).vector().array()

# Find the |dominant eigenvalue of the Hessian|
m = p.array([Hxx_a,Hxy_a,Hxy_a,Hyy_a])
m = p.array( [max( abs( p.eig( m[:,i].reshape((2,2)))[0] )) for i in \
        range(len(Hxx_a))] )
m = p.log(m)

# Read into function space
H_lambda = Function(V)
H_lambda.vector().set_local(m)

# Output data
output = File('hess_py_output.pvd')
output << H_lambda
