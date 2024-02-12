# -------------------------------------------------
# The Weyl semimetal model 
# -------------------------------------------------
# Mario Solís
# PhD student in Physics
# Grupo de Partículas y Campos
# Instituto Balseiro - Centro Atómico Bariloche
# m.f.solis.benites@gmail.com
# ===========================================================================================================================
#
import kwant
import math
#from math import pi, sqrt, tanh
import cmath
import scipy.sparse.linalg as sla
import numpy
from matplotlib import pyplot # For plotting
import tinyarray # For matrix construction
#
# ===========================================================================================================================
# Define the Pauli matrices
#
sigma_0 = tinyarray.array([[1, 0], [0, 1]]) # Identity!
sigma_x = tinyarray.array([[0, 1], [1, 0]])
sigma_y = tinyarray.array([[0, -1j], [1j, 0]])
sigma_z = tinyarray.array([[1, 0], [0, -1]])
#
# We define values of variables
#
t = 1.0
v = 1.0
L = 13    # length y
W = 13   # length x and  z
k0 = math.pi/10.
#k0 = -0.4636476090008061
m0 = math.cos(k0)
vtilde = v/math.sin(k0)
h = vtilde*m0 + 2*v
a=1
#
# ===========================================================================================================================
#
def make_system(t, v, L, W, k0, m0, vtilde, h, a):
    lat = kwant.lattice.cubic(a)
    
    syst = kwant.Builder()
#
# ===========================================================================================================================
#
for i in range(W):
    for j in range(L):
       for k in range(W):

             # On-site Hamiltonian
            syst[lat(i, j, k)] = h*sigma_x

             # Hopping in x-direction
            if i > 0:
                syst[lat(i, j, k), lat(i - 1, j, k)] = -(vtilde/2)*sigma_x

              # Hopping in y-direction
            if j > 0:
               syst[lat(i, j, k), lat(i, j - 1, k)] = -(v/2)*((1j)*sigma_y + sigma_x)

              # Hopping in z-direction
            if k > 0:
                 syst[lat(i, j, k), lat(i, j, k - 1)] =-(v/2)*((1j)*sigma_z + sigma_x)


#  I close hopping in x-direction, connecting W-1 and 0
for j in range(L):
    for k in range(W):
        syst[lat(0, j, k), lat(W-1, j, k)] = -(vtilde/2)*sigma_x

# I close the hopping in z-direction, connecting W-1 and 0
for i in range(W):
    for j in range(L):
        syst[lat(i, j, 0), lat(i, j, W-1)] = -(v/2)*((1j)*sigma_z + sigma_x)

#
# ===========================================================================================================================
# Now construct the leads
# We define the right lead, along y
#sym_right_lead = kwant.TranslationalSymmetry((0, a, 0))
#right_lead = kwant.Builder(sym_right_lead)


#for i in range(W):
  #  for k in range(W):


    #    right_lead[lat(i, 0, k)] = m0*sigma_0

      #  if i > 0:
          #  right_lead[lat(i, 0, k), lat(i-1, 0, k)] = -t*sigma_0

        #if k > 0:
        #    right_lead[lat(i, 0, k), lat(i, 0, k-1)] = -t*sigma_0


        #right_lead[lat(i, 1, k), lat(i, 0, k)] = -t*sigma_0


# I close hopping in x-direction, connecting  W-1 and 0
#for k in range(W):
  #  right_lead[lat(0, 0, k), lat(W-1, 0, k)] = -t*sigma_0

 # I close hopping in z-direction, connecting  W-1
#for i in range(W):
    #right_lead[lat(i, 0, 0), lat(i, 0, W-1)] = -t*sigma_0
    
    #syst.attach_lead(right_lead)
    #left_lead = right_lead.reversed()
    #syst.attach_lead(left_lead)

     #return syst # return system
#
# ===========================================================================================================================
#
def main():
    #lead = make_lead().finalized()
    #kwant.plotter.bands(lead, show=False) # To use the bands computation
    #pyplot.xlabel("Momentum [(lattice constant)^-1")
    #pyplot.ylabel("Energy [t]")
    #pyplot.show()
    syst = make_system()

    kwant.plot(syst)

    syst=syst.finalized()

if __name__ == '__main__':
    main()
