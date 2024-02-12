# Examples of "Tutorial Kwant 1.4.0"
# -------------------------------------------------
# The quantum wire model
# -------------------------------------------------
# Mario Solís
# PhD student in Physics
# Grupo de Partículas y Campos
# Instituto Balseiro - Centro Atómico Bariloche
# m.f.solis.benites@gmail.com
# ===========================================================================================================================
# Call Kwant!
#
import kwant
#
# Import the Plot lybrary of matplotlib, obviously for plotting!
#
from matplotlib import pyplot
#
# ===========================================================================================================================
# Now we define the system "syst", to this objective, we use the "kwant.Builder"
#
syst = kwant.Builder()
#
# define the system: square lattice "lat" with lattice spacing "a" equal to 1, by simplicity,
#
a = 1
lat = kwant.lattice.square(a)
#
# ===========================================================================================================================
# Now we define the size of the system
#
t = 1.0 # h\bar^2 /(2 m a^2)
W = 10 # width
L = 30 #long
for i in range(L):
    for j in range(W):
        # On site Hamiltonian
        syst[lat(i,j)] = 4*t
        
        # Hopping in y-direction
        if j>0:
            syst[lat(i,j), lat(i,j-1)] = -t

        # Hopping in x-direction
        if i>0:
            syst[lat(i,j), lat(i-1,j)] = -t
#
# ===========================================================================================================================
# Now we define the Leads, they will construct with Builder comand, we start with left lead.
# REMIND: the system have transaltional symmentry!
#
#sym_left_lead = kwant.TranslationalSymmetry((-a, 0)) # the value "-a" must point in a direction away from scattering region
#left_lead = kwant.Builder(sym_left_lead) # construct the lead.
#
# For a lead in y-direction
#
#for j in range(W):
#   left_lead[lat(0, j)] = 4*t
#   if j > 0:
#        left_lead[lat(0, j), lat(0, j-1)] = -t
#   left_lead[lat(1, j), lat(0,j)] = -t
#
# Now we use "syst.attach_lead()" to join the left lead to the system in correct position
#
#syst.attach_lead(left_lead)
#
# Thus, in similar way, we construct the right lead. Be aware about that in this case the vector must point in the other direction "a"
#
#sym_right_lead = kwant.TranslationalSymmetry((a, 0)) # the value "a" must point in a direction away from scattering region
#right_lead = kwant.Builder(sym_right_lead) # construct the lead.
#
# For a lead in y-direction
#
#for j in range(W):
#   right_lead[lat(0, j)] = 4*t
#   if j > 0:
#        right_lead[lat(0, j), lat(0, j-1)] = -t
#   right_lead[lat(1, j), lat(0,j)] = -t
#
#syst.attach_lead(right_lead)
#
# ===========================================================================================================================
# To plot construction the system! "kwant.plot(syst)"
#
kwant.plot(syst)
#
# ===========================================================================================================================
# To FINISH the construction of our system we use "syst.finalized()". Then we allow to start the transport computation
#
#syst = syst.finalized()
#
# ===========================================================================================================================
# ===========================================================================================================================
# ===========================================================================================================================
# Now compute the conductance!
#
#energies = [] # define the empty value of energies
#data = [] # empty variable for the data
#for ie  in range(100):
#    energy = ie*0.01

    # compute the scattering matrix at a given energy
#    smatrix = kwant.smatrix(syst, energy) # compute the scattering matrix and solve the sparse linear system.

    #compute the transmission probability from lead 0 to lead 1
#    energies.append(energy)
#    data.append(smatrix.transmission(1, 0)) # coment below!
    # Warning! the number on "smatrix.transmission", must be the same when we defined the "attach.lead" in previous part.
#
# ===========================================================================================================================  
# Finally, plot the result of the conductance
#
#pyplot.figure()
#pyplot.plot(energies, data)
#pyplot.xlabel("Energy [t]")
#pyplot.ylabel("Conductance [e^2/h]")
#pyplot.show()











