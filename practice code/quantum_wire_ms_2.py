# Examples of "Tutorial Kwant 1.4.0"
#----------------------------------------------------------
# The quantum wire model - More compact!
# ---------------------------------------------------------
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
# Now we start define the system with this values, an "empty square-lattice"
#
def make_system(a=1, t=1.0, W=6, L=10):
    lat = kwant.lattice.square(a)
    
    syst = kwant.Builder()
    
    syst[(lat(x, y) for x in range(L) for y in range(W))] = 4*t # this construct the lattice points in one line!
    syst[lat.neighbors()] = -t # Hopping!
# ===========================================================================================================================
# Now construct the lead similar and before code.
#
#    lead = kwant.Builder(kwant.TranslationalSymmetry((-a, 0)))
#    lead[(lat(0, j) for j in range(W))] = 4*t
#    lead[lat.neighbors()] = -t
#    syst.attach_lead(lead) # paste the lead in left
#    syst.attach_lead(lead.reversed()) # By reversed and translational sym, we join the right lead.
    return syst # return system

def make_system(a=1, t=2.0, W=6, L=10):
    lat=kwant.lattice.square(a)

    syst=kwant.Builder()

    syst[(lat(x, y) for x in range(L) for y in range(W))] = 4*t # this construct the lattice points in one line!
    syst[lat.neighbors()] = -t # Hopping!
    return syst
# ===========================================================================================================================
# In this step, compute the conductance
#
#def plot_conductance(syst, energies):
#    data = []
#    for energy in energies:
#        smatrix =  kwant.smatrix(syst, energy)
#        data.append(smatrix.transmission(1, 0))

#    pyplot.figure()
#    pyplot.plot(energies, data)
#    pyplot.xlabel("Energy [t]")
#    pyplot.ylabel("Conductance [e^2/h]")
#    pyplot.show()

def main():
    syst = make_system()

    kwant.plot(syst)
    syst = syst.finalized()
#    plot_conductance(syst, energies = [0.01*i for i in range(100)])
    
if __name__ == '__main__':
    main()
