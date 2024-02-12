# Examples of "Tutorial Kwant 1.4.0"
# -------------------------------------------------
# The quantum wire model: Spin - Orbit
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
# Import for Matrix construction
#
import tinyarray
#
# Define the Pauli matrices
sigma_0 = tinyarray.array([[1, 0], [0, 1]]) # Identity!!!
sigma_x = tinyarray.array([[0, 1], [1,0]])
sigma_y = tinyarray.array([[0, -1j], [1j, 0]])
sigma_z = tinyarray.array([[1, 0], [0, -1]])
#
# ===========================================================================================================================
# Now we start define the system with this values, an "empty square-lattice"
#
def make_system(t=1.0, alpha=0.5, e_z=0.08, W=10, L=30):
    lat = kwant.lattice.square() # see that it is not a square lattice of "a", a=1
    
    syst = kwant.Builder()
#
# ===========================================================================================================================
# As difference with respect to Quantum wire model_2, we define in "kwant.Builder" with matrix term and not numbers as before:
#    
    syst[(lat(x, y) for x in range(L) for y in range(W))] = \
        4 * t * sigma_0 + e_z * sigma_z
# The hopping in x-direction
    syst[kwant.builder.HoppingKind((1, 0), lat, lat)] = \
        -t * sigma_0 + 1j * alpha * sigma_y / 2
# The hopping in y-direction
    syst[kwant.builder.HoppingKind((0, 1), lat, lat)] = \
        -t * sigma_0 - 1j * alpha * sigma_x / 2
# ===========================================================================================================================
# Now construct the lead
#
    lead = kwant.Builder(kwant.TranslationalSymmetry((-1, 0)))
    
    lead[(lat(0, j) for j in range(W))] = 4 * t * sigma_0 + e_z * sigma_z
# Hopping in x-direction
    lead[kwant.builder.HoppingKind((1, 0), lat, lat)] = \
        -t * sigma_0 + 1j * alpha * sigma_y / 2
# Hopping in y-direction
    lead[kwant.builder.HoppingKind((0,1), lat, lat)] = \
        -t * sigma_0 - 1j * alpha * sigma_x / 2
    
    syst.attach_lead(lead) # paste the lead in left
    syst.attach_lead(lead.reversed()) # By reversed and translational sym, we join the right lead.
    
    return syst # return system

# In this step, compute the conductance
def plot_conductance(syst, energies):
    data = []
    for energy in energies:
        smatrix =  kwant.smatrix(syst, energy)
        data.append(smatrix.transmission(1, 0))
        
# To show the plot of conductance
    pyplot.figure()
    pyplot.plot(energies, data)
    pyplot.xlabel("Energy [t]")
    pyplot.ylabel("Conductance [e^2/h]")
    pyplot.show()

def main():
    syst = make_system()

    kwant.plot(syst) # to appreciate the lattice
    
    syst = syst.finalized() #finalize the system
    
    plot_conductance(syst, energies = [0.01 * i - 0.3 for i in range(100)])
    
    
if __name__ == '__main__':
    main()
