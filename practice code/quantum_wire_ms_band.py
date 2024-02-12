# Examples of "Tutorial Kwant 1.4.0"
# -------------------------------------------------
# The quantum wire model - Band structure
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
# Now we start define the system with this values, an empty square-lattice 
#

def make_lead(a=1, t=1.0, W=5):
    lat = kwant.lattice.square(a)
    
    sym_lead = kwant.TranslationalSymmetry((-a,0)) # We working with one Lead, there is not the scattering region.
    lead = kwant.Builder(sym_lead)
    
    # Define the lead and the next hopping term
    for j in range(W):
        lead[lat(0, j)] = 4*t

        if j>0:
            lead[lat(0, j), lat(0, j-1)] = -t

        lead[lat(1, j), lat(0, j)] = -t

    return lead

# Then, we finalize
def main():
    lead = make_lead().finalized()
    kwant.plotter.bands(lead, show=False) # To use the bands computation
    pyplot.xlabel("Momentum [(lattice constant)^-1")
    pyplot.ylabel("Energy [t]")
    pyplot.show()

    
if __name__ == '__main__':
    main()
