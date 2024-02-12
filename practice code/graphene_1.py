# Examples of "Tutorial Kwant 1.4.0"
# -------------------------------------------------
# The Graphene model 
# -------------------------------------------------
# Mario Solís
# PhD student in Physics
# Grupo de Partículas y Campos
# Instituto Balseiro - Centro Atómico Bariloche
# m.f.solis.benites@gmail.com
# ===========================================================================================================================
# Call math function comands
from math import pi, sqrt, tanh
#
# Call Kwant!
import kwant
#
# ===========================================================================================================================
# For computing eigenvalues IMPORTANT!!!
import scipy.sparse.linalg as sla
# ===========================================================================================================================
#
# Import the Plot lybrary of matplotlib, obviously for plotting!
#
from matplotlib import pyplot
#
# ===========================================================================================================================
# Now we start define the primitive vectors of the lattice, to this end, we use "kwant.lattice.general" comand
#
sin_30, cos_30 = (1 / 2, sqrt(3) / 2)
graphene = kwant.lattice.general([(1, 0), (sin_30, cos_30)],
                                 [(0, 0), (0, 1 / sqrt(3))])
a, b = graphene.sublattices # The honeycomb lattice has two basis atoms: a and b.
#
# ===========================================================================================================================
# Now we define the scattering region as circular shape!
# REMIND: "w" is the width and "pot" is the maximun potential of the p-n junction
#
def make_system(r=10, w=2.0, pot=0.1):

    def circle(pos):
        x, y = pos
        return x ** 2 + y ** 2 < r ** 2 # clearly is just a circle without boundary
    syst = kwant.Builder()

    def potential(site):
        (x, y) = site.pos
        d = y * cos_30 + x * sin_30
        return pot * tanh(d/ w)

    syst[graphene.shape(circle, (0, 0))] = potential
#
# ===========================================================================================================================
# Then, is time for the hoppings!
#
    hoppings = (((0, 0), a, b), ((0, 1), a, b), ((-1, 1), a, b))# Remind, the directions are not more orthogonal anymore, this due by the primitive vectors.
# 
# Thus, we add the hoppings using "kwant.builder.HoppingKind",
    syst[[kwant.builder.HoppingKind(*hopping) for hopping in hoppings]] = -1
#
# ===========================================================================================================================
# IMPORTANT!!! We remove one atom of the lattice A (that also remove the hoppings from/to this site).
#
    del syst[a(0, 0)]
    syst[a(-2, 1), b(2, 2)] = -1
#
# ===========================================================================================================================
# Is time for the Leads!
#
    sym0 = kwant.TranslationalSymmetry(graphene.vec((-1, 0)))

    def lead0_shape(pos):
        x, y = pos
        return (-0.4 * r < y < 0.4 * r)

    lead0 = kwant.Builder(sym0)
    lead0[graphene.shape(lead0_shape, (0, 0))] = - pot
    lead0[[kwant.builder.HoppingKind(*hopping) for hopping in hoppings]] = -1

    sym1 = kwant.TranslationalSymmetry(graphene.vec((0, 1)))

    def lead1_shape(pos):
        v = pos[1] * sin_30 - pos[0] * cos_30
        return (-0.4 * r < v < 0.4 * r)

    lead1 = kwant.Builder(sym1)
    lead1[graphene.shape(lead1_shape, (0, 0))] = pot
    lead1[[kwant.builder.HoppingKind(*hopping) for hopping in hoppings]] = -1

    return syst, [lead0, lead1]
#
# ===========================================================================================================================
# Compute the eigenvalues of the  CLOSED system
#
def compute_evs(syst):
    # Compute some eigenvalues of the closed system
    sparse_mat = syst.hamiltonian_submatrix(sparse=True)

    evs = sla.eigs(sparse_mat, 2)[0]
    print(evs.real)


def plot_conductance(syst, energies):
    # Compute transmission as a function of energy
    data = []
    for energy in energies:
        smatrix = kwant.smatrix(syst, energy)
        data.append(smatrix.transmission(0, 1))

    pyplot.figure()
    pyplot.plot(energies, data)
    pyplot.xlabel("Energy [t]")
    pyplot.ylabel("Conductance [e^2/h]")
    pyplot.show()


def plot_bandstructure(flead, momenta):
    bands = kwant.physics.Bands(flead)
    energies = [bands(k) for k in momenta]

    pyplot.figure()
    pyplot.plot(momenta, energies)
    pyplot.xlabel("Momentum [(lattice constant)^-1]")
    pyplot.ylabel("Energy [t]")
    pyplot.show()
#
# ===========================================================================================================================
# Define "main"  function
#
def main():
    pot = 0.1
    syst, leads = make_system(pot=pot)

    # To highlight the two sublattices of graphene, we plot one with
    # a filled, and the other one with an open circle:
    def family_colors(site):
        return 0 if site.family == a else 1

    # Plot the closed system without leads.
    kwant.plot(syst, site_color=family_colors, site_lw=0.1, colorbar=False)

    # Compute some eigenvalues.
    compute_evs(syst.finalized())

    # Attach the leads to the system.
    for lead in leads:
        syst.attach_lead(lead)

    # Then, plot the system with leads.
    kwant.plot(syst, site_color=family_colors, site_lw=0.1,
               lead_site_lw=0, colorbar=False)

    # Finalize the system.
    syst = syst.finalized()

    # Compute the band structure of lead 0.
    momenta = [-pi + 0.02 * pi * i for i in range(101)]
    plot_bandstructure(syst.leads[0], momenta)

    # Plot conductance.
    energies = [-2 * pot + 4. / 50. * pot * i for i in range(51)]
    plot_conductance(syst, energies)


# Call the main function if the script gets executed (as opposed to imported).
if __name__ == '__main__':
    main()

