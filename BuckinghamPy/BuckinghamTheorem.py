import numpy
import BuckinghamPy.buckinghampy.buckinghampi as bpy
import BuckinghamPy.buckinghampy.buckinghampigui as bpygui

def main():
    # initialize BuckinghamPy object
    PBE = bpy.BuckinghamPi()
    # add all variables with units to the system
    # volume scale
    PBE.add_variable(name="M1", dimensions="L^3*P")
    # particle unit scale
    PBE.add_variable(name="M0", dimensions="P")
    # aggregation rate
    PBE.add_variable(name="a", dimensions="T^(-1)*P^(-1)")
    # growth rate
    PBE.add_variable(name="G", dimensions="L^(3)*T^(-1)")
    # breakage rate
    PBE.add_variable(name="beta", dimensions="T^(-1)")
    # time
    PBE.add_variable(name="t", dimensions="T")

    # generate dimensionless pi terms and delete all repeating sets for all combinations of the dimensional matrix R
    PBE.generate_pi_terms()

    # print the remaining pi terms in Jupyter cell and/or console format
    # Jupyter cell
    # PBE.print_all()
    # console output (preferred as it is more accessible)
    print(PBE.pi_terms)

if __name__ == '__main__':
    main()  