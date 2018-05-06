from __future__ import print_function, division
import matplotlib.pyplot as plt
import numpy as np

__version__ = "2017-08-11"


class Paramfit(object):
    """
    Class to automatize paramfit fitting procedure.
    """
    def __init__(self, nstructures):
        self.nstructures = nstructures

    def write_set_params_in(self, set_params_in, set_params_out="set_params.out"):
        """
        Write job control file to set the parameters
        """
        with open(set_params_in, "w") as f_out:
            f_out.write(
                "RUNTYPE=SET_PARAMS\n" +
                "PARAMETER_FILE_NAME={}\n".format(set_params_out)
            )

    def write_fit_k(self, fit_k_in):
        """
        Write job control file to fit K (initial fit)
        """
        with open(fit_k_in, "w") as f_out:
            f_out.write(
                "RUNTYPE=FIT\n"
                "PARAMETERS_TO_FIT=K_ONLY\n" +
                "NSTRUCTURES={}\n".format(self.nstructures) +
                "COORDINATE_FORMAT=TRAJECTORY\n" +
                "FUNC_TO_FIT=SUM_SQUARES_AMBER_STANDARD\n" +
                "QM_ENERGY_UNITS=HARTREE\n"
            )

    def read_fit_k(self, fit_k_out):
        # extract K from fit_K.out
        with open(fit_k_out, "r") as f_in:
            for line in f_in:

                if "FINAL PARAMETERS" in line:
                    # skip the following two lines
                    f_in.next()
                    f_in.next()
                    force_const_K = float(f_in.next().split()[2])
                    print("***Paramfit-Info: Found value of K: ", force_const_K)
                    break
        return force_const_K

    def write_fit_params(self,
                         fit_params_in,
                         set_params_out,
                         force_const_K,
                         optimizations=200,
                         max_generations=20000,
                         generations_to_conv=20,
                         search_space=-1.000000,
                         file_fitted_energies="fit_output_energy.dat",
                         fitted_frcmod="fitted_params.frcmod"):
        with open(fit_params_in, "w") as f_out:
            f_out.write(
                "RUNTYPE=FIT\n" +
                "PARAMETERS_TO_FIT=LOAD\n" +
                "PARAMETER_FILE_NAME={}\n".format(set_params_out) +
                "COORDINATE_FORMAT=TRAJECTORY\n" +
                "NSTRUCTURES={}\n".format(self.nstructures) +
                "K={}\n".format(force_const_K) +
                "FUNC_TO_FIT=SUM_SQUARES_AMBER_STANDARD\n" +
                "QM_ENERGY_UNITS=HARTREE\n" +
                "ALGORITHM=BOTH\n" +
                "OPTIMIZATIONS={}\n".format(optimizations) +
                "MAX_GENERATIONS={}\n".format(max_generations) +
                "GENERATIONS_TO_CONV={}\n".format(generations_to_conv) +
                "GENERATIONS_TO_SIMPLEX=5\n" +
                "GENERATIONS_WITHOUT_SIMPLEX=5\n" +
                "MUTATION_RATE=0.100000\n" +
                "PARENT_PERCENT=0.250000\n" +
                "SEARCH_SPACE={}\n".format(search_space) +
                "SORT_MDCRDS=OFF\n" +
                "WRITE_ENERGY={}\n".format(file_fitted_energies) +
                "WRITE_FRCMOD={}\n".format(fitted_frcmod)
            )

    def plot_energy(self, energy_dat):
        """
        Plot paramfit output file.
        """
        num = []
        amber_k = []
        quantum = []
        initial_amber_k = []

        with open(energy_dat) as energy_dat_in:
            energy_dat_in.next()
            energy_dat_in.next().split()
            for line in energy_dat_in:
                line = line.split()
                num.append(int(line[0]))
                amber_k.append(float(line[1]))
                quantum.append(float(line[2]))
                initial_amber_k.append(float(line[3]))

        custom_fontsize = 10
        mymarkersize = 3

        xyfig = plt.figure()
        xyfig.canvas.set_window_title("Fit energy")
        plt.title("Fit energy")
        plt.xlabel("Structure", fontsize=custom_fontsize)
        plt.ylabel("Energy / kcal*mol-1", fontsize=custom_fontsize)
        plt.xticks(np.arange(min(num), max(num)+1, 5))
        plt.plot(num, amber_k, "r-^", markersize=mymarkersize, antialiased=True, label="Fit Amber")
        plt.plot(num, initial_amber_k, "b-*", markersize=mymarkersize, antialiased=True, label="Initial Amber")
        plt.plot(num, quantum, "k->", markersize=mymarkersize, antialiased=True, label="Quantum")
        plt.legend(bbox_to_anchor=(0., 1., 1., .1), loc="upper right",
                   borderaxespad=0., frameon=True, shadow=False, numpoints=1,
                   prop={'size': 8})
        plt.show()

    def write_tleap_in(self,
                       mol2,
                       filename="tleap.in",
                       forcefield="leaprc.gaff2",
                       prmtop_out="refined.prmtop",
                       inpcrd_out="refined.inpcrd",
                       *frcmods):
        """
        Write input file for tleap.
        """
        with open(filename, "w") as f_out:
            f_out.write("source {}\n".format(forcefield) +
                        "MOL = loadmol2 {}\n".format(mol2))

            for cfrcmod in frcmods:
                f_out.write("loadamberparams {}\n".format(cfrcmod))

            f_out.write("check MOL\n" +
                        "saveamberparm MOL {} {}\n".format(prmtop_out, inpcrd_out) +
                        "quit\n")
