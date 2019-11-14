
import os
import pdb
import numpy as np
import argparse
import ag_unify_md as agum

AU_TO_EV = 27.211396


def find_gau_outs(maindir, method, basset):
    method = method.upper()
    basset = basset.upper()
    gaussian_output_files = []

    for subdir in os.listdir(maindir):
        subdir = maindir + subdir
        #pdb.set_trace()

        if os.path.isdir(subdir):
            subfiles = os.listdir(subdir)
            for subfile in subfiles:
                subfile = subdir + "/" + subfile

                if subfile.endswith(".gau.out") and method in subfile.upper() and basset in subfile.upper():
                    print(subfile)
                    gaussian_output_files.append(subfile)
                    break

    return gaussian_output_files


def read_gauout_energy(gau_out_file):
    with open(gau_out_file) as opened_gau_out_file:
        line = opened_gau_out_file.readline()

        while line != "":
            if "Counterpoise corrected energy" in line:
                # energy in eV
                counterpoise_corrected_energy = float(line.split()[4])
                break

            line = opened_gau_out_file.readline()

    return counterpoise_corrected_energy


def calculate_distance(gau_out_file, idxs_atm1, idxs_atm2):
    dimer_sys = agum.Unification()
    dimer_sys.read_gau_log(gau_out_file, read_summary=True)
    #pdb.set_trace()

    if len(idxs_atm1) == 1:
        cog1 = dimer_sys.ts_coords[-1][idxs_atm1[0]]
    else:
        cog1 = dimer_sys.get_cog(-1, *idxs_atm1)

    if len(idxs_atm2) == 1:
        cog2 = dimer_sys.ts_coords[-1][idxs_atm2[0]]
    else:
        cog2 = dimer_sys.get_cog(-1, *idxs_atm2)

    distance = np.linalg.norm(cog1 - cog2)
    return distance


def get_results(maindir, method, basset, index_atm1, index_atm2):
    output_files = find_gau_outs(maindir, method, basset)
    dists_and_energies = []
    #pdb.set_trace()

    for output_file in output_files:
        try:
            cdist = calculate_distance(output_file, index_atm1, index_atm2)
            cenergy = read_gauout_energy(output_file)
            dists_and_energies.append([cdist, cenergy])
        except IndexError:
            print(("Gaussian calculation has not finished yet or properly; filename is {}".format(output_file)))

    return dists_and_energies


def norm_results(results):
    # sort by dist
    results.sort(key=lambda x: float(x[0]))
    lowest_energy = sorted([i[1] for i in results])[0]
    #pdb.set_trace()
    results = [[i[0], (i[1] - lowest_energy) * AU_TO_EV] for i in results]
    return results


def write_results(results, filename="default_energies.txt"):
    output_str = "{:> 4.4f}{:> 12.10f}\n"
    with open(filename, "w") as opened_gau_out_file:
        opened_gau_out_file.write("Energies in eV\n")
        for dist_energy in results:
            opened_gau_out_file.write(output_str.format(dist_energy[0], dist_energy[1]))


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser()
    PARSER.add_argument("maindir",
                        help="Directory with all gaussian output files")
    PARSER.add_argument("method",
                        help="Name of the method, e.g. mp2, wb97xd which occurs in the filename")
    PARSER.add_argument("basset",
                        help="Name of the basis set, e.g. pvtz, pvdz which occurs in the filename")
    PARSER.add_argument("-a1",
                        nargs="*",
                        type=int,
                        help="Name of the atoms to measure the distance from; if more than one atom index is given, the center of geometry will be calculated")
    PARSER.add_argument("-a2",
                        nargs="*",
                        type=int,
                        help="Name of the atoms to measure the distance to; if more than one atom index is given, the center of geometry will be calculated")

    ARGS = PARSER.parse_args()
    MY_RESULTS = get_results(ARGS.maindir, ARGS.method, ARGS.basset, ARGS.a1, ARGS.a2)
    MY_RESULTS = norm_results(MY_RESULTS)
    write_results(MY_RESULTS, "{}_{}_scan_results.txt".format(ARGS.method, ARGS.basset))
