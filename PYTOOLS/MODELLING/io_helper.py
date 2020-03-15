"""
Functions that often needed for little I/O things.

Since some functions come quite in handy, reusing them is just reasonable. This
is what this module tries to achieve.
"""


from collections import OrderedDict

# ==============================================================================#
# Helping functions
# ==============================================================================#


def get_scanned_geometry(gau_log):
    """
    Bla.

    Search for keyword scan in gaussian output file and return type of geometry
    (e.g. bond, angle, dihedral) that was scanned.
    """
    with open(gau_log) as opened_gau_log:
        line = opened_gau_log.readline()
        while line != "":
            if "Initial Parameters" in line:
                while "Scan" not in line:
                    line = opened_gau_log.readline()
                else:
                    split_line = line.split()
                    scanned_geometry = split_line[2]
                    geometry_value = float(split_line[3])
                    break

            line = opened_gau_log.readline()

    return (scanned_geometry, geometry_value)


def write_energies(filename):
    """
    Bla.

    Write header for energies file.
    """
    with open(filename, "w") as opened_filename:
        opened_filename.write(
            "{:>20s} {:>20s}\n".format("Angle [Degrees]", "Energy [eV]")
        )


def norm_energy(energy_file_in, energy_file_out):
    """
    Bla.

    Norm energies by smallest value and sort by angle.
    """
    keys = []
    original_values = []

    with open(energy_file_in) as opened_energy_file:
        # skip header
        opened_energy_file.readline()
        line = opened_energy_file.readline()

        while line != "":
            split_line = line.split()
            keys.append(float(split_line[0]))
            original_values.append(float(split_line[1]))
            line = opened_energy_file.readline()

    min_value = min(original_values)
    normed_values = []

    for val in original_values:
        val -= min_value
        normed_values.append(val)

    keys_and_values = dict(list(zip(keys, normed_values)))
    keys_and_values = OrderedDict(sorted(keys_and_values.items()))

    write_energies(energy_file_out)

    with open(energy_file_out, "a") as opened_energy_file:

        for key, value in keys_and_values.items():
            opened_energy_file.write("{:> 20.8f} {:> 20.8f}\n".format(key, value))


def read_data(data_file):
    """
    Bla.

    Read data from a file which has a header and two columns.
    """
    x_values = []
    y_values = []

    with open(data_file) as opened_data_file:
        # header
        header_line = opened_data_file.readline()
        # first row with data
        line = opened_data_file.readline()

        while line != "":
            split_line = [float(i) for i in line.split()]
            x_values.append(split_line[0])
            y_values.append(split_line[1])
            line = opened_data_file.readline()

    return (x_values, y_values)


def shift_data(data_points, value):
    """
    Bla.

    Shift data by value
    """
    for idx, cdata in enumerate(data_points):
        if cdata < 0:
            data_points[idx] = cdata + value
