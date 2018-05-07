import os

filename = "/home/user1/SSHFS/zahnserver/Research/MODELLING/CBZ_Polymorphs/built_structures/CBZ_II/3.CBZII.xyz"
filename_out = os.path.basename(filename)

with open(filename) as f:
    num_atoms = int(f.next())
    f.next()

    for i in xrange(num_atoms/30):
        with open("{}_{}.xyz".format(filename_out, i), "w") as fout:
            fout.write("30")
            fout.write("\n\n")

            for j in xrange(30):
                line = f.next()
                fout.write(line)

