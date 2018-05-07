from __future__ import print_function
import sys
import collections
import matplotlib.pyplot as plt

script, energy_file = sys.argv

energies = {}

with open(energy_file) as en_in:
    for line in en_in:
        line = line.split()
        line[0] = int(float(line[0]))
        line[1] = float(line[1])
        energies[line[0]] = line[1]

# order by angles
od = collections.OrderedDict(sorted(energies.items()))

xyfig = plt.figure()

xvals = od.keys()
yvals = od.values()

plt.plot(xvals, yvals, "r-", label="%r vs. %r" % ("Angle", "Energy"),
         antialiased=True)

ax = plt.gca()
ax.ticklabel_format(useOffset=False, style='plain')
# Set x- and y-labels
plt.xlabel("Angle", fontsize=12, fontweight='bold')
plt.ylabel("Energy", fontsize=12, fontweight='bold')
# Add legend outside
plt.legend(bbox_to_anchor=(0., 1., 1., .1),
           loc="upper right", borderaxespad=0., frameon=True,
           shadow=True, numpoints=1, prop={'size': 10})
plt.title("Scatter-Plot: %r vs. %r" % ("Angle", "Energy"), fontweight='bold')
plt.show()
