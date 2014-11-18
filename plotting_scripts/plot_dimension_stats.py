import sys
import parse

import numpy as np
import matplotlib.pyplot as plt

low_name, avg_name, high_name = sys.argv[1:4]
for (fname, line_label) in [(low_name, "Lowest Density"),
                            (avg_name, "Average Density"),
                            (high_name, "Highest Density")]:
    rs, r_freqs = parse.get_col(fname, 0, 1)
    vs, v_freqs = parse.get_col(fname, 2, 3)

    # for i in xrange(2, len(rs) + 1):
    #     r_freqs[-i] += r_freqs[-i + 1]
    #     v_freqs[-i] += v_freqs[-i + 1]

    # r_freqs /= max(r_freqs)
    # v_freqs /= max(v_freqs)

    plt.figure(1)
    plt.plot(np.log10(rs), r_freqs, linewidth=2, label=line_label)

    plt.figure(2)
    plt.plot(np.log10(vs), v_freqs, linewidth=2, label=line_label)

plt.figure(1)
plt.title("$L$ = 125, $z$ = 0")
plt.ylabel(r"$N(\log_{10} R) / N / \Delta \log_{10} R$")
plt.xlabel(r"$\log_{10}R$")
plt.legend(loc=3)

plt.figure(2)
plt.title("$L$ = 125, $z$ = 0")
plt.ylabel(r"$N(\log_{10} V) / N / \Delta \log_{10} V$")
plt.xlabel(r"$\log_{10}V$")

b_cell_v = np.log10(125.0**3 / (1<<30))
init_v = np.log10(125.0**3 / (1<<30) / 6)
min_y, max_y = plt.ylim()

plt.plot([b_cell_v, b_cell_v],  [min_y, max_y], "--k")
plt.plot([init_v, init_v],  [min_y, max_y], "--k")

plt.legend(loc=3)
plt.show()
