#!/usr/bin/python
import sys
from itertools import izip

def main():
    if len(sys.argv) < 2 or len(sys.argv) == 2 and sys.argv[1] == "-h":
        print "Usage: ./plot_output.py filename [surf | contour]"
        return
    filename = sys.argv[1]
    drawing_type = "surf"
    if len(sys.argv) > 2:
        drawing_type = sys.argv[2]
        if drawing_type not in ["surf", "contour"]:
            print "Usage: ./plot_output.py filename [surf | contour]"
            return

    alphas = set()
    betas = set()
    errors = {}
    with open(filename, "r") as file:
        for line in file:
            alpha, beta, error = [float(value) for value in line.strip().split()]
            errors.setdefault(alpha, {})
            errors[alpha][beta] = error

            alphas.add(alpha)
            betas.add(beta)

    import matplotlib.pyplot as plt
    from matplotlib import cm
    from mpl_toolkits.mplot3d import Axes3D
    import numpy as np

    alphas = sorted(list(alphas))
    betas = sorted(list(betas))

    alpha_indexes = dict(izip(alphas, xrange(len(alphas))))
    beta_indexes = dict(izip(betas, xrange(len(betas))))

    alphas, betas = np.meshgrid(alphas, betas)

    plot_errors = np.zeros((len(alphas), len(betas)))
    for alpha, betas_dict in errors.iteritems():
        alpha_index = alpha_indexes[alpha]
        for beta, error in betas_dict.iteritems():
            beta_index = beta_indexes[beta]
            plot_errors[alpha_index][beta_index] = error * 2

    print alphas.shape
    print betas.shape
    print plot_errors.shape

    fig = plt.figure()

    if drawing_type == "surf":
        ax = fig.gca(projection='3d')
        surf = ax.plot_wireframe(betas, alphas, plot_errors,
            rstride=10, cstride=10, cmap=cm.coolwarm,
            linewidth=1, antialiased=False,
            colors="k")
        plt.xlabel("beta")
        plt.ylabel("alpha")
    else:
        CS = plt.contour(betas, alphas, plot_errors, range(0, 300, 30), colors='k')
        plt.xlabel("beta")
        plt.ylabel("alpha")
#        plt.clabel(CS, fontsize=9, inline=1)

    plt.show()

if __name__ == "__main__":
    main()

