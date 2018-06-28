import argparse
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np


def plot_loglog(data, output_filename):
    plt.plot(range(len(data)), data, 'b*')
    plt.yscale('log')
    plt.xscale('log')
    plt.ylabel('Frequency')
    plt.xlabel('Degree')
    plt.savefig(output_filename)
    plt.clf()
    print 'Plot saved as %s.' % output_filename


def main(input_snap_graph, is_directed, output_filename):
    fh = open(input_snap_graph, 'rb')

    G = None
    if is_directed:
        G = nx.read_edgelist(fh, create_using=nx.DiGraph())
    else:
        G = nx.read_edgelist(fh)

    fh.close()

    degrees = [deg for (_, deg) in G.degree()]
    degree_hist = np.bincount(degrees)
    plot_loglog(degree_hist, output_filename)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Plots degree distribution in log-log normal with log binning.')
    parser.add_argument('-i', '--input_snap_graph',
                        help='Relative path of input graph in SNAP format.',
                        required=True)
    parser.add_argument('--is_directed',
                        help='Treats input graph as directed (default False).',
                        dest='is_directed',
                        action='store_true')
    parser.set_defaults(is_directed=False)
    parser.add_argument('-o', '--output_filename',
                        help='Relative path of output plot (png, pdf, or eps).',
                        required=True)
    args = parser.parse_args()

    main(**vars(args))
