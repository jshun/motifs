import argparse
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np


def plot_cdfs(data, output_suffix):
    # CDF
    output_filename = 'cdf-' + output_suffix
    s = float(data.sum())
    cdf = data.cumsum(0)/s
    plt.plot(range(len(cdf)),cdf,'bo')
    plt.xscale('log')
    plt.ylim([0,1])
    plt.ylabel('CDF')
    plt.xlabel('Degree')
    plt.savefig(output_filename)
    plt.clf()
    print 'CDF saved as %s' % output_filename

    # CCDF
    output_filename = 'ccdf-' + output_suffix
    ccdf = 1-cdf
    plt.plot(range(len(ccdf)),ccdf,'bo')
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim([0,1])
    plt.ylabel('CCDF')
    plt.xlabel('Degree')
    plt.savefig(output_filename)
    plt.clf()
    print 'CCDF saved as %s' % output_filename


def plot_loglog(data, output_suffix):
    output_filename = 'loglog-' + output_suffix
    plt.plot(range(len(data)), data, 'b*')
    plt.yscale('log')
    plt.xscale('log')
    plt.ylabel('Frequency')
    plt.xlabel('Degree')
    plt.savefig(output_filename)
    plt.clf()
    print 'Log-log saved as %s.' % output_filename


def main(input_snap_graph, is_directed, output_suffix):
    fh = open(input_snap_graph, 'rb')

    G = None
    if is_directed:
        G = nx.read_edgelist(fh, create_using=nx.DiGraph())
    else:
        G = nx.read_edgelist(fh)

    fh.close()

    # TODO: Maximum likehood fitting if we want to fit the exponent:
    # https://cs.brynmawr.edu/Courses/cs380/spring2013/section02/slides/10_ScaleFreeNetworks.pdf
    keys, degrees = zip(*G.degree())
    degree_dist = np.bincount(degrees)

    plot_loglog(degree_dist, output_suffix)
    plot_cdfs(degree_dist, output_suffix)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Plots degree distribution CDF, CCDF and log-log.')
    parser.add_argument('-i', '--input_snap_graph',
                        help='Relative path of input graph in SNAP format.',
                        required=True)
    parser.add_argument('--is_directed',
                        help='Treats input graph as directed (default False).',
                        dest='is_directed',
                        action='store_true')
    parser.set_defaults(is_directed=False)
    parser.add_argument('-o', '--output_suffix',
                        help='Suffix of output filename, saved in current dir.',
                        required=True)
    args = parser.parse_args()

    main(**vars(args))
