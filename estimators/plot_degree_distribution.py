import argparse
import matplotlib #; matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import seaborn as sns

# Seaborn global settings.
sns.set(context='paper',
        font='serif',
        font_scale=2,
        rc={
            'lines.linewidth': 3,
            'lines.markeredgewidth': 1,
            'patch.linewidth': 2  # bar border
        })
sns.set_style('white', {'axes.linewidth': 1.5,
                        'axes.edgecolor': '.0',
                        'axes.labelcolor': '.0',
                        'font.family': 'serif',
                        'font.serif': 'Times New Roman',
                        'grid.color': '.0',
                        'text.color': '.0',
                        'xtick.color': '.0',
                        'ytick.color': '.0'})
sns.despine(top=False, bottom=False)

FIG_HEIGHT = 4
FIG_WIDTH = 6.5


def plot_cdfs(data, output_suffix):
    # CDF
    output_filename = 'cdf-' + output_suffix
    s = float(data.sum())
    cdf = data.cumsum(0)/s
    plt.figure()
    plt.plot(range(len(cdf)),cdf,'bo',markeredgecolor='k')
    plt.xscale('log')
    plt.ylim([0,1])
    plt.ylabel('CDF')
    plt.xlabel('Degree')
    fig = plt.gcf()
    fig.set_figheight(FIG_HEIGHT)
    fig.set_figwidth(FIG_WIDTH)
    fig.savefig(output_filename, format='pdf', bbox_inches='tight')
    plt.clf()
    print 'CDF saved as %s' % output_filename

    # CCDF
    output_filename = 'ccdf-' + output_suffix
    ccdf = 1-cdf
    print 'ccdf: %s' % ccdf
    freq = ccdf * s
    print 'freq: %s' % freq
    plt.figure()
    plt.plot(range(len(freq)),freq,'bo',markeredgecolor='k')

    # Fit exponential for power-law (linear on log-log CCDF)
    # starting at 2nd rank of X log.
    log_x = np.log10(np.array(range(len(freq))) + 1)
    log_y = np.log10(np.array(freq) + 1)
    first_rank = 2
    log_x = log_x[first_rank:]
    log_y = log_y[first_rank:]

    # Don't try to fit if we don't have at least 2 points
    if len(log_x) > 2:
        m, b = np.polyfit(log_x, log_y, 1)  # log(y) = m*log(x) + b
        y_fit = np.power(10, m * log_x + b)
        plt.plot(range(first_rank, len(log_x) + first_rank), y_fit, ':')

    plt.xscale('log')
    plt.yscale('log')
#    plt.ylim([0,1])
    plt.ylabel('Freq. sample > x')
    plt.xlabel('Degree')
    fig = plt.gcf()
    fig.set_figheight(FIG_HEIGHT)
    fig.set_figwidth(FIG_WIDTH)
    fig.savefig(output_filename, format='pdf', bbox_inches='tight')
    plt.clf()
    print 'CCDF (frequency) saved as %s' % output_filename


def plot_loglog(data, output_suffix):
    output_filename = 'loglog-' + output_suffix
    plt.figure()
    plt.plot(range(len(data)), data, 'bo', markeredgecolor='k')
    plt.yscale('log')
    plt.xscale('log')
    plt.ylabel('Frequency')
    plt.xlabel('Degree')
    fig = plt.gcf()
    fig.set_figheight(FIG_HEIGHT)
    fig.set_figwidth(FIG_WIDTH)
    fig.savefig(output_filename, format='pdf', bbox_inches='tight')
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
