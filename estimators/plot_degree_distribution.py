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

FIG_HEIGHT = 4  # 6.5
FIG_WIDTH = 4


# Round up to nearest power of 10.
def round_up(x):
  return int(10**(np.ceil(np.log10(x))))


def plot_cdfs(data, title):
    # CDF
    output_filename = 'cdf-%s.pdf' % title
    s = float(data.sum())
    cdf = data.cumsum(0)/s
    plt.figure()
    x = range(len(cdf))
    plt.plot(x, cdf, 'bo',markeredgecolor='k')
    max_x = round_up(np.max(x))
    plt.xlim([-0.01, max_x])
    plt.ylim([-0.01, 1.01])  # CDF.
    plt.xscale('symlog')
    ax = plt.gca()
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    plt.title(title)
    plt.ylabel('CDF')
    plt.xlabel('Degree')
    fig = plt.gcf()
    fig.set_figheight(FIG_HEIGHT)
    fig.set_figwidth(FIG_WIDTH)
    fig.savefig(output_filename, format='pdf', bbox_inches='tight')
    plt.clf()
    print 'CDF saved as %s' % output_filename

    # CCDF
    output_filename = 'ccdf-%s.pdf' % title
    ccdf = 1-cdf
    freq = ccdf * s
    plt.figure()
    plt.plot(range(len(freq)),freq,'bo',markeredgecolor='k')

    # Fit exponential for power-law (linear on log-log CCDF)
    # starting at 2nd rank of X log.
    x = np.array(range(len(freq))) + 1
    y = np.array(freq) + 1
    log_x = np.log10(x)
    log_y = np.log10(y)
    first_rank = 2
    log_x = log_x[first_rank:]
    log_y = log_y[first_rank:]
    # Don't try to fit if we don't have at least 2 points
    if len(log_x) > 2:
        m, b = np.polyfit(log_x, log_y, 1)  # log(y) = m*log(x) + b
        y_fit = np.power(10, m * log_x + b)
        plt.plot(range(first_rank, len(log_x) + first_rank),
            y_fit, '--', color='red', linewidth=3.0, alpha=1.0)
    max_y = round_up(np.max(y))
    max_x = round_up(np.max(x))
    plt.title(title)
    plt.xlim([-0.01, max_x])
    plt.ylim([-0.01, max_y])
    plt.yscale('symlog')
    plt.xscale('symlog')
    ax = plt.gca()
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    plt.ylabel('Freq. deg. > x')
    plt.xlabel('Degree')
    fig = plt.gcf()
    fig.set_figheight(FIG_HEIGHT)
    fig.set_figwidth(FIG_WIDTH)
    fig.savefig(output_filename, format='pdf', bbox_inches='tight')
    plt.clf()
    print 'CCDF (frequency) saved as %s' % output_filename


def plot_loglog(data, title):
    output_filename = 'loglog-%s.pdf'
    plt.figure()
    x = range(len(data))
    y = data
    plt.plot(x, y, 'bo', markeredgecolor='k')
    max_y = round_up(np.max(y))
    max_x = round_up(np.max(x))
    plt.title(title)
    plt.xlim([-0.01, max_x])
    plt.ylim([-0.01, max_y])
    plt.yscale('symlog')
    plt.xscale('symlog')
    ax = plt.gca()
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    plt.ylabel('Frequency')
    plt.xlabel('Degree')
    fig = plt.gcf()
    fig.set_figheight(FIG_HEIGHT)
    fig.set_figwidth(FIG_WIDTH)
    fig.savefig(output_filename, format='pdf', bbox_inches='tight')
    plt.clf()
    print 'Log-log saved as %s.' % output_filename


def main(input_file, file_type, is_directed, title):
    if file_type == 'SNAP':
      fh = open(input_file, 'rb')

      G = None
      if is_directed:
          G = nx.read_edgelist(fh, create_using=nx.DiGraph())
      else:
          G = nx.read_edgelist(fh)

      fh.close()
      keys, degrees = zip(*G.degree())
    else:
      data = np.loadtxt(input_file, delimiter=' ', dtype='i8')
      degrees = data[:,1]

    degree_dist = np.bincount(degrees)

#    plot_loglog(degree_dist, title)
    plot_cdfs(degree_dist, title)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Plots degree distribution CDF, CCDF and log-log for a given input SNAP graph, or input degrees file.')
    parser.add_argument('-i', '--input_file',
          help='Input file, which is either a graph in SNAP format (default) or a degrees file',
          required=True)
    parser.add_argument('--is_directed',
          help='Treats graph as directed (default False).',
          dest='is_directed',
          action='store_true')
    parser.set_defaults(is_directed=False)
    parser.add_argument('--title',
          help='Plot title and also used as output file suffix.',
          required=True)
    parser.add_argument('--file_type',
          const='SNAP',
          nargs='?',
          choices=['SNAP', 'degrees'],
          help='Whether input file is either SNAP or a space separated list of "nodeid outdegree" entries',
          required=True)

    args = parser.parse_args()

    main(**vars(args))
