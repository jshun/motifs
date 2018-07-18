import argparse
import networkx as nx
import numpy as np


def get_components(graph):
    ccs = None

    if nx.is_directed(graph):
        return nx.strongly_connected_components(graph)

    return nx.connected_components(graph)


# Heuristic: we use only the degrees of nodes that fall in connected
# components which are of size at least k.
def get_degrees(k, graph):
    print 'Computing connected components...'
    ccs = get_components(graph)
    print 'Done computing CCs!'

    degrees = []
    for cc in ccs:
        # Ignore components of size less than k, as they lack k-hop paths.
        if len(cc) < k:
            continue

        # Get degrees of vertices in remaining components.
        degrees += [deg for (node, deg) in graph.degree() if node in cc]

    return degrees


# Estimates number of k-length paths using CC heuristic.
def k_length_paths_CC(k, degrees):
    if len(degrees) == 0:
      return 0

    # Using average.
    #    deg_summary = sum(degrees) / (float(len(degrees)))
    # Using median.
    #    deg_summary = np.percentile(degrees, 50)
    # Using 95th %ile.
    deg_summary = np.percentile(degrees, 95)

    # Number of nodes.
    n = len(degrees)

    return int(n * pow(deg_summary, k))


def main(k, input_snap_graph, is_directed):
    fh = open(input_snap_graph, 'rb')

    G = None
    if is_directed:
        G = nx.read_edgelist(fh, create_using=nx.DiGraph())
    else:
        G = nx.read_edgelist(fh)

    fh.close()

    degrees = get_degrees(k, G)
    k_length_paths = k_length_paths_CC(k, degrees)

    print 'Using CC heuristic=%s' % k_length_paths


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Estimator of number of k-length paths.')
    parser.add_argument('-k', help='Path length',
                        required=True,
                        type=int)
    parser.add_argument('-i', '--input_snap_graph',
                        help='Input graph in SNAP format.')
    parser.add_argument('--is_directed',
                        help='Treats input graph as directed (default False).',
                        dest='is_directed',
                        action='store_true')
    parser.set_defaults(is_directed=False)
    args = parser.parse_args()

    main(**vars(args))
