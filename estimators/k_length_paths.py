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


# The number of potentially overlapping k-hop spanners is at most the diameter
# of the graph, divided by k (need to show this concept visually in the paper).
def estimate_number_k_hop_spanners(k, diameter):
    return float(diameter) / float(k)


def estimate_size_k_hop_spanner(k, degrees, diameter):
    number_k_hop_views = estimate_number_k_hop_spanners(k, diameter)

    # This one is already a loose upper bound:
    number_k_length_paths = estimate_number_k_length_paths(k, degrees)

    # The total number of k-length paths is how many edges you have,
    # so we can estimate the size of an "average" k-hop spanner in the
    # graph by dividing total number of edges by how many views you
    # can have.

    return float(number_k_length_paths) / float(number_k_hop_views)


# NOTE: These are actually combining multiple k-hop spanner views,
# since the k-length paths are overlapping (and not "stitched
# together", as we have them in the paper).
def estimate_number_k_length_paths(k, degrees):
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

    print 'Computing diameter...'
    # FIXME: Computing the diameter of a large graph is quite slow.
    # There are approximation algorithms which use sketches like
    # HyperLogLog to approximate the neighborhood function of a node.
    # They then use that as input to estimate path lengths in the graph,
    # which allows you to approximate the diameter. Example algos include
    # HyperANF, ANF (implemented in SNAP), HADI and others. Good overview
    # at Tej Chahed's MS thesis:
    # http://dprg.cs.uiuc.edu/docs/tej_thesis/graph-sampling-2.pdf
    #
    # For now we just assume a large diameter.
    diameter = 100  # nx.diameter(G)
    print 'Done computing diameter!'

    k_length_paths = estimate_number_k_length_paths(k, degrees)
    spanner_size = estimate_size_k_hop_spanner(k, degrees, diameter)

    print 'estimated number of overlapping k-length paths=%s' % k_length_paths
    print 'estimated average size of a k-hop spanner=%s' % spanner_size


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Unbiased estimator of number of k-length paths.')
    parser.add_argument('-k', help='Length of paths of which we want'
                        ' to estimate the likelihood.',
                        required=True,
                        type=int)
    parser.add_argument('-i', '--input_snap_graph',
                        help='Relative path of input graph in SNAP format.',
                        required=True)
    parser.add_argument('--is_directed',
                        help='Treats input graph as directed (default False).',
                        dest='is_directed',
                        action='store_true')
    parser.set_defaults(is_directed=False)
    args = parser.parse_args()

    main(**vars(args))
