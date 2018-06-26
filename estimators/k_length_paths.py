import argparse
import networkx as nx


def get_components(graph):
    ccs = None

    if nx.is_directed(graph):
        return nx.strongly_connected_components(graph)

    return nx.connected_components(graph)


def estimate_number_k_length_paths(k, graph):
    ccs = get_components(graph)

    degrees = []
    for cc in ccs:
        # Ignore components of size less than k, as they lack k-hop paths.
        if len(cc) < k:
            continue

        # Get degrees of vertices in remaining components.
        cc_degrees = [deg for (node, deg) in graph.degree() if node in cc]
        degrees += cc_degrees

    # Estimate number of k-length paths in the remaining vertices.
    n = len(degrees)
    deg_avg = sum(degrees) / (float(len(degrees)))
#    print 'n=%s, deg_avg=%s, k=%s' % (n, deg_avg, k)

    return int(n * pow(deg_avg, k))


def main(k, input_snap_graph, is_directed):
    fh = open(input_snap_graph, 'rb')

    G = None
    if is_directed:
        G = nx.read_edgelist(fh, create_using=nx.DiGraph())
    else:
        G = nx.read_edgelist(fh)

    fh.close()

    estimate = estimate_number_k_length_paths(k, G)
    print 'estimated number of k-length paths=%s' % estimate


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
