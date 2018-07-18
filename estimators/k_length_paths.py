import argparse
import json
import networkx as nx
import numpy as np
import scipy
from scipy.special import binom

# (n choose (k+1)) * (m/(n choose 2))^k
def lower_bound(n, m, k):
  num_possible_k_hop_paths = binom(n, k+1)
  num_possible_edges = binom(n, 2)
  num_1_hop_paths = float(m) / float(num_possible_edges)

  return int(num_possible_k_hop_paths * pow(num_1_hop_paths, k))


# Estimates number of k-length paths.
def k_length_paths(n, m, k_hop_schema, config):
    k = len(k_hop_schema)
    k_length_paths = {}  # One for each value of alpha.

    for vertex_type in k_hop_schema:
      outdegrees_file = config[vertex_type]
      print 'loading %s degrees...' % vertex_type
      data = np.loadtxt(outdegrees_file, delimiter=' ')
      degrees = data[:,1]

      # Number of nodes
      n = len(degrees)

      for alpha in [50,95,'avg']:
          if alpha not in k_length_paths:
              k_length_paths[alpha] = 0

          if alpha == 'avg':
              deg_summary = np.average(degrees)
              factor = 0.5
          else:
              deg_summary = np.percentile(degrees, alpha)
              factor = float(alpha)/100.0

          if alpha not in k_length_paths:
              k_length_paths[alpha] = 0

          k_length_paths[alpha] += int(float(n) * factor
              * pow(deg_summary, k))

    # Also estimate using lower bound formula.
#    k_length_paths['lower_bound'] = lower_bound(n, m, k)

    return k_length_paths


def main(num_nodes, num_edges, k_hop_schema, input_outdegree_files):
    print k_hop_schema
    config=json.loads(input_outdegree_files)

    print 'Using het. network estimator=%s' % k_length_paths(
        num_nodes, num_edges, k_hop_schema, config)


# e.g.,
# python k_length_paths.py -m 10000 -s author article author
# -i '{"author": \
#   "/big_fast_drive/jfon/nets/1graphdblp/outdeg_author-10k.txt", \
#   "article": \
#   "/big_fast_drive/jfon/nets/1graphdblp/outdeg_article-10k.txt"}'
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Estimate number of k-length paths in heterogeneous net.')
    parser.add_argument('-m', '--num_edges', type=int, required=True)
    parser.add_argument('-n', '--num_nodes', type=int, required=True)
    parser.add_argument('-s', '--k_hop_schema',
        help='Space separated list, e.g., author article author',
        nargs='+',
        required=True)
    parser.add_argument('-i', '--input_outdegree_files',
        help='One for each node type in a JSON object.',
        required=True)
    args = parser.parse_args()

    main(**vars(args))
