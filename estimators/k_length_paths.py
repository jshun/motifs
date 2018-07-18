import argparse
import json
import networkx as nx
import numpy as np


# Estimates number of k-length paths.
def k_length_paths(k_hop_schema, config):
    k = len(k_hop_schema)
    k_length_paths = 1

    for vertex_type in k_hop_schema:
      outdegrees_file = config[vertex_type]
      print 'loading degrees...'
      data = np.loadtxt(outdegrees_file, delimiter=' ')
      degrees = data[:,1]

      # Number of nodes
      n = len(degrees)

      # Using average.
      #    deg_summary = sum(degrees) / (float(len(degrees)))
      # Using median.
      #    deg_summary = np.percentile(degrees, 50)
      # Using 95th %ile.
      deg_summary = np.percentile(degrees, 50)

      k_length_paths *= n * pow(deg_summary, k)

    return int(k_length_paths)


def main(k_hop_schema, input_outdegree_files):
    print k_hop_schema
    config=json.loads(input_outdegree_files)

    print 'Using het. network estimator=%s' % k_length_paths(
        k_hop_schema, config)


# e.g.,
# python k_length_paths.py -s author \
#   -i '{"author": \
#   "/big_fast_drive/jfon/nets/1graphdblp/outdeg_author.txt", \
#   "article": "/big_fast_drive/jfon/nets/1graphdblp/outdeg_author.txt"}'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Estimate number of k-length paths in heterogeneous net.')
    parser.add_argument('-s', '--k_hop_schema',
        help='Space separated list, e.g., author article author',
        nargs='+',
        required=True)
    parser.add_argument('-i', '--input_outdegree_files',
        help='One for each node type in a JSON object.',
        required=True)
    args = parser.parse_args()

    main(**vars(args))
