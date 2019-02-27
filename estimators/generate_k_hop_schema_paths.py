import itertools


# Same thing we do in Prolog, but now in Python.
def generate_k_hop_schema_paths(edges, paths_so_far, k, curr_k):
    if curr_k == 0:
        return [p for p in paths_so_far if len(p) == k]

    if k == curr_k:
        new_paths = [[e] for e in edges]
        return generate_k_hop_schema_paths(edges, new_paths, k, k - 1)

    new_paths = []
    for i, path in enumerate(paths_so_far):
        src = path[0][0]
        dst = path[-1][1]

        for j, edge in enumerate(edges):
            # Add edge to the end of the path.
            if dst == edge[0]:
                new_paths.append(path + [edge])

            # Add edge to the front of the path.
            if src == edge[1]:
                new_paths.append([edge] + path)
    paths_so_far = new_paths

    # Remove duplicates.
    paths_so_far.sort()
    paths_so_far = list(i for i, _ in itertools.groupby(paths_so_far))

    # Remove paths that did not grow any longer
    paths_so_far = [p for p in paths_so_far if len(p) == (k - curr_k + 1)]
    return generate_k_hop_schema_paths(edges, paths_so_far, k, curr_k - 1)


def to_dot_graph(view, view_name):
    filename = 'dot_graphs/schema_k_hop_%s.dot' % view_name
    with open(filename, 'w') as f:
        f.write('digraph %s {\n' % view_name)
        # left-to-right layout
        f.write(' rankdir="LR";\n')

        # all k-hop schema views are simple directed paths, so number of node is k+1
        k = len(view)
        num_nodes = k + 1
        node_idx = 0
        for i in xrange(k):
            edge = view[i]
            f.write(' %d [label="%s"];\n' % (node_idx + 1, edge[0]))

            # print dst node if this is the last edge in the path
            if i == k - 1:
                f.write(' %d [label="%s"];\n' % (node_idx + 2, edge[1]))
            node_idx += 1

        f.write('\n')

        # second pass over the edges to print them out
        node_idx = 0
        for edge in view:
            f.write(' %d -> %d [label="%s"];\n' % (node_idx + 1, node_idx + 2,
                                                   edge[2]))
            node_idx += 1

        # close digraph
        f.write('}\n')


def count_views(schema, k_lower, k_upper):
    print '\n%s:' % schema['name']
    edges = schema['edges']

    for schema_edge in edges:
        print schema_edge

    nodes = list(itertools.chain(*edges))
    nodes = list(set(nodes))

    for k in xrange(k_lower, k_upper + 1):
        views = generate_k_hop_schema_paths(edges, [], k, k)
        N = schema['N']
        M = schema['M']  # same as len(schema['edges'])
        print 'k=%d, N=%d, M=%d, #views=%d' % (k, N, M, len(views))

        views_slice = views[:3]
        if len(views_slice) > 0:
            print 'saving first %d views as dot graph...' % len(views_slice)
            for i in xrange(len(views_slice)):
                view = views_slice[i]
                view_name = 'N%d_M%d_k%d_view%d' % (N, M, k, i)
                to_dot_graph(view, view_name)


# yapf: disable
schemas = [{
    'name': 'Homogeneous Single Layer (N=1,M=1)',
    'edges': [('VType_1', 'VType_1', 'EType_1')],
    'N': 1,
    'M': 1
}, {
    'name':
    'Homogeneous Multilayer (N=1,M=2)',
    'edges': [('VType_1', 'VType_1', 'EType_1'),
              ('VType_1', 'VType_1', 'EType_2')],
    'N': 1,
    'M': 2
}, {
    'name':
    'Homogeneous Multilayer (N=1,M=3)',
    'edges': [('VType_1', 'VType_1', 'EType_1'),
              ('VType_1', 'VType_1', 'EType_2'),
              ('VType_1', 'VType_1', 'EType_3')],
    'N': 1,
    'M': 3
}, {
    'name':
    'Homogeneous Multilayer (N=1,M=4)',
    'edges': [('VType_1', 'VType_1', 'EType_1'),
              ('VType_1', 'VType_1', 'EType_2'),
              ('VType_1', 'VType_1', 'EType_3'),
              ('VType_1', 'VType_1', 'EType_4')],
    'N': 1,
    'M': 4
}, {
    'name':
    'Heterogeneous Multilayer w/ Cycle (N=2,M=2)',
    'edges': [('VType_1', 'VType_2', 'EType_1'),
              ('VType_2', 'VType_1', 'EType_2')],
    'N': 2,
    'M': 2
}, {
    'name':
    'Heterogeneous Multilayer w/ Cycle (N=3,M=3)',
    'edges': [('VType_1', 'VType_2', 'EType_1'),
              ('VType_2', 'VType_3', 'EType_2'),
              ('VType_3', 'VType_1', 'EType_3')],
    'N': 3,
    'M': 3
}, {
    'name':
    'Heterogeneous Multilayer w/o Cycle (N=3,M=2)',
    'edges': [('VType_1', 'VType_2', 'EType_1'),
              ('VType_2', 'VType_3', 'EType_2')],
    'N': 3,
    'M': 2
}, {
    'name':
    'Heterogeneous Multilayer w/o Cycle (N=4,M=3)',
    'edges': [('VType_1', 'VType_2', 'EType_1'),
              ('VType_2', 'VType_3', 'EType_2'),
              ('VType_3', 'VType_4', 'EType_3')],
    'N': 4,
    'M': 3
}]
# yapf: enable

for schema in schemas:
    count_views(schema, 2, 8)
