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


def count_views(schema, k_lower, k_upper):
    print '\n%s:' % schema['name']
    edges = schema['edges']

    for schema_edge in edges:
        print schema_edge

    nodes = list(itertools.chain(*edges))
    nodes = list(set(nodes))

    for k in xrange(k_lower, k_upper + 1):
        views = generate_k_hop_schema_paths(edges, [], k, k)
        print 'k=%d, N=%d, M=%d, #views=%d' % (k, len(nodes), len(edges),
                                               len(views))

        subset = views[:3]
        print 'first %d views:' % len(subset)
        for view in subset:
            print view


schemas = [{
    'name': 'Homogeneous Single Layer (N=1,M=1)',
    'edges': [('Node', 'Node', 'Edge')]
}, {
    'name':
    'Homogeneous Multilayer (N=1,M=2)',
    'edges': [('Node', 'Node', 'EdgeType1'), ('Node', 'Node', 'EdgeType2')]
}, {
    'name':
    'Homogeneous Multilayer (N=1,M=3)',
    'edges': [('Node', 'Node', 'EdgeType1'), ('Node', 'Node', 'EdgeType2'),
              ('Node', 'Node', 'EdgeType3')]
}, {
    'name':
    'Homogeneous Multilayer (N=1,M=4)',
    'edges': [('Node', 'Node', 'EdgeType1'), ('Node', 'Node', 'EdgeType2'),
              ('Node', 'Node', 'EdgeType3'), ('Node', 'Node', 'EdgeType4')]
}, {
    'name':
    'Heterogeneous Multilayer w/ Cycle (N=2,M=2)',
    'edges': [('NodeType1', 'NodeType2', 'EdgeType1'),
              ('NodeType2', 'NodeType1', 'EdgeType2')]
}, {
    'name':
    'Heterogeneous Multilayer w/ Cycle (N=3,M=3)',
    'edges': [('NodeType1', 'NodeType2',
               'EdgeType1'), ('NodeType2', 'NodeType3', 'EdgeType2'),
              ('NodeType3', 'NodeType1', 'EdgeType3')]
}, {
    'name':
    'Heterogeneous Multilayer w/o Cycle (N=3,M=2)',
    'edges': [('NodeType1', 'NodeType2', 'EdgeType1'),
              ('NodeType2', 'NodeType3', 'EdgeType2')]
}, {
    'name':
    'Heterogeneous Multilayer w/o Cycle (N=4,M=3)',
    'edges': [('NodeType1', 'NodeType2',
               'EdgeType1'), ('NodeType2', 'NodeType3', 'EdgeType2'),
              ('NodeType3', 'NodeType4', 'EdgeType3')]
}]

for schema in schemas:
    count_views(schema, 2, 8)
