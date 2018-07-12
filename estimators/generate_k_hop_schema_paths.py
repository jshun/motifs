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

#edges = [
#    ('Job', 'File', 'writesTo'),
#    ('File', 'Job', 'isReadBy')
#]
edges = [
    ('Node', 'Node', 'EdgeType1'),
    ('Node', 'Node', 'EdgeType2'),
    ('Node', 'Node', 'EdgeType3')
]
k = 2
print 'paths:'
for i, path in enumerate(generate_k_hop_schema_paths(edges, [], k, k)):
    print "%d:\n%s\n" % (i, path)
