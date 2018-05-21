/**
 * From:
 * https://www.cc.gatech.edu/dimacs10/archive/task.shtml
 * https://www.cc.gatech.edu/dimacs10/archive/stream.cpp
 */

/*This program generates graphs that match the basic properties observed for the
 * computational task graphs of stream processing systems. In these graphs, the
 * vertices represent the kernels and the edges represent a continous
 * data-stream movement from one kernel to another.*/

#include <algorithm>
#include <climits>
#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <queue>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <vector>
using namespace std;

#define SEED 135792468

class edge {
private:
  unsigned int src;
  unsigned int tgt;

public:
  edge(unsigned int _source, unsigned int _target)
      : src(_source), tgt(_target) {}
  unsigned int source() { return src; }
  unsigned int target() { return tgt; }
  void set_source(unsigned int _source) { src = _source; }
};

class edgeWithWeight {
private:
  unsigned int src;
  unsigned int tgt;
  unsigned int src_label;
  double wgt;

public:
  edgeWithWeight(unsigned int _source, unsigned int _target,
                 unsigned int _source_label)
      : src(_source), tgt(_target), src_label(_source_label) {}
  edgeWithWeight(unsigned int _source, unsigned int _target, double _wgt)
      : src(_source), tgt(_target), wgt(_wgt) {}
  edgeWithWeight() {
    src = 0;
    tgt = 0;
    src_label = 0;
    wgt = 0;
  }
  double weight() { return wgt; }
  unsigned int source() { return src; }
  unsigned int target() { return tgt; }
  unsigned int source_label() { return src_label; }
  void set_weight(double _wgt) { wgt = _wgt; }
  void set_source_label(unsigned int _src_label) { src_label = _src_label; }
};

bool edgeWithWeightLabel(edgeWithWeight e1, edgeWithWeight e2) {
  if (e1.source_label() < e2.source_label())
    return true;
  if (e1.source_label() > e2.source_label())
    return false;
  if (e1.source() < e2.source())
    return true;
  if (e1.source() > e2.source())
    return false;
  if (e1.target() < e2.target())
    return true;
  return false;
}

bool edgeWithWeightSource(edgeWithWeight e1, edgeWithWeight e2) {
  if (e1.source() < e2.source())
    return true;
  if (e1.source() > e2.source())
    return false;
  if (e1.target() < e2.target())
    return true;
  return false;
}

bool edgeTarget(edge e1, edge e2) {
  if (e1.target() < e2.target())
    return true;
  if (e1.target() > e2.target())
    return false;
  if (e1.source() < e2.source())
    return true;
  return false;
}

bool edgeSource(edge e1, edge e2) {
  if (e1.source() < e2.source())
    return true;
  if (e1.source() > e2.source())
    return false;
  if (e1.target() < e2.target())
    return true;
  return false;
}

bool edge_equal(edge e1, edge e2) {
  return ((e1.source() == e2.source()) && (e1.target() == e2.target()));
}

void compute_inout_degree(unsigned int n_nodes, vector<edge> &graph,
                          vector<unsigned int> &indegree,
                          vector<unsigned int> &outdegree) {
  outdegree.clear();
  indegree.clear();
  sort(graph.begin(), graph.end(), edgeSource);
  vector<edge>::iterator iter_ver = graph.begin();
  for (int i = 0; i < n_nodes; i++) {
    unsigned int i_outedges = 0;
    while ((iter_ver != graph.end()) && ((*iter_ver).source() == i)) {
      i_outedges++;
      ++iter_ver;
    }
    outdegree.push_back(i_outedges);
  }
  sort(graph.begin(), graph.end(), edgeTarget);

  iter_ver = graph.begin();
  for (int i = 0; i < n_nodes; i++) {
    unsigned int i_inedges = 0;
    while ((iter_ver != graph.end()) && ((*iter_ver).target() == i)) {
      i_inedges++;
      ++iter_ver;
    }
    indegree.push_back(i_inedges);
  }
}

void compute_longest_paths(unsigned int n_nodes,
                           vector<edge> &edge_vector_unique,
                           vector<unsigned int> &indegree,
                           vector<unsigned int> &longest_path) {
  // Recompute indegree
  sort(edge_vector_unique.begin(), edge_vector_unique.end(), edgeTarget);
  vector<edge>::iterator iter = edge_vector_unique.begin();
  for (int i = 0; i < n_nodes; i++) {
    unsigned int i_inedges = 0;
    while ((iter != edge_vector_unique.end()) && ((*iter).target() == i)) {
      i_inedges++;
      ++iter;
    }
    indegree.push_back(i_inedges);
  }
  // cout<<"Indegree re-computed"<<endl;cout.flush();

  //  for(int i=0;i<indegree.size();i++)
  //    cout<<"Indegree of vertex "<<i<<" is "<<indegree[i]<<endl;
  //  cout.flush();

  // Compute Topological Order = Longest Path of vertices
  sort(edge_vector_unique.begin(), edge_vector_unique.end(), edgeSource);

  vector<vector<int>> neighb;

  iter = edge_vector_unique.begin();
  for (int i = 0; i < n_nodes; i++) {
    vector<int> i_neighb;
    while ((iter != edge_vector_unique.end()) && ((*iter).source() == i)) {
      i_neighb.push_back((*iter).target());
      ++iter;
    }
    neighb.push_back(i_neighb);
  }
  // cout<<"Neighbors pushed back"<<endl;cout.flush();

  for (int i = 0; i < n_nodes; i++)
    longest_path.push_back(0);

  queue<int> sources;
  for (int i = 0; i < n_nodes; i++) {
    if (indegree[i] == 0)
      sources.push(i);
  }
  // cout<<"Sources pushed back"<<endl;cout.flush();
  while (!sources.empty()) {
    int v = sources.front();
    sources.pop();
    for (vector<int>::iterator iter1 = neighb[v].begin();
         iter1 != neighb[v].end(); ++iter1) {
      int w = *iter1;
      if (longest_path[w] < longest_path[v] + 1)
        longest_path[w] = longest_path[v] + 1;
      indegree[w]--;
      if (indegree[w] == 0)
        sources.push(w);
    }
  }
  // cout<<"Top sort computed"<<endl;cout.flush();
}

char hex_char(int e) {
  if (e == 0)
    return '0';
  if (e == 1)
    return '1';
  if (e == 2)
    return '2';
  if (e == 3)
    return '3';
  if (e == 4)
    return '4';
  if (e == 5)
    return '5';
  if (e == 6)
    return '6';
  if (e == 7)
    return '7';
  if (e == 8)
    return '8';
  if (e == 9)
    return '9';
  if (e == 10)
    return 'A';
  if (e == 11)
    return 'B';
  if (e == 12)
    return 'C';
  if (e == 13)
    return 'D';
  if (e == 14)
    return 'E';
  if (e == 15)
    return 'F';
}

char *print_color(double col) {
  string out_str;
  out_str = " [color=\"#";
  char color_str[6];
  int f = floor(col * 255.99);
  int f1 = f / 16;
  int f2 = f % 16;
  int f3 = (255 - f) / 16;
  int f4 = (255 - f) % 16;
  color_str[0] = hex_char(f1);
  color_str[1] = hex_char(f2);
  color_str[2] = hex_char(f3);
  color_str[4] = hex_char(f4);
  color_str[3] = 'F';
  color_str[5] = 'F';
  for (int i = 0; i < 6; i++)
    out_str += color_str[i];
  out_str += "\"]";
  return strdup(out_str.c_str());
}

int main(int argv, char *argc[]) {

  // char filename[100];
  string file_base;
  unsigned int n;
  unsigned int m;

  if (argv < 3) {
    cout << "Usage: ./generate n_nodes filename" << endl;
    cout << "The generated graph with n_nodes vertices is written as directed "
            "graph in <filename>.dgr, as undirected graph in <filename>.gr"
         << endl;
    cout.flush();
    return -1;
  } else {
    n = atoi(argc[1]);
    file_base = argc[2];
  }

  vector<edge> edge_vector;

  timeval tt;
  gettimeofday(&tt, NULL);
  // srandom(tt.tv_usec + tt.tv_sec*1000000);
  // changed to allow for reproducibility
  srandom(SEED);

  // First generate the core of the graph with O(n^{1/3}) vertices and
  // O(n^{2/3}) edges
  // The vertices in this graph will later be divided into splits and joins. The
  // filters will be added even later.
  int k = 1;
  unsigned int n_nodes = ceil(2 * k * pow(n, 0.3333));
  unsigned int m_edges = ceil(pow(n, 0.6666));

  // Part 1 of generating core
  // Generate series-parallel directed acyclic multi-graph
  unsigned int curr_n_nodes = 2;
  unsigned int curr_edges = random() % (n_nodes / 2 - 1) + 1;
  for (int i = 0; i < curr_edges; ++i)
    edge_vector.push_back(edge(0, 1));
  while (curr_n_nodes < n_nodes) {
    unsigned int select_edge = random() % edge_vector.size();
    edge_vector.push_back(
        edge(edge_vector[select_edge].source(), curr_n_nodes));
    edge_vector.push_back(
        edge(curr_n_nodes + 1, edge_vector[select_edge].target()));
    curr_edges = random() % (n_nodes / 2 - 1) + 1;
    for (int i = 0; i < curr_edges; ++i)
      edge_vector.push_back(edge(curr_n_nodes, curr_n_nodes + 1));
    edge_vector.erase(edge_vector.begin() + select_edge);
    curr_n_nodes += 2;
  }
  n_nodes = curr_n_nodes;

  vector<unsigned int> indegree;
  vector<unsigned int> outdegree;

  vector<unsigned int> longest_path;

  compute_longest_paths(n_nodes, edge_vector, indegree, longest_path);

  // Add a random sequence of edges to the series-parallel core
  for (int i = 0; i < n_nodes; i++) {
    unsigned int first = random() % n_nodes;
    unsigned int second = random() % n_nodes;
    while (longest_path[first] == longest_path[second]) {
      first = random() % n_nodes;
      second = random() % n_nodes;
    }
    if (longest_path[first] < longest_path[second])
      edge_vector.push_back(edge(first, second));
    else
      edge_vector.push_back(edge(second, first));
  }

  compute_inout_degree(n_nodes, edge_vector, indegree, outdegree);
  sort(edge_vector.begin(), edge_vector.end(), edgeSource);

  // Decompose vertices with multi-in and multi-out edges into split-join pair
  vector<edge>::iterator iter = edge_vector.begin();
  vector<edge> add_edges;
  int n_nodes_fix = n_nodes;

  for (int i = 0; i < n_nodes_fix; i++) {
    if ((indegree[i] > 1) && (outdegree[i] > 1)) {
      while ((iter != edge_vector.end()) && ((*iter).source() < i))
        ++iter;
      while ((iter != edge_vector.end()) && ((*iter).source() == i)) {
        (*iter).set_source(n_nodes);
        ++iter;
      }
      add_edges.push_back(edge(i, n_nodes));
      n_nodes++;
    }
  }

  //  cout<<"Nodes with multi-in and multi-out divided"<<endl;cout.flush();

  // Add the edges between the splitted vertices
  for (vector<edge>::iterator add_iter = add_edges.begin();
       add_iter != add_edges.end(); ++add_iter)
    edge_vector.push_back(*add_iter);

  // cout<<"After adding the in to out edges"<<endl;cout.flush();

  // This can be the place to add some high degree splits in the beginning and
  // high degree joins in the end.

  indegree.clear();
  longest_path.clear();

  compute_longest_paths(n_nodes, edge_vector, indegree, longest_path);

  // for(int i=0;i<n_nodes;i++)
  // cout<<"Longest path to "<<i<<" is "<<longest_path[i]<<endl;

  // Replace edges by paths ensuring that path-difference is very small, if at
  // all
  long sum_edge_diff = 0;

  for (vector<edge>::iterator iter1 = edge_vector.begin();
       iter1 != edge_vector.end(); ++iter1) {
    sum_edge_diff +=
        longest_path[(*iter1).target()] - longest_path[(*iter1).source()];
  }

  sort(edge_vector.begin(), edge_vector.end(), edgeSource);

  // cout<<"sum_edge_diff = "<<sum_edge_diff<<endl;cout.flush();
  double filters_per_edge = ((double)(n - n_nodes)) / sum_edge_diff;
  // cout<<"filters_per_edge = "<<filters_per_edge<<endl;cout.flush();
  vector<edge> unw_graph;
  for (vector<edge>::iterator iter1 = edge_vector.begin();
       iter1 != edge_vector.end(); ++iter1) {
    int this_edge_diff =
        longest_path[(*iter1).target()] - longest_path[(*iter1).source()];
    double filters_this_edge = this_edge_diff * filters_per_edge;
    int add_filters_this_edge = (int)floor(filters_this_edge + 0.5);
    if (add_filters_this_edge == 0)
      unw_graph.push_back(edge((*iter1).source(), (*iter1).target()));
    else {
      unw_graph.push_back(edge((*iter1).source(), n_nodes));
      n_nodes++;
      for (int i = 1; i < add_filters_this_edge; ++i, ++n_nodes)
        unw_graph.push_back(edge(n_nodes - 1, n_nodes));
      unw_graph.push_back(edge(n_nodes - 1, (*iter1).target()));
    }
  }

  // Remove duplicates from unweighted graph
  sort(unw_graph.begin(), unw_graph.end(), edgeTarget);
  vector<edge>::iterator end_edge_vector =
      std::unique(unw_graph.begin(), unw_graph.end(), edge_equal);
  vector<edge> unw_graph_unique =
      vector<edge>(unw_graph.begin(), end_edge_vector);
  edge_vector.clear();
  unw_graph.clear();

  // Divide splits and joins into sub-categories and assign weights to edges. We
  // make 35% of splits as copying splits, rest divide the input equally among
  // the output. All joins add the weights of input to output. 10% of filters
  // reduce the weights by half, mostly in earlier part of the DAG. StreamIT
  // benchmark only says that 66% of filters have same multiplicity, which means
  // that 2/3 of filters get the same input rate.

  // Compute a new topological order of graph and start assigning weights based
  // on the above criteria.
  longest_path.clear();
  indegree.clear();
  compute_longest_paths(n_nodes, unw_graph_unique, indegree, longest_path);

  compute_inout_degree(n_nodes, unw_graph_unique, indegree, outdegree);

  sort(unw_graph_unique.begin(), unw_graph_unique.end(), edgeTarget);

  vector<edgeWithWeight> graph;
  for (vector<edge>::iterator graph_iter = unw_graph_unique.begin();
       graph_iter != unw_graph_unique.end(); ++graph_iter) {
    edgeWithWeight e((*graph_iter).source(), (*graph_iter).target(),
                     longest_path[(*graph_iter).source()]);
    graph.push_back(e);
  }

  sort(graph.begin(), graph.end(), edgeWithWeightLabel);
  multimap<int, edgeWithWeight> target_tree;

  vector<edgeWithWeight>::iterator wgraph_iter;
  for (wgraph_iter = graph.begin(); wgraph_iter != graph.end(); ++wgraph_iter) {
    int curr_vertex = (*wgraph_iter).source();
    if ((indegree[curr_vertex] <= 1) && (outdegree[curr_vertex] == 1)) {
      double w = 1.0;
      if (indegree[curr_vertex] == 1) {
        multimap<int, edgeWithWeight>::iterator curr_vertex_inedge =
            target_tree.find(curr_vertex);
        edgeWithWeight e1;
        e1 = (*curr_vertex_inedge).second;
        w = e1.weight();
        // cout<<"Received weight "<<w<<endl;cout.flush();
        target_tree.erase(curr_vertex_inedge);
      }
      if (random() % (longest_path[curr_vertex] / 4 + 1) == 0)
        (*wgraph_iter).set_weight(w / 2);
      else
        (*wgraph_iter).set_weight(w);
      edgeWithWeight e = *wgraph_iter;
      target_tree.insert(pair<int, edgeWithWeight>((*wgraph_iter).target(), e));
    }
    if ((indegree[curr_vertex] <= 1) && (outdegree[curr_vertex] > 1)) {
      double w = 1;
      // cout<<"w = "<<w<<" outdegree[curr_vertex] =
      // "<<outdegree[curr_vertex]<<" ratio =
      // "<<w/outdegree[curr_vertex]<<endl;cout.flush();
      if (indegree[curr_vertex] == 1) {
        multimap<int, edgeWithWeight>::iterator curr_vertex_inedge =
            target_tree.find(curr_vertex);
        edgeWithWeight e1;
        e1 = (*curr_vertex_inedge).second;
        w = e1.weight();
        target_tree.erase(curr_vertex_inedge);
      }
      bool IsCopy;
      if (random() % 100 < 35)
        IsCopy = true;
      else
        IsCopy = false;
      for (int i = 0; i < outdegree[curr_vertex]; i++) {
        if (IsCopy)
          (*wgraph_iter).set_weight(w);
        else
          (*wgraph_iter).set_weight(w / outdegree[curr_vertex]);
        edgeWithWeight e = *wgraph_iter;
        target_tree.insert(
            pair<int, edgeWithWeight>((*wgraph_iter).target(), e));
        if (wgraph_iter != graph.end())
          ++wgraph_iter;
      }
      --wgraph_iter;
    }
    if (indegree[curr_vertex] > 1) {
      double w = 0;
      pair<multimap<int, edgeWithWeight>::iterator,
           multimap<int, edgeWithWeight>::iterator>
          piter = target_tree.equal_range(curr_vertex);
      for (multimap<int, edgeWithWeight>::iterator curr_vertex_inedge =
               piter.first;
           curr_vertex_inedge != piter.second;) {
        edgeWithWeight e1;
        e1 = (*curr_vertex_inedge).second;
        w += e1.weight();
        target_tree.erase(curr_vertex_inedge++);
      }
      (*wgraph_iter).set_weight(w);
      edgeWithWeight e = *wgraph_iter;
      target_tree.insert(pair<int, edgeWithWeight>((*wgraph_iter).target(), e));
    }
  }

  double min_weight = INT_MAX;
  for (wgraph_iter = graph.begin(); wgraph_iter != graph.end(); ++wgraph_iter) {
    if (min_weight > (*wgraph_iter).weight())
      min_weight = (*wgraph_iter).weight();
  }

  double make_int = log(1 / min_weight) / log(2);
  int max_label = 0;
  // cout<<"Min weight is "<<min_weight<<" and make_int is
  // "<<make_int<<endl;cout.flush();

  // Cap it up so its solvable by partitioners like METIS
  if (make_int > 20.0)
    make_int = 20.0;

  for (wgraph_iter = graph.begin(); wgraph_iter != graph.end(); ++wgraph_iter) {
    double new_wgt = pow(2, make_int) * ((*wgraph_iter).weight());
    (*wgraph_iter).set_source_label(floor(new_wgt + 0.5));
    if ((*wgraph_iter).source_label() == 0)
      (*wgraph_iter).set_source_label(1);
    //  cout<<(*wgraph_iter).weight()<<"
    //  "<<(*wgraph_iter).source_label()<<endl;cout.flush();
    if (max_label < (*wgraph_iter).source_label())
      max_label = (*wgraph_iter).source_label();
  }

  // Output in Dot format
  string dot_dir_file_name = file_base + ".dot.dir.gr";
  ofstream dot_dir_graph_file;
  dot_dir_graph_file.open(dot_dir_file_name.c_str());
  //    cout<<endl<<endl<<"Creating Dot File"<<endl<<endl;
  dot_dir_graph_file << "digraph G {" << endl;
  dot_dir_graph_file << "rotate=90;" << endl;
  for (wgraph_iter = graph.begin(); wgraph_iter != graph.end(); ++wgraph_iter) {
    dot_dir_graph_file << (*wgraph_iter).source() << " -> "
                       << (*wgraph_iter).target();
    dot_dir_graph_file << print_color(((double)(*wgraph_iter).source_label()) /
                                      max_label);
    dot_dir_graph_file << ";" << endl;
  }
  dot_dir_graph_file << "}" << endl;
  dot_dir_graph_file.close();

  m = 0;
  vector<edgeWithWeight> undGraph;
  wgraph_iter = graph.begin();
  for (; wgraph_iter != graph.end(); ++wgraph_iter) {
    undGraph.push_back(edgeWithWeight((*wgraph_iter).source(),
                                      (*wgraph_iter).target(),
                                      (*wgraph_iter).source_label()));
    undGraph.push_back(edgeWithWeight((*wgraph_iter).target(),
                                      (*wgraph_iter).source(),
                                      (*wgraph_iter).source_label()));
    ++m;
  }

  sort(graph.begin(), graph.end(), edgeWithWeightSource);

  // Output as Directed graph
  string dir_file_name = file_base + ".dir.gr";
  ofstream dir_graph_file;
  dir_graph_file.open(dir_file_name.c_str());

  dir_graph_file << n_nodes << " " << m << " 11 1" << endl;
  cout.flush();
  wgraph_iter = graph.begin();
  for (int i = 0; i < n_nodes; i++) {
    int node_wgt = 1 + pow(longest_path[i], 2);
    dir_graph_file << node_wgt << " ";
    while ((wgraph_iter != graph.end()) && ((*wgraph_iter).source() == i)) {
      // int edge_wgt = 1+((int)
      // (((*wgraph_iter).weight()-min_weight)/(max_weight-min_weight)*(INT_MAX-1)));
      dir_graph_file << " " << (*wgraph_iter).target() + 1 << " "
                     << (*wgraph_iter).source_label();
      ++wgraph_iter;
    }
    dir_graph_file << endl;
    cout.flush();
  }

  dir_graph_file.close();

  // Output in METIS format (also DIMACS challenge format)
  sort(undGraph.begin(), undGraph.end(), edgeWithWeightSource);
  string und_file_name = file_base + ".gr";
  ofstream und_graph_file;
  und_graph_file.open(und_file_name.c_str());

  und_graph_file << n_nodes << " " << m << " 11 1" << endl;
  cout.flush();
  wgraph_iter = undGraph.begin();
  for (int i = 0; i < n_nodes; i++) {
    int node_wgt = 1 + pow(longest_path[i], 2);
    und_graph_file << node_wgt << " ";
    while ((wgraph_iter != undGraph.end()) && ((*wgraph_iter).source() == i)) {
      // int edge_wgt = 1+((int)
      // (((*wgraph_iter).weight()-min_weight)/(max_weight-min_weight)*(INT_MAX-1)));
      und_graph_file << " " << (*wgraph_iter).target() + 1 << " "
                     << (*wgraph_iter).source_label();
      ++wgraph_iter;
    }
    und_graph_file << endl;
    cout.flush();
  }

  und_graph_file.close();
  // printf("%s", filename);

  return 0;
}
