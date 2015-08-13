#include <string>
#include <fstream>
#include <istream>
#include <sstream>
#include <unordered_map>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/iteration_macros.hpp>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/variables_map.hpp>

#include <glpk.h>

// The following function comes from stackoverflow. Code by user sehe, adapted
// to also read vertex weights.
// http://stackoverflow.com/questions/30415388/how-to-read-dimacs-vertex-coloring-graphs-in-c
template <typename Graph>
bool read_graph(std::istream& dimacs, Graph& g) {
        size_t vertices = 0, edges = 0;

        std::string line;
        while (getline(dimacs, line)) {
                std::istringstream iss(line);
                char ch;

                if (iss >> ch) {
                        size_t from, to;
                        std::string format;

                        switch(ch) {
                        case 'c': break;
                        case 'p':
                                if (vertices||edges) return false;
                                if (iss >> format >> vertices >> edges) {
                                        if ("edge" != format) return false;
                                }
                                break;
                        case 'n':
                                if (!vertices) return false;
                                size_t v, weight;
                                if (iss >> v >> weight)
                                        g[v - 1].weight = weight;
                                break;
                        case 'e':
                                if (edges-- && (iss >> from >> to) &&
                                    (add_edge(from-1, to-1, g).second))
                                        break;
                        default:
                                return false;
                        }
                }
        }

        return !(edges || !dimacs.eof());
}

template<typename Graph>
class MaxWeightIndependentSet
{
        typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
        typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
        typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIndex;

        typedef boost::property<boost::edge_weight_t, double> EdgeWeightProperty;
        typedef boost::adjacency_list<boost::setS,
                                      boost::vecS,
                                      boost::undirectedS,
                                      boost::no_property,
                                      EdgeWeightProperty> EdgeWeightGraph;
        typedef typename boost::graph_traits<EdgeWeightGraph>::vertex_descriptor EWGVertex;
        typedef typename boost::graph_traits<EdgeWeightGraph>::edge_descriptor EWGEdge;

public:
        MaxWeightIndependentSet(const Graph& G, bool separate_odd_cycles = true)
         : G_(G), separate_odd_cycles_(separate_odd_cycles) {
                vertex_index_ = get(boost::vertex_index, G_);
                build_bipartite_graph();
                setup_lp_problem();
        }

        void solve() {
                run_cut_plane_method();
        }

        void print_solution() {
                std::cout << "Found stable set with weight: "
                          << glp_mip_obj_val(lp_) << std::endl;

                int n = glp_get_num_cols(lp_);

                std::cout << "Vertices: " << std::endl;
                for (int i = 1; i <= n; i++)
                        if (glp_mip_col_val(lp_, i))
                                std::cout << i - 1 << ", ";
                std::cout << std::endl;
        }

        void print_edge_weigths() {
                std::cout << "Edge weights:" << std::endl;
                BGL_FORALL_EDGES_T(e, H_, EdgeWeightGraph)
                        std::cout << get(weightmap_, e) << ", ";
                std::cout << std::endl;
        }

private:
        void build_bipartite_graph() {
                BGL_FORALL_VERTICES_T(v, G_, Graph) {
                        Vertex v1, v2;

                        map_to_H1_[v] = add_vertex(H_);
                        map_to_H2_[v] = add_vertex(H_);
                        map_to_G_[map_to_H1_[v]] = v;
                        map_to_G_[map_to_H2_[v]] = v;
                }

                BGL_FORALL_VERTICES_T(v, G_, Graph) {
                        typename Graph::adjacency_iterator it, end;

                        boost::tie(it, end) = adjacent_vertices(v, G_);
                        for (; it != end; ++it) {
                                const Vertex& u = *it;

                                add_edge(map_to_H1_[v], map_to_H2_[u], H_);
                                add_edge(map_to_H2_[v], map_to_H1_[u], H_);
                        }
                }
        }

        void update_bipartite_graph_weights() {
                BGL_FORALL_EDGES_T(e, H_, EdgeWeightGraph) {
                        double weight, x_u, x_v;
                        auto u = source(e, H_);
                        auto v = target(e, H_);

                        x_u = glp_get_col_prim(lp_, vertex_index_[map_to_G_[u]] + 1);
                        x_v = glp_get_col_prim(lp_, vertex_index_[map_to_G_[v]] + 1);

                        weight = (1 - x_u - x_v) / 2;
                        if (weight < 0)
                                weight = 0;

                        put(weightmap_, e, weight);
                }
        }

        void setup_lp_problem() {
                lp_ = glp_create_prob();
                std::vector<int> row_index;
                std::vector<double> row_coefficients;

                glp_add_rows(lp_, num_edges(G_));
                glp_add_cols(lp_, num_vertices(G_));
                glp_set_obj_dir(lp_, GLP_MAX);

                // GLPK indexes start at 1, and it access elements 1..n of the coefficient
                // vectors, leaving element zero unused. So we need to allocate space
                // for n + 1 elements.
                row_index.reserve(num_vertices(G_) + 1);
                row_coefficients.reserve(num_vertices(G_) + 1);

                // Setup vertex variables
                BGL_FORALL_VERTICES_T(v, G_, Graph) {
                        glp_set_obj_coef(lp_, vertex_index_[v] + 1, G_[v].weight);
                        glp_set_col_bnds(lp_, vertex_index_[v] + 1, GLP_DB, 0.0, 1.0);
                        glp_set_col_kind(lp_, vertex_index_[v] + 1, GLP_BV);
                }

                int eidx=1;
                BGL_FORALL_EDGES_T(e, G_, Graph) {
                        Vertex u, v;

                        BGL_FORALL_VERTICES_T(v, G_, Graph) {
                                row_index[vertex_index_[v] + 1] = vertex_index_[v] + 1;
                                row_coefficients[vertex_index_[v] + 1] = 0.0;
                        }

                        u = source(e, G_);
                        v = target(e, G_);

                        row_coefficients[vertex_index_[u] + 1] = 1.0;
                        row_coefficients[vertex_index_[v] + 1] = 1.0;

                        glp_set_mat_row(lp_, eidx, num_vertices(G_),
                                        &row_index[0], &row_coefficients[0]);
                        glp_set_row_bnds(lp_, eidx, GLP_UP, 0.0, 1.0);
                        eidx++;
                }
        }

        void find_initial_solution() {
                // TODO: use a maximal stable set or clique covering
                glp_simplex(lp_, NULL);
        }

        void add_odd_cycle_inquality(Vertex v, std::vector<EWGVertex> p,
                                     std::unordered_map<Vertex, bool>& covered) {
                int new_row = glp_add_rows(lp_, 1);
                std::vector<int> indexes(num_vertices(G_) + 1);
                std::vector<double> coefficients(num_vertices(G_) + 1, 0.0);
                std::list<Vertex> cycle;

                BGL_FORALL_VERTICES_T(v, G_, Graph)
                        indexes[vertex_index_[v] + 1] = vertex_index_[v] + 1;

                EWGVertex u = map_to_H2_[v];
                while (u != map_to_H1_[v]) {
                        cycle.push_back(map_to_G_[u]);
                        covered[map_to_G_[u]] = true;
                        u = p[u];
                }

                // Check for repeated vertices. This can happen since we are
                // not taking only the shortest path, but any path that with
                // the appropriate distance. In that case, it is possible that
                // cycle is not actually a cycle.
                for (auto it1 = cycle.begin(); it1 != cycle.end(); ++it1) {
                        auto it2 = it1;
                        for (++it2; it2 != cycle.end(); ++it2) {
                                if (*it1 == *it2) {
                                        cycle.erase(cycle.begin(), it1);
                                        cycle.erase(it2, cycle.end());
                                        break;
                                }
                        }
                }

                for (auto u: cycle)
                        coefficients[vertex_index_[u] + 1] = 1.0;

                assert(cycle.size() % 2);

                glp_set_mat_row(lp_, new_row, num_vertices(G_),
                                &indexes[0], &coefficients[0]);
                glp_set_row_bnds(lp_, new_row, GLP_UP, 0.0,
                                 (double) (cycle.size() - 1) / 2);
        }

        void separate_odd_cycles() {
                int added_inequalities = 0;
                std::unordered_map<Vertex, bool> covered;

                update_bipartite_graph_weights();

                BGL_FORALL_VERTICES_T(v, G_, Graph)
                        covered[v] = false;

                BGL_FORALL_VERTICES_T(v, G_, Graph) {
                        std::vector<EWGVertex> p(num_vertices(H_));
                        std::vector<double> d(num_vertices(H_));

                        dijkstra_shortest_paths(H_, map_to_H1_[v],
                                                predecessor_map(boost::make_iterator_property_map(p.begin(), get(boost::vertex_index, H_))).
                                                distance_map(boost::make_iterator_property_map(d.begin(), get(boost::vertex_index, H_))));

                        // If the weight of a cycle is < 0.5, then it's odd
                        // cycle inequality is violated.
                        if (!covered[v] && d[map_to_H2_[v]] + 0.000001 < 0.5) {
                                add_odd_cycle_inquality(v, p, covered);
                                added_inequalities++;
                        }
                }

                std::cout << "Added " << added_inequalities << " inequalities." << std::endl;

                if (added_inequalities == 0)
                        separate_odd_cycles_ = false;
        }

        void mip_callback(glp_tree *T) {
                switch (glp_ios_reason(T)) {
                        case GLP_ICUTGEN:
                                if (separate_odd_cycles_)
                                        separate_odd_cycles();
                                break;
                        default:
                                break;
                }
        }

        static void static_mip_callback(glp_tree *T, void *data) {
                MaxWeightIndependentSet<Graph> *solver =
                        static_cast<MaxWeightIndependentSet<Graph> *>(data);
                solver->mip_callback(T);
        }

        void run_cut_plane_method() {
                glp_iocp parm;

                glp_init_iocp(&parm);
                parm.cb_func = MaxWeightIndependentSet<Graph>::static_mip_callback;
                parm.cb_info = this;

                find_initial_solution();

                glp_intopt(lp_, &parm);
        }

        const Graph& G_;
        VertexIndex vertex_index_;

        // A bipartite graph where each partition is a copy of the vertices
        // of G, and two vertices are adjacent if their corresponding vertices
        // in the original graph form an edge.
        EdgeWeightGraph H_;

        // Map vertices from G to H
        std::unordered_map<Vertex, EWGVertex> map_to_H1_, map_to_H2_;
        // Map vertices from H back to G
        std::unordered_map<EWGVertex, Vertex> map_to_G_;

        boost::property_map<EdgeWeightGraph, boost::edge_weight_t>::type weightmap_;

        glp_prob *lp_;

        // Algorithm knobs
        bool separate_odd_cycles_;
};

struct vertex_properties {
        double weight;
};

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, struct vertex_properties> Graph;

int main(int argc, char *argv[])
{
        Graph G;
        std::ifstream input_file;

        namespace po = boost::program_options;

        po::options_description desc("Allowed options");
        desc.add_options()
                ("odd-cycle", "use odd cycle separation")
                ("input-file", po::value<std::string>(), "input file")
                ;

        po::positional_options_description p;
        p.add("input-file", -1);

        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
        po::notify(vm);

        bool separate_odd_cycles = vm.count("odd-cycle");

        if (!vm.count("input-file")) {
                std::cout << desc << std::endl;
                return 1;
        }

        input_file.open(vm["input-file"].as<std::string>());
        if (!read_graph(input_file, G)) {
                std::cerr << "Failed to parse input file." << std::endl;
                return 1;
        }

        MaxWeightIndependentSet<Graph> solver(G, separate_odd_cycles);
        solver.solve();
        solver.print_solution();

        return 0;
}
