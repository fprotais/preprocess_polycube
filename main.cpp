#include <string>
#include <iostream>
#include <fstream>
#include <ultimaille/all.h>
#include <OpenNL_psm/OpenNL_psm.h>


using namespace UM;
#define FOR(i, n) for(int i = 0; i < n; i++)

void preprocess(Tetrahedra& m) {
    std::cerr << "Split for overconstrained boundary faces " << std::endl;
    int last_checked = 0;
    int nb_of_change = 0;
    while (last_checked != m.ncells()) {
        VolumeConnectivity vec(m);
        std::cerr << "change " << nb_of_change++ << " : " << last_checked << " / "  << m.ncells() << std::endl;
        for (; last_checked < m.ncells(); last_checked++) {
            bool change = true;
            int cf0 = -1;
            int cf1 = -1;
            int c = last_checked;
            FOR(cf, 4) if (vec.adjacent[4 * c + cf] != -1) {
                if (cf0 == -1) cf0 = cf;
                else if (cf1 == -1) cf1 = cf;
                else change = false;
            }
            if (!change) continue;

            if (cf1 == -1) { // 3 facet on bnd
                std::cerr << "3 facets on boundary, splitting the last one " << std::endl;
                int opp_cf = vec.adjacent[4 * c + cf0];
                int c2 = opp_cf / 4;
                int cf2 = opp_cf % 4;
                vec3 new_vec = m.util.bary_facet(c, cf0);
                int i_n = m.points.create_points(1);
                m.points[i_n] = new_vec;

                std::vector<std::array<int, 4>> new_cells(6);
                FOR(j, 3) {
                    FOR(cfv, 3) new_cells[j][cfv + 1] = m.facet_vert(c, (cf0 + 1 + j) % 4, cfv);
                    new_cells[j][0] = i_n;
                    FOR(cfv, 3) new_cells[j + 3][cfv + 1] = m.facet_vert(c2, (cf2 + 1 + j) % 4, cfv);
                    new_cells[3 + j][0] = i_n;
                }

                int start = m.create_cells(4);
                FOR(cv, 4) {
                    m.vert(c, cv) = new_cells[0][cv];
                    m.vert(c2, cv) = new_cells[1][cv];
                    FOR(j, 4) m.vert(start + j, cv) = new_cells[2 + j][cv];
                }

            }
            else { // 2 facet on bnd
                std::cerr << "2 facets on boundary, splitting opp edge " << std::endl;
                std::array<int, 3> T1, T2;
                FOR(cfv, 3) T1[cfv] = m.facet_vert(c, cf0, cfv);
                FOR(cfv, 3) T2[cfv] = m.facet_vert(c, cf1, cfv);

                int v0 = -1;
                int v1 = -1;
                FOR(cf0v, 3) FOR(cf1v, 3) {
                    if (T1[cf0v] == T2[cf1v]) {
                        if (v0 == -1) v0 = T2[cf1v];
                        else if (T2[cf1v] != v0 && v1 == -1) v1 = T2[cf1v];
                    }
                }

                int he = vec.halfedge_from_verts(c, v0, v1);

                std::vector<int> cells_to_split = vec.halfedges_around_edge(he);

                int n_i = m.points.create_points(1);
                m.points[n_i] = 0.5 * (m.points[v0] + m.points[v1]);

                int nb = 0;
                std::vector<std::array<int, 4>> new_cells(2 * cells_to_split.size());
                std::vector<int> cells_changed;
                for (int cell_he : cells_to_split) {
                    int f0 = vec.cell_facet(vec.opposite_f(vec.next(cell_he)));
                    int f1 = vec.cell_facet(vec.opposite_f(vec.prev(cell_he)));
                    int cell = vec.cell(cell_he);
                    cells_changed.push_back(cell);
                    FOR(cfv, 3) new_cells[2 * nb][cfv + 1] = m.facet_vert(cell, f0, cfv);
                    new_cells[2 * nb][0] = n_i;
                    FOR(cfv, 3) new_cells[2 * nb + 1][cfv + 1] = m.facet_vert(cell, f1, cfv);
                    new_cells[2 * nb + 1][0] = n_i;
                    nb++;
                }
                int start = m.create_cells(cells_to_split.size());

                FOR(i, cells_to_split.size()) {
                    FOR(cv, 4) {
                        m.vert(cells_changed[i], cv) = new_cells[2 * i][cv];
                        m.vert(start + i, cv) = new_cells[2 * i + 1][cv];

                    }
                }
            }
            break;
        }


    }


    std::cerr << "Split for overconstrained edges " << std::endl;



    nb_of_change = 0;

    last_checked = 0;
    while (last_checked != m.ncells()) {
        VolumeConnectivity vec(m);
        std::cerr << "change " << nb_of_change << " : " << last_checked << " / " << m.ncells() << std::endl;

        std::vector<bool> bnd_vert(m.nverts(), false);
        FOR(c, m.ncells()) FOR(cf, 4) if (vec.adjacent[4 * c + cf] == -1) FOR(cfv, 3) {
            bnd_vert[m.facet_vert(c, cf, cfv)] = true;
        }

        for (; last_checked < m.ncells(); last_checked++) FOR(cf, 4) FOR(cfv, 3) {

            int c = last_checked;
            int he = vec.halfedge(c, cf, cfv);

            bool is_boundary = false;
            int around_e_cir = he;

            do { // rewind if boundary
                if (vec.opposite_c(around_e_cir) < 0) {
                    is_boundary = true;
                    break;
                }
                around_e_cir = vec.opposite_f(vec.opposite_c(around_e_cir));
            } while (around_e_cir != he);
            if (is_boundary) continue;
            if (!bnd_vert[vec.to(he)] || !bnd_vert[vec.from(he)]) continue;

            std::vector<int> cells_to_split = vec.halfedges_around_edge(he);

            int n_i = m.points.create_points(1);
            m.points[n_i] = 0.5 * (m.points[vec.to(he)] + m.points[vec.from(he)]);

            int nb = 0;
            std::vector<std::array<int, 4>> new_cells(2 * cells_to_split.size());
            std::vector<int> cells_changed;
            for (int cell_he : cells_to_split) {
                int f0 = vec.cell_facet(vec.opposite_f(vec.next(cell_he)));
                int f1 = vec.cell_facet(vec.opposite_f(vec.prev(cell_he)));
                int cell = vec.cell(cell_he);
                cells_changed.push_back(cell);
                FOR(cfv, 3) new_cells[2 * nb][cfv + 1] = m.facet_vert(cell, f0, cfv);
                new_cells[2 * nb][0] = n_i;
                FOR(cfv, 3) new_cells[2 * nb + 1][cfv + 1] = m.facet_vert(cell, f1, cfv);
                new_cells[2 * nb + 1][0] = n_i;
                nb++;
            }
            int start = m.create_cells(cells_to_split.size());

            FOR(i, cells_to_split.size()) {
                FOR(cv, 4) {
                    m.vert(cells_changed[i], cv) = new_cells[2 * i][cv];
                    m.vert(start + i, cv) = new_cells[2 * i + 1][cv];

                }
            }

            goto breaking;
        }

    breaking:
        nb_of_change++;
    }
}

void naive_flagging(const Tetrahedra& m, CellFacetAttribute<int>& flags) {
	std::cerr << "Computing naive flagging ...";
    const std::array<UM::vec3, 6> flag2normal = { UM::vec3(1,0,0), UM::vec3(-1,0,0), UM::vec3(0,1,0), UM::vec3(0,-1,0), UM::vec3(0,0,1), UM::vec3(0,0,-1) };
    VolumeConnectivity vec(m);
    FOR(cf, m.ncells() * 4) if (vec.adjacent[cf] == -1) {
        vec3 n = m.util.facet_normal(cf / 4, cf % 4);
        flags[cf] = 0;
        double best = flag2normal[0] * n;
        FOR(i, 5) {
            double d = flag2normal[i + 1] * n;
            if (d > best) {
                best = d;
                flags[cf] = i + 1;
            }
        }
    }
	std::cerr << " Done.\n";
}


int main(int argc, char** argv) {

	if (argc < 2) {
		std::cout << "Usage is: " << argv[0] << " tetmesh.ext resmesh.ext flagfile(optionnal)" << std::endl;
		std::cout << "exemple: " << argv[0] << " ../S1.mesh res.mesh res.flags" << std::endl;
		return 1;
	}
	std::string inputfile = argv[1];
	std::string outputfile = argv[2];
	
	Tetrahedra m;
	read_by_extension(inputfile, m);
	
	preprocess(m);
	write_by_extension(outputfile, m);

	
	if (argc == 4) {
        
		CellFacetAttribute<int> flag(m, -1);
		naive_flagging(m, flag);
		std::ofstream ofs(argv[3]);
		if (!ofs.is_open()) {
			std::cerr << "Failed opening of flags at : " << argv[3] << std::endl;
			abort();
		}
		FOR(cf, m.ncells() * 4) ofs << flag[cf] << "\n";
		ofs.close();
	}
	
	return 0;
}

