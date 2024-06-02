/*
* This code constructs a BSP (binary space partitioning) tree from geometrical room data.
* It translates the current geometry representation into the representation used by the 
* geometry splitting algorithm. It then recursively builds the tree by calculating the 
* Ranta-Eskola criterion and using the wall selected by the criterion to split the geometry 
* in the subtree. It then calculates the geometrically motivated data structures 
* direct_reflectables for the newly split up walls, harmonises identifiers and returns the
* resulting tree structure.
* This software was written in the course of my Master's Thesis at the  Professorship of Audio 
* Information Processing at TUM (Technical University of Munich).
*/

#pragma once
#include <iterator>
#include <map>
#include <utility>
#include <vector>
#include <armadillo>

#ifdef _DEBUG
#include <iostream>
#endif

#include "material.h"
#include "read_obj.h"
#include "wall.h"

// space partitioning includes (https://github.com/erich666/GraphicsGems/blob/master/gemsv/ch7-4/)
#include "spatial-partitioning/polygon.h"
#include "spatial-partitioning/dedge.h"
#include "spatial-partitioning/vector.h"
#include "spatial-partitioning/list.h"

namespace rts {

	struct BSPNode {
		const std::vector<rts::wall*> node_walls;
		const BSPNode* front;
		const BSPNode* back;
		const bool leaf_node = false;
	};

	class room_model {
	public:
		// Used only during program initialization: load parameters from config files, disable walls accordingly.
		std::vector<unsigned int> walls_to_disable;

		// Walls used by the algorithms
		std::vector<rts::wall*> walls;

		// Walls created by the BSP algorithm
		std::vector<rts::wall*> pwalls_BSP;
		std::vector<rts::wall> walls_BSP;

		// Plane-Polygon Map as per: PhD_Thesis_Schröder_Physically_based_real_time_auralization.pdf
		std::vector<std::vector<rts::wall*>> plane_polygon_map;

		// BSP-Tree + resulting height
		const BSPNode* bsp_tree;
		int bsp_tree_height = 0;

		/*!
			Construct and return Binary partitioning tree, to accelerate IS wall lookup / intersection

			/param polygons
			The walls used to create the PolygonSpatial data structure used by the spatial partitioning algorithm
			/param threshold 
			The threshold for splitting up the room geometry, see Ranta-Eskola criterion
		*/ 
		const BSPNode* set_up_room_model(std::vector<rts::wall> polygons, const double threshold) {
			
			// Transmogrify the rts::wall data structure to PolygonSpatial data structure for use in the algorithm
			std::vector<PolygonSpatial*> polygonSpatialPartitioning = construct_polygonspatial_model(polygons);
			for (int i = 0; i < polygons.size(); i++) {
				walls.push_back(&(polygons[i]));
			}

			// Use the newly built walls (PolygonSpatial) to create the Binary tree structure 
			bsp_tree = build_BSP(polygonSpatialPartitioning, pwalls_BSP, threshold);

			// Give new IDs + harmonise identifiers
			for (int i = 0; i < walls.size(); i++) {
				int monotone_id_increment = 0;
				for (auto& j : pwalls_BSP) {
					if (j->parent_id == i) {
						j->setID(j->parent_id * 1000 + monotone_id_increment);
						monotone_id_increment++;
					}
				}
			}

			// Construct plane-polygon map (critical operation to see which (new) walls are coplanar
			create_plane_polygon_map(pwalls_BSP);

			// Update blockables with new walls
			for (size_t i = 0; i < pwalls_BSP.size(); ++i)
				pwalls_BSP[i]->init_wall_state(&pwalls_BSP);

			update_blockable_walls(&pwalls_BSP);

			// Update direct_reflectables with the information from the plane polygon map 
			// E.g. for each (coplanar) wall which other walls are able to reflect this walls sources
			for (size_t i = 0; i < pwalls_BSP.size(); ++i)
				pwalls_BSP[i]->apply_plane_polygon_map(plane_polygon_map);

			for (auto& j : pwalls_BSP)
				walls_BSP.push_back(*j);

			// Make unify direct_reflectables for each plane polygon map entry and make them unique
			for (auto& i : plane_polygon_map) {
				std::vector<rts::wall*> direct_reflectables_united;
				for (auto& k : i) {
					for (auto& j : k->direct_reflectables) {
						direct_reflectables_united.push_back(j);
					}
				}
				// Make direct_reflectables for new walls unique and sort them in lexicographic order
				std::vector<rts::wall*>::iterator ip;
				std::sort(direct_reflectables_united.begin(), direct_reflectables_united.end());
				ip = std::unique(direct_reflectables_united.begin(), direct_reflectables_united.end());
				direct_reflectables_united.resize(std::distance(direct_reflectables_united.begin(), ip));
				i[0]->direct_reflectables = direct_reflectables_united;
			}

			// Make sure tree structure is somewhat well formed and find out tree height, used for visited nodes buffer vector size estimation in backtracking_with_BSP
			bsp_tree_height = traverseTree(bsp_tree);
			BOOST_LOG_TRIVIAL(info) << "Done building BSP tree, tree height: " << bsp_tree_height << std::endl;

#ifdef _DEBUG
			// Check BSP-tree && correctness of splitting algorithm
			write_obj("Object_" + std::to_string(1), pwalls_BSP);
			printBT(bsp_tree);
#endif // DEBUG

			return bsp_tree;
		}

		const BSPNode* build_BSP(std::vector<PolygonSpatial*> polygons, std::vector<rts::wall*>& new_walls, const double threshold) {
			// Check for convexity, then terminate, otherwise continue to build
			if (polygons.empty()) {
				return nullptr;
			}

			// Check convexity of polygons, if subspace spanned by polygons convex -> return 
			bool subspace_convex = true;
			for (auto i : polygons) {
				for (auto j : polygons) {
					if (i == j)
						continue;
					DEdge* edge_to_check = j->first();
					for (int k = 0; k < j->nPoints(); k++) {
						if (i->plane().whichSide(edge_to_check->srcPoint()) == BELOW) {
							subspace_convex = false;
							break;
						}
						edge_to_check = edge_to_check->next();
					}
					if (!subspace_convex)
						break;
				}
				if (!subspace_convex)
					break;
			}
			if (subspace_convex) {
				const BSPNode* node = new BSPNode{ construct_rtswall_model(polygons), nullptr, nullptr, true };
				for (auto j : node->node_walls)
					new_walls.push_back(j);
				return node;
			}
			// Check, which wall is best for splitting:
			// using r(p) && r(s) criterion Real-Time Processing of Image Sources Using Binary Space Partitioning on pg. 5/608
			// calculate r(p) with criterion of Ranta-Eskola
			//vector counting in (wall-id; behind; front; crosses/r(s)) for polygons
			std::vector<std::tuple<int, int, int, int>> measures;
			int temp_id = 0;
			for (auto& i : polygons) {
				measures.push_back(std::tuple<int, int, int, int>(temp_id, 0, 0, 0));
				for (auto& j : polygons) {
					if (i == j)
						continue;
					bool behind = false;
					bool in_front = false;
					DEdge* edge_to_check = j->first();
					for (int k = 0; k < j->nPoints(); k++) {
						if (i->plane().whichSide(edge_to_check->srcPoint()) == BELOW)
							behind = true;
						if (i->plane().whichSide(edge_to_check->srcPoint()) == ABOVE)
							in_front = true;
						edge_to_check = edge_to_check->next();
					}
					// if i in front of j -> infront += 1
					if (behind && !in_front)
						std::get<1>(measures.back()) += 1;
					// if i behind of j -> behind += 1
					if (!behind && in_front)
						std::get<2>(measures.back()) += 1;
					// if i crosses j -> crosses += 1
					if (behind && in_front)
						std::get<3>(measures.back()) += 1;
					// if i is ON j -> default case: in_front += 1
					if (!behind && !in_front)
						std::get<2>(measures.back()) += 1;
				}
				temp_id += 1;
			}
			// select splitting plane 
			std::vector<std::tuple<int, double, int>> r_p;
			for (auto& i : measures) {
				r_p.push_back(std::tuple<int, double, int>(std::get<0>(i), 0.0, 0));
				int behind = std::get<1>(i);
				int in_front = std::get<2>(i);
				if (behind == 0)
					std::get<1>(r_p.back()) = 0.0f;
				else if (in_front == 0)
					std::get<1>(r_p.back()) = 0.0f;
				else {
					double b = std::min((double)behind / (double)in_front, (double)in_front / (double)behind);
					std::get<1>(r_p.back()) = b;
				}
				std::get<2>(r_p.back()) = std::get<3>(i);
			}
			// (sort measures table by the required parameters)
			std::sort(r_p.begin(), r_p.end(), [](const std::tuple<int, double, int>& x, const std::tuple<int, double, int>& y)
				{
					return std::get<2>(x) < std::get<2>(y);
				});
			int wall_id = -1;
			for (auto& i : r_p) {
				// If number of crosses minimal and threshold criterion met -> use as partition wall
				if (std::get<1>(i) >= threshold) {
					// Found the wall id to use as partition node!
					wall_id = std::get<0>(i);
					break;
				}
			}
			// If no wall manages to satisfy the given Ranta-Eskola criterion threshold 
			// we default to to the element with the smallest difference to the criterion
			if (wall_id == -1) {
				for (auto& i : r_p) {
					std::get<1>(i) = abs(std::get<1>(i) - threshold);
				}
				std::sort(r_p.begin(), r_p.end(), [](const std::tuple<int, double, int>& x, const std::tuple<int, double, int>& y)
					{
						return std::get<1>(x) < std::get<1>(y);
					});
				wall_id = std::get<0>(r_p.front());
			}
			
			List<PolygonSpatial> above;
			List<PolygonSpatial> on;
			List<PolygonSpatial> below;

			PolygonSpatial* partition_wall = polygons[wall_id];

			for (int i = 0; i < polygons.size(); i++) {
				if (i == wall_id)
					continue;
				split(polygons[i], partition_wall->plane(), above, on, below);
			}
			// Collect walls, above, on and below the partition_wall to send off to construct the respective subtrees
			std::vector<PolygonSpatial*> above_polys, on_polys, below_polys;
			on_polys.push_back(partition_wall);
			forEachItemOnList(above)
				above_polys.push_back(getItem(PolygonSpatial));
			forEachItemOnList(on)
				on_polys.push_back(getItem(PolygonSpatial));
			forEachItemOnList(below)
				below_polys.push_back(getItem(PolygonSpatial));

			// Transmogrify the PolygonSpatial data structure back to rts::wall for use in the algorithm
			std::vector<rts::wall*> sortedWalls = construct_rtswall_model(on_polys);
			update_blockable_walls(&sortedWalls);
			const std::vector<rts::wall*> myWalls = sortedWalls;
			// Create binary tree node recursively
			const BSPNode* node = new BSPNode{ myWalls, build_BSP(above_polys, new_walls, threshold), build_BSP(below_polys, new_walls, threshold), false };

			for (auto j : node->node_walls)
				new_walls.push_back(j);
			// Initialise walls like normal
			for (size_t i = 0; i < node->node_walls.size(); ++i) {
				node->node_walls[i]->init_wall_state(&node->node_walls);
			}
			for (size_t i = 0; i < node->node_walls.size(); ++i) {
				node->node_walls[i]->sort_blockables_by_dist();
			}
			return node;
		}

		/*
		* disable walls marked in the config file; used only during initialization
		*/
		void disable_selected_walls(std::vector<rts::wall>& walls) {
			for (auto wall_id : walls_to_disable) {
				walls.at(wall_id).enabled = false;
			}
		}

		// Plane-Polygon Map as per Schroeder_Dirk_Diss_Physically_based_real_time_auralization.pdf (pg. 116)
		void create_plane_polygon_map(std::vector<rts::wall*> walls_to_check) {
			// Sort according to id for easier analysis down the line && backwards compatibility
			struct {
				bool operator()(rts::wall* a, rts::wall* b) const { return a->id < b->id; }
			} comparator;
			std::sort(walls_to_check.begin(), walls_to_check.end(), comparator);
			plane_polygon_map.push_back({ walls_to_check[0] });
			plane_polygon_map[0][0]->plane_polygon_map_id = 0;
			for (auto& one_wall : walls_to_check) {
				if (!one_wall->enabled) continue;
				for (int i = 0; i < plane_polygon_map.size(); i++) {
					if (one_wall->id == plane_polygon_map[i][0]->id)
						continue;
					bool same_normal = true;
					// Checking coplanarity by checking all points/corners of the polygon in question
					for (int k = 0; k < one_wall->n.size(); k++) {
						if ((plane_polygon_map[i][0]->n.at(k) != one_wall->n.at(k)) || (signbit(one_wall->n.at(k)) != signbit(plane_polygon_map[i][0]->n.at(k)))) {
							same_normal = false;
							break;
						}
					}
					// If distance from origin different, trivially non-coplanar
					if (one_wall->d != plane_polygon_map[i][0]->d)
						same_normal = false;
					// If distance and normals are coincident, and there is already a wall with these parameters add to already existing plane polygon map entry
					if (same_normal) {
						plane_polygon_map[i].push_back(one_wall);
						one_wall->plane_polygon_map_id = i;
						break;
					}
					// If at end of array and plane does not fit anywhere: make new entry in the plane-polygon map
					else if (i == plane_polygon_map.size() - 1) {
						plane_polygon_map.push_back({ one_wall });
						one_wall->plane_polygon_map_id = i + 1;
					}
				}
			}
			// Request to minimise memory footprint of plane_polygon_map
			for (auto& i : plane_polygon_map)
				i.shrink_to_fit();
			plane_polygon_map.shrink_to_fit();
		};

		// [...]

		// translating .obj data into polygonal model of the spatial partitioning code
		std::vector<PolygonSpatial*> construct_polygonspatial_model(std::vector<rts::wall> polygons) {
			std::vector<PolygonSpatial*> polygonVec;
			for (auto& i : polygons) {
				if (i.enabled) {
					std::vector<Point> pointVec;
					for (int j = 0; j < i.corners.size(); j++) {
						Point myPoint = Point(i.corners[j][0], i.corners[j][1], i.corners[j][2]);
						pointVec.push_back(myPoint);
					}
					polygonVec.push_back(new PolygonSpatial(pointVec, (int)i.id, i.n.at(0), i.n.at(1), i.n.at(2), i.d));
				}
			}
			return polygonVec;
		};

		std::vector<rts::wall*> construct_rtswall_model(std::vector<PolygonSpatial*> polygons) {
			std::vector<rts::wall*> wall_model;
			int id = 0;
			for (auto& i : polygons) {
				// Construct vector of arma:vecs
				std::vector<arma::fvec3> arma_vecs;
				DEdge* startEdge = i->first();
				arma_vecs.push_back(arma::fvec3{ startEdge->srcPoint().x(), startEdge->srcPoint().y(), startEdge->srcPoint().z() });
				DEdge* endEdge = i->first()->next();
				while (startEdge != endEdge) {
					arma_vecs.push_back(arma::fvec3{ endEdge->srcPoint().x(), endEdge->srcPoint().y(), endEdge->srcPoint().z() });
					endEdge = endEdge->next();
				}
				const unsigned int myID = (i->m_parentID * 10000 + id);
				// We force load these walls, since mostly numeric inaccuracies occur, which are not relevant to the validity of the geometry
				// Also, sometimes the usual wall calculation of the normal vectors are faulty, when edges are flipped/inserted by the spatial partitioning algorithm
				// In any case, we use the parent planes orientation, see below
				rts::wall* myWall = new rts::wall{ myID, arma_vecs, walls[i->m_parentID]->material, true, true };
				
				// The usual normal derivation algorithm by the spatial partitioning algorithm can be wonky, so we force to use the parents normal in all cases
				myWall->n = walls[i->m_parentID]->n;
				myWall->double_n = walls[i->m_parentID]->double_n;
				myWall->d = walls[i->m_parentID]->d;
				myWall->setParentID(i->m_parentID);
				myWall->setParent(walls[i->m_parentID]);
				wall_model.push_back(myWall);
				i++;
			}
			return wall_model;
		};

		int traverseTree(const BSPNode* node) {
			// Get the height of the tree
			if (!node)
				return 0;
			else {
				// Find the height of both subtrees
				// and use the larger one
				int left_height = traverseTree(node->front);
				int right_height = traverseTree(node->back);
				if (left_height >= right_height)
					return left_height + 1;
				else
					return right_height + 1;
			}
		}

// Some DEBUG functions, mainly to visualise space partitioning output
#ifdef _DEBUG
		void printBT(const std::string& prefix, const BSPNode* node, bool isLeft)
		{
			if (node != nullptr)
			{
				std::cout << prefix;

				std::cout << (isLeft ? "|---" : "'---");

				// print the value of the node
				for (auto i : node->node_walls)
					std::cout << i->id << " ";
				std::cout << std::endl;

				// enter the next tree level - left and right branch
				printBT(prefix + (isLeft ? "|   " : "    "), node->front, true);
				printBT(prefix + (isLeft ? "|   " : "    "), node->back, false);
			}
		}

		void printBT(const BSPNode* node)
		{
			printBT("", node, false);
		}

		static void printPolys(const char* const label, const List<PolygonSpatial>& pL)
		{
			if (pL.size()) {
				std::cout << "----------" << std::endl
					<< pL.size() << " polygon(s) " << label << std::endl;
				forEachItemOnList(pL) {
					std::cout << "  PolygonSpatial:" << std::endl;
					const PolygonSpatial* const g = getItem(PolygonSpatial);
					forEachDEdgeOfPoly(d1, g) {
						const Point& p = d1->srcPoint();
						std::cout << "  " << p.x() << ' ' << p.y() << ' ' << p.z() << std::endl;
					}
				}
			}
		}
#endif // DEBUG

	};
}
