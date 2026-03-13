#include "preprocess.h"
#include "make_edge_conformal.hpp"
#include <vector>
#include <iostream>
#include <array>
#include <cassert>
#include <algorithm>
#define VERBOSE 0



void my_assert(bool fine, std::string message = std::string()){
  if(!fine){
    std::cout << message << std::endl;
    throw std::runtime_error {message};
  }
}

struct less_than_key {
  inline bool operator() (const std::array<int,2>& edge1,const std::array<int,2>& edge2){
    std::array<int,2> sedge1 = edge1;
    std::array<int,2> sedge2 = edge2;
    std::sort(sedge1.begin(),sedge1.end());
    std::sort(sedge2.begin(),sedge2.end());
    if(sedge1[0] < sedge2[0]){
      return true;
    }else if(sedge1[0]==sedge2[0]){
      if(sedge1[1]< sedge2[1]){
	return true;
      }else{
	return false;
      }
    }else{
      return false;
    }
    my_assert(false,"compare");
  }
};
struct oposite {
  inline bool operator() (const std::array<int,2>& edge1,const std::array<int,2>& edge2) const{
    if((edge1[0] == edge2[1]) && (edge1[1] == edge2[0])){
      return true;
    }else{
      return false;
    }
  }
};


bool next_edge(const std::array<int,2>& left,const std::array<int,2>& right){
	    return left[1] == right[0];
};

std::vector<int> sorted_outer_boundary(const struct processed_grid& grid,
			   const std::vector<std::vector<int>>& dir_faces,
			   const std::vector<std::vector<int>>& dir_hfaces,
			   int dir,
			   int cell){    
  std::vector<std::array<int,2>> edges;
for(int locind = 0; locind < dir_faces[dir].size(); locind++){
	int hface = dir_hfaces[dir][locind];
	int face = grid.cell_faces[hface];
        std::vector<std::array<int,2>> face_edges;
	// add edges
	for(int fnodePos = grid.face_ptr[face]; fnodePos < grid.face_ptr[face+1]-1; ++fnodePos){
	  std::array<int,2> new_edge ={grid.face_nodes[fnodePos],grid.face_nodes[fnodePos+1]};
	  face_edges.push_back( new_edge);
	}
	std::array<int,2> last_edge ={ grid.face_nodes[grid.face_ptr[face+1]-1], grid.face_nodes[grid.face_ptr[face]]};
	
	face_edges.push_back( last_edge);
	if(cell == grid.face_neighbors[2*face+ 1]){
	  // reverte edges
	  for(std::array<int,2>& edge: face_edges){
	    int tmp = edge[0];
	    //face_nodes[face].push_back(tmp);// just fill face_nodes her should be same as in processed grid;
	    edge[0] = edge[1];
	    edge[1] = tmp;
	  }
	}
	edges.insert(edges.end(),face_edges.begin(), face_edges.end());
      }
            
      // sort so intenal edges is after each other
      std::sort(edges.begin(),edges.end(),less_than_key());
#if VERBOSE
      std::cout << " Edges " << std::endl;
      for(const auto& e: edges){
	std::cout << e[0] << " " << e[1] << std::endl;
      }
#endif
      // // auto ip = std::unique(edges.begin(),edges.end(),oposite());
      // // edges.resize(std::distance(edges.begin(), ip));
      // oposite is_oposite;
      // std::vector<bool> remove_index(edges.size(),false);
      
      // for(int i=0; i < edges_size()-1; ++i){
      // 	auto e = edges[i];
      // 	auto en = edges[i+1];
      // 	if(is_oposite(e,en)){
      // 	  remove_index[i] = true;
      // 	  remove_index[i+1] = true;
      // 	}	
      // }

      std::vector<std::array<int,2 >> new_edges;
      // remove internal edges
      auto iter = edges.begin();
      auto iternext = iter;
      ++iternext;
      oposite is_oposite;
      for(; iternext != edges.end(); ){
	if(is_oposite(*iter,*iternext)){
	  iter = iternext;
	  ++iter;
	  if(iter != edges.end()){
	    iternext = iter;
	    ++iternext;
	  }else{
	    iternext = iter;
	  }
	}else{
	  new_edges.push_back(*iter);
	  ++iter;
	  ++iternext;
	}
      }
      if(iter !=edges.end()){
	new_edges.push_back(*iter);
      }
      // sort so intenal edges is after each other
      edges = new_edges;
#if VERBOSE    
      std::cout << " Edges outer" << std::endl;
      for(const auto& e: edges){
	std::cout << e[0] << " " << e[1] << std::endl;
      }
#endif
      // now edges should only contain oriented outer edges
      // Now put them after each other
      std::vector<std::array<int,2>> sedges;// = edges;
      sedges.push_back(*edges.begin());  
      edges.erase(edges.begin());
      while(edges.size()>0){
	auto cedges = sedges.back();//std::prev(sedges.end());
	auto nextelem = std::find_if(edges.begin(),
				  edges.end(), [cedges](std::array<int,2> &s){return s[0] == cedges[1];});
	sedges.push_back(*nextelem);
	edges.erase(nextelem);
      // 	std::cout << " *** " <<std::endl;
      //  for(const auto& e: edges){
      // 	std::cout << e[0] << " " << e[1] << std::endl;
      // }
      }
#if VERBOSE         
      std::cout << " Edges ordered" << std::endl;
      for(const auto& e: sedges){
	std::cout << e[0] << " " << e[1] << std::endl;
      }
#endif
     std::vector<int> sedge;
     for(const auto& edge: sedges){
	  sedge.push_back(edge[0]);
     }
     return sedge;
}

std::vector<int> new_tb(const struct processed_grid& grid,
			const std::vector<int>& bfnodes_in,
			const std::vector<int>& sedge, std::array<int,2>& bedge,
			int fsign,
			int fsigntb){
	std::vector<int> newedge;
	std::vector<int> bfnodes = bfnodes_in;
	auto ind2 = find_if(sedge.begin(),sedge.end(), [bedge](const int &s){return s == bedge[1];});
	auto ind1 = find_if(sedge.begin(),sedge.end(), [bedge](const int &s){return s == bedge[0];});
	int addnode = 0;
	my_assert(ind1 != sedge.end(),"ind1 not found");
	my_assert(ind2 != sedge.end(),"ind2 not found");  
	if( (ind1 == (std::prev(sedge.end()))) && (ind2 == sedge.begin()) ){
	  //std::cout << "end to first" << std::endl;
	    addnode = 0;
	}else if( (ind2 - ind1) < 0 ){
	  //std::cout << "Cyclic case nodes " << std::endl;
	    newedge.insert(newedge.end(),ind1,sedge.end());
	    newedge.insert(newedge.end(),sedge.begin(),std::next(ind2));
	    addnode = newedge.size() -2;
	}else if ( (ind2 - ind1) > 0) {
	  //std::cout << "Ordered case nodes " << std::endl;
	    newedge.insert(newedge.end(), ind1, std::next(ind2));
	    addnode = newedge.size() -2;
	}else{
	  //std::cout << "Do not exist" << std::endl;
	  my_assert(false,"Do not exist");
	}
#if VERBOSE    	
	std::cout << " addnode " << addnode << " ";
	std::cout << " newedge ";
	for(const int& n: newedge){
	  std::cout << n << " ";
	}
#endif
	if(fsigntb == 1){
              std::reverse(newedge.begin(),newedge.end());
        }
#if VERBOSE
	std::cout << std::endl << " newedge ordered" << std::endl;
	for(const int& n: newedge){
	  std::cout << n << " ";
	}

	std::cout << std::endl;
	std::cout << "BT Face nodes ";
	for(const int& n: bfnodes){
	  std::cout << n << " ";
	}
	std::cout << std::endl;
#endif
	  // ett possibly modified nodes	  
	  if( addnode > 0 ){
	    // if(fsign){
	    //   std::reverse(newedge.begin(),newedge.end());
	    // }
	    // iterator to start and end of new nodes
	    
	    auto iterstart = newedge.begin();
	    ++iterstart;
	    auto iterend = newedge.end();
	    --iterend;
	    auto node2 = std::find_if(bfnodes.begin(), bfnodes.end(),
				      [newedge](int s){return newedge.back() == s;});
	    auto node1 = std::find_if(bfnodes.begin(), bfnodes.end(),
				      [newedge](int s){return *(newedge.begin()) == s;});
	    if(node1 == std::prev(bfnodes.end()) &&  node2 == bfnodes.begin()){
	      
	      bfnodes.insert(bfnodes.end(),iterstart,iterend);
	    }else if( node2 - node1 == 1){
	      bfnodes.insert(std::next(node1),iterstart, iterend);
	    }else{
	      // nodes should al readdy be added
	      auto node_end = node2;
	      auto it1 = node1;
	      auto it2 = newedge.begin();
	      while(it1 != node_end && it1 != bfnodes.end()){
		//std::cout << *it1 << " " << *it2 << std::endl;
	       	my_assert(*it1 == *it2,"Existing face wrong");
	       	++it1;++it2;
	      }
	      if(it1 == bfnodes.end() && node2 != std::prev(bfnodes.end())){
		it1 = bfnodes.begin();
		while(it1 != node_end && it1 != bfnodes.end()){
		  //std::cout << *it1 << " " << *it2 << std::endl;
		  my_assert(*it1 == *it2,"Existing face wrong cyclic case");
		  ++it1;++it2;
		}
	      }
	      addnode = 0;
	    }
	  }
          return bfnodes;
}

void fix_edges_at_top(const struct processed_grid& grid,
                      std::vector<int>& nodes,
                      std::vector<int>& nodePos){    

  // are going to be the new face nodes
  std::vector<std::vector<int>> face_nodes(grid.number_of_faces);
  for(int i=0; i < grid.number_of_faces; ++i){
    for(int fpos = grid.face_ptr[i]; fpos < grid.face_ptr[i+1]; fpos++){
      face_nodes[i].push_back(grid.face_nodes[fpos]);
    }
  }
  int nhf=grid.cell_facePos[grid.number_of_cells];
  for(int cell=0; cell < grid.number_of_cells; ++cell){
#if VERBOSE    
    std::cout<< "*********************" << std::endl;
    std::cout<< "Processing cell " << cell << std::endl;
#endif
    // process a cells
    // and top and bottom face number
    std::vector<std::vector<int>> dir_faces(6);
    std::vector<std::vector<int>> dir_hfaces(6);

   for(int hface = grid.cell_facePos[cell]; hface < grid.cell_facePos[cell+1]; ++hface){
	int face = grid.cell_faces[hface];
	int hface_tag = grid.cell_faces[hface+nhf];
	dir_faces[hface_tag].push_back(face);
        dir_hfaces[hface_tag].push_back(hface);
   }
#if VERBOSE       
   for(int dir = 0; dir < 6; ++dir){
     std::cout << "Face size dir " <<  dir << " size " <<  dir_faces[dir].size()  << std::endl;
   }
#endif
   my_assert(dir_faces[4].size() == 1,"face size wrong bottom");
   my_assert(dir_faces[5].size() == 1,"face size wrong top");

    
    for(int dir = 0; dir < 4; ++dir){
 #if VERBOSE         
      std::cout<< "*********************" << std::endl;
      std::cout << "Processing dir " << dir << std::endl;
      std::cout << "Number faces " << dir_faces[dir].size() << std::endl;
 #endif
      if(dir_faces[dir].size()>1){// can not be intersection with only one face
      // process each vertical direction
      // find all oriented edges in each direction
      std::vector<int> sedge = sorted_outer_boundary(grid,
						     dir_faces,
						     dir_hfaces,
						     dir,
						     cell);
#if VERBOSE
      std::cout << "Sedge ";
      for(const int& e: sedge){
	std::cout << e << " ";
      }
      std::cout << std::endl;
#endif      
      // now sedge should be orient orderd list of edges
      // find top/bottom edge to be considered
      std::array<int,2> bedge;	
      // loop over top and bootm
      for(int tb=4; tb<6; ++tb){
#if VERBOSE    	
	std::cout<< "********" << std::endl;
	std::cout << "Look at " << tb << std::endl;
#endif
	my_assert(dir_faces[tb].size() == 1,"face size wrong tb/bottom");
	int bface = dir_faces[tb][0];
	std::vector<int> org_bfnodes;
	for(int fnodePos = grid.face_ptr[bface]; fnodePos < grid.face_ptr[bface+1]; ++fnodePos){
	  org_bfnodes.push_back(grid.face_nodes[fnodePos]);
	}

	std::vector<int> bfnodes=face_nodes[bface];//bd
#if VERBOSE    	
	std::cout << "BT Face nodes ";
	for(const int& n: bfnodes){
	  std::cout << n << " ";
	}
	std::cout << std::endl;
	std::cout << "BT Org nodes ";
	for(const int& n: org_bfnodes){
	  std::cout << n << " ";
	}
	std::cout << std::endl;
#endif
	std::array<int,4> odir = {0,2,1,3};
	if(odir[dir] == 0){
	  bedge[0] = org_bfnodes[3];
	  bedge[1] = org_bfnodes[0];
	}else{
	  bedge[0] = org_bfnodes[odir[dir]-1];
	  bedge[1] = org_bfnodes[odir[dir]];
	}
	int fsign = 1;
	int fsigntb = 1;
	if(grid.face_neighbors[2*bface+1] == cell){
          fsigntb=-1;                    
        }	

        if(fsigntb == 1){
              std::reverse(bedge.begin(),bedge.end());
        }
#if VERBOSE    	
	std::cout << "bedge " << bedge[0] << " " << bedge[1]
		  << " fsign " << fsign << " fbsigntb " << fsigntb << std::endl;
#endif    
       auto bfnodes_new = new_tb(grid,
				 bfnodes,
				 sedge,
				 bedge,
				 fsign,
				 fsigntb);
#if VERBOSE    	       
       std::cout << "new face nodes ";
	for(const int& n: bfnodes_new){
	  std::cout << n << " ";
	}
	std::cout << std::endl;
#endif
       face_nodes[bface] = bfnodes_new;
      }
      }//end tb
    }// end dir
  }// end cell
  // add the others
  for(size_t face=0; face < face_nodes.size(); ++face){
    if(face_nodes[face].size() == 0){
      for(int fnodePos = grid.face_ptr[face]; fnodePos < grid.face_ptr[face+1]; ++fnodePos){
	face_nodes[face].push_back(grid.face_nodes[fnodePos]);
      }
    }
  }
  
  nodePos.resize(face_nodes.size()+1);
  nodePos[0]=0;
  int count = 0;
  for(const auto& face: face_nodes){
    nodes.insert(nodes.end(),face.begin(),face.end());
    nodePos[count+1]=nodePos[count]+face.size();
    ++count;
  }
}//end cell

extern "C" void make_edge_conformal(struct processed_grid* grid){
    std::cout << "Fixing edge grid to be edge conformal" << std::endl;
    std::vector<int> nodes;
    std::vector<int> nodePos;
    fix_edges_at_top(*grid, nodes, nodePos);
    int nf = grid->number_of_faces;
    for (int i = 0; i <= nf; i++) {
      grid->face_ptr[i] = nodePos[i];
    }
    int nfn = grid->face_ptr[nf];  /* Total number of face nodes */
    for (int i = 0; i < nfn; i++) {
      grid->face_nodes[i] = nodes[i];
    }
}
