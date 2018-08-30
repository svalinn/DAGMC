#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/GeomTopoTool.hpp"

class Network {
   // Link class for use in the network
   class Link {
      public:
      // make new link
      Link(moab::EntityHandle A, moab::EntityHandle B) {
          eA = A;
          eB = B;
      }
      // get the entities of the link
      void GetLink(moab::EntityHandle &A, moab::EntityHandle &B) {
        A = eA;
        B = eB;
      }
      private:
      moab::EntityHandle eA;
      moab::EntityHandle eB;  
   };

   // convenience link type
   typedef std::pair<moab::EntityHandle,moab::EntityHandle> link_t ;

   public:
    Network(){}
   ~Network(){}

   public:
   // Add link to the network
   void AddLink(moab::EntityHandle a, moab::EntityHandle b) {
     link_t link(a,b);
     network[link] = new Link(a,b); // add link to network
   }

   // remove link from the network
   void RemoveLink(moab::EntityHandle a, moab::EntityHandle b) {
       link_t link(a,b);
       network.erase(link); // remove from the map
   }

   // find the nodes that are adjacent to the current entity
   std::vector<moab::EntityHandle> FindNeighbours(moab::EntityHandle parent) {
     std::vector<moab::EntityHandle> neighbours;
     std::map<link_t,Link*>::iterator link_it = network.begin();
     for ( ; link_it != network.end() ; link_it++ )
       if ( (link_it->first).first == parent) neighbours.push_back((link_it->first).second);
     return neighbours;   
   }

   moab::EntityHandle FindParent(moab::EntityHandle child) {
    std::map<link_t,Link*>::iterator link_it = network.begin();
    for ( ; link_it != network.end() ; link_it++ )
      if ( (link_it->first).second == child) return (link_it->first).first;
    return 0;
  }
  
   // assuming that the top node is indeed the top node
   // navigate the network creating parent child links
   // (i.e. delting child->parent links)
   void Directionalise(moab::EntityHandle top) {

       std::vector<moab::EntityHandle> top_volumes; // volumes at the current level
       std::vector<moab::EntityHandle> children; // the children of the current volume
       
       // loop over the top volumes
       std::vector<moab::EntityHandle> new_top;
       new_top.push_back(top); // populate with initial value
       
       // iterators
       std::vector<moab::EntityHandle>::iterator child_it;
       std::vector<moab::EntityHandle>::iterator vec_it;

       while (new_top.size() != 0) { // while there are nodes to process
         top_volumes = new_top; // set the top_volumes to be the children of the last round
         new_top.clear(); // clear out history
         for ( vec_it = top_volumes.begin() ; vec_it != top_volumes.end() ; ++vec_it ) {
            children = FindNeighbours(*vec_it); // get children 
            for ( child_it = children.begin() ; child_it != children.end() ; ++child_it ) {
              RemoveLink(*child_it,*vec_it); // remove backwards link
            }
            // new level of links to remove
            new_top.insert(new_top.end(), children.begin(), children.end());
         }
       } // all done
   }

   // Print the Network
   void Print() {
      std::cout << "digraph {" << std::endl;
      std::map<std::pair<moab::EntityHandle, moab::EntityHandle>, Link* >::iterator map_it;
      for ( map_it = network.begin() ; map_it != network.end() ; ++map_it ) {
          std::cout << (map_it->first).first; 
          std::cout << "->";
          std::cout << (map_it->first).second;
          std::cout << std::endl;
      }
      std::cout << "}" << std::endl;
   }

   private:
   // the network;
   std::map<std::pair<moab::EntityHandle, moab::EntityHandle>, Link* > network;
};

typedef moab::Range::iterator range_it;

class DetermineHierarchy {
    public:
    DetermineHierarchy(moab::Interface* moab, moab::GeomTopoTool *geomtool = NULL);
   ~DetermineHierarchy();

   public:
   // Determine Hierarchy for the whole instance
   moab::ErrorCode DetermineTheHierarchy(bool loaded = false);
   // determine child volumes of a given volume
   moab::ErrorCode DetermineVolume(moab::EntityHandle volume);
   // Print moab range
   void PrintRange(moab::Range item);
   // Print Tree 
   void PrintTree(moab::EntityHandle vol, moab::Range item);
   void PrintTree();
   // Get the children of the the current parent
   moab::ErrorCode GetChildren(moab::EntityHandle parent,
			       std::vector<moab::EntityHandle> &children);

   moab::ErrorCode GetChildren(int parent_id,
			       std::vector<int> &children_id);

  moab::ErrorCode FindParent(moab::EntityHandle child, moab::EntityHandle &parent);
  moab::ErrorCode FindParent(int child, int &parent);


   private:
   // Given a volume and surface set, return the other unique volumnes
   // that share the surface
   moab::ErrorCode GetParentSets(moab::EntityHandle volume, moab::EntityHandle surface, 
                   moab::Range &parent_volumes);

   private:
   moab::Interface *mbi;
   moab::GeomTopoTool *gtt;
   Network *network;
};
