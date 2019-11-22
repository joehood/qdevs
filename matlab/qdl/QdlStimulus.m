

classdef QdlStimulus < QdlDevice
   
    properties
        
        connected_nodes
        connected_branches
        nnode
        nbranch
        x0
        xdc
        xa
        x1
        x2
        source_type
        
    end
    
    methods
        
        function self = QdlStimulus(name)
            
            self@QdlDevice(name);
            
            self.connected_nodes = QdlNode.empty(0);
            self.connected_branches = QdlBranch.empty(0);
            self.nnode = 0;
            self.nbranch = 0;
            
            self.source_type = QdlSystem.SourceNone;
            self.x0 = 0;
            self.xdc = 0;
            self.xa = 0;
            self.x1 = 0;
            self.x2 = 0;
            
        end
        
        function add_node(self, node)
            
            self.nnode = self.nnode + 1;
            self.connected_nodes(self.nnode) = node;
            
        end
        
        function add_branch(self, branch)
            
            self.nbranch = self.nsbranch + 1;
            self.connected_branches(self.nbranch) = branch;
            
        end
        
        function add_bgain(self, branch)
            
            self.nbranch = self.nsbranch + 1;
            self.connected_branches(self.nbranch) = branch;
            
        end
    
    end
 
end