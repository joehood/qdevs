

classdef QdlNode < QdlDevice
   
    properties
        
        C
        G
        H
        B
        S
        bnodes
        sbranches
        nbnode
        nsbranch
        v0
        vdc
        va
        v1
        v2
        source_type
        
    end
    
    methods
        
        function self = QdlNode(name, C, G, H)
            
            self@QdlDevice(name);
            
            self.C = C;
            self.G = G;
            self.H = H;
            
            self.B = double.empty(0);
            self.S = double.empty(0);
            self.bnodes = QdlNode.empty(0);
            self.sbranches = QdlBranch.empty(0);
            self.nbnode = 0;
            self.nsbranch = 0;
            
            self.source_type = QdlSystem.SourceNone;
            self.v0 = 0;
            self.vdc = 0;
            self.va = 0;
            self.v1 = 0;
            self.v2 = 0;
            
        end
        
        function add_bnode(self, node, gain)
            
            self.nbnode = self.nbnode + 1;
            self.B(self.nbnode) = gain;
            self.nbnodes(self.nbnode) = node;
            
        end
        
        function add_sbranch(self, branch, gain)
            
            self.nsbranch = self.nsbranch + 1;
            self.S(self.nsbranch) = gain;
            self.sbranches(self.nsbranch) = branch;
            
        end
    
    end
 
end