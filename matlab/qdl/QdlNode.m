

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
        V0
        Vdc
        Va
        V1
        V2
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
            self.bnodes = double.empty(0);
            self.sbranches = double.empty(0);
            self.nbnode = 0;
            self.nsbranch = 0;
            
            self.source_type = QdlSystem.SourceNone;
            self.V0 = 0;
            self.Vdc = 0;
            self.Va = 0;
            self.V1 = 0;
            self.V2 = 0;
            
        end
        
        function add_bnode(self, gain, node)
            
            self.nbnode = self.nbnode + 1;
            self.T(self.nbnode) = gain;
            self.nbnodes(self.nbnode) = node;
            
        end
        
        function add_sbranch(self, gain, branch)
            
            self.nsbranch = self.nsbranch + 1;
            self.Z(self.nsbranch) = gain;
            self.sbranches(self.nsbranch) = branch;
            
        end
    
    end
 
end