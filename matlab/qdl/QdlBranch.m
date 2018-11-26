

classdef QdlBranch < QdlDevice
   
    properties
        
        inode
        jnode
        L
        R
        E
        T
        Z
        tnodes
        zbranches
        ntnode
        nzbranch
        I0
        Idc
        Ia
        I1
        I2
        source_type
        
    end
    
    methods
        
        function self = QdlBranch(name, inode, jnode, L, R, E)
            
            self@QdlDevice(name);
            
            self.inode = inode;
            self.jnode = jnode;
            self.L = L;
            self.R = R;
            self.E = E;
            
            self.T = double.empty(0);
            self.Z = double.empty(0);
            self.tnodes = double.empty(0);
            self.zbranches = double.empty(0);
            self.ntnode = 0;
            self.nzbranch = 0;
            
            self.source_type = QdlSystem.SourceNone;
            self.I0 = 0;
            self.Idc = 0;
            self.Ia = 0;
            self.I1 = 0;
            self.I2 = 0;

        end
        
        function add_tnode(self, gain, node)
            
            self.ntnode = self.ntnode + 1;
            self.T(self.ntnode) = gain;
            self.tnodes(self.ntnode) = node;
            
        end
        
        function add_zbranch(self, gain, branch)
            
            self.nzbranch = self.nzbranch + 1;
            self.Z(self.nzbranch) = gain;
            self.zbranches(self.nzbranch) = branch;
            
        end
    
    end
 
end