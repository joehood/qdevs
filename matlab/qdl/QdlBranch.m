

classdef QdlBranch < QdlDevice
   
    properties
        
        inode
        jnode
        bindex
        L
        R
        E
        T
        Z
        tnodes
        zbranches
        ntnode
        nzbranch
        i0
        idc
        ia
        i1
        i2
        stim_type
        source_type
        signal_type
        
    end
    
    methods
        
        function self = QdlBranch(name, L, R, E)
            
            self@QdlDevice(name);

            self.L = L;
            self.R = R;
            self.E = E;
            
            self.T = double.empty(0);
            self.Z = double.empty(0);
            self.tnodes = QdlNode.empty(0);
            self.zbranches = QdlBranch.empty(0);
            self.ntnode = 0;
            self.nzbranch = 0;
            
            self.stim_type = QdlSystem.StimNone;
            self.source_type = QdlSystem.SourceNone;
            self.signal_type = QdlSystem.SignalNone;
            self.i0 = 0;
            self.idc = 0;
            self.ia = 0;
            self.i1 = 0;
            self.i2 = 0;

        end
        
        function connect(self, inode, jnode)
            
            self.inode = inode;
            self.jnode = jnode;
            
        end
        
        function add_tnode(self, node, gain)
            
            self.ntnode = self.ntnode + 1;
            self.T(self.ntnode) = gain;
            self.tnodes(self.ntnode) = node;
            
        end
        
        function add_zbranch(self, branch, gain)
            
            self.nzbranch = self.nzbranch + 1;
            self.Z(self.nzbranch) = gain;
            self.zbranches(self.nzbranch) = branch;
            
        end
    
    end
 
end