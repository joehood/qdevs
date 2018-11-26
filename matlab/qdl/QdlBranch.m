

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
            
            self.source_type = QdlSystem.SourceNone;
            self.I0 = 0;
            self.Idc = 0;
            self.Ia = 0;
            self.I1 = 0;
            self.I2 = 0;

        end
    
    end
 
end