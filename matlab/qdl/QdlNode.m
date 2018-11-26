

classdef QdlNode < QdlDevice
   
    properties
        
        C
        G
        H
        B
        S
        bnodes
        sbranches
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
            
            self.source_type = QdlSystem.SourceNone;
            self.V0 = 0;
            self.Vdc = 0;
            self.Va = 0;
            self.V1 = 0;
            self.V2 = 0;
            
        end
    
    end
 
end