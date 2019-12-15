
classdef QdlGenSet < handle 
    
    properties
        
        R
        C
        L
        G
        B
        J
        
    end
    
    methods
        
        function self = QdlGenSet(R, C, L, G, B, J)    
            self.R = R;
            self.C = C;
            self.L = L;
            self.G = G;
            self.B = B;
            self.J = J;
        end
            
    end
    
end