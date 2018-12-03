
classdef LiqssLimRCI < LiqssLimDevice
  
    methods
  
        function self = LiqssLimRCI(name, R, C, I, v0, dQvoltage, dQcurrent, eps_factor)
            
            self@LiqssLimDevice(name, dQcurrent, dQvoltage, eps_factor);
            
            node = LiqssLimNode2(name, I, 1/R, C, v0, dQvoltage, eps_factor * dQvoltage);

            self.add_node(node);
            
        end

    end
    
end