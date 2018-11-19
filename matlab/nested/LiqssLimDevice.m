

classdef LiqssLimDevice < handle

    properties
        
        name;       % device name
        dQcurrent;  % quantization resolution
        dQvoltage;  % quantization resolution
        eps_factor; % hysterisis factor (eps = eps_factor * dQ)
        branches;   % atomic qss branch models
        nodes;      % atomic qss node models
        ground;     % ground node
        nnode;      % number of lim nodes
        nbranch;    % number of lim branches
        time;       % current simulation time
        tlast;      % last update time
        tstop;      % simulation stop time
        min_step;   % minimum time step
        
    end  
  
    methods
  
        function self = LiqssLimDevice(name, dQvoltage, dQcurrent, eps_factor)
            
            self.name = name;
            self.dQvoltage = dQvoltage;
            self.dQcurrent = dQcurrent;
            self.eps_factor = eps_factor;
            
            self.nnode = 0;
            self.nbranch = 0;
            self.branches = LiqssLimBranch.empty(0);
            self.nodes = LiqssLimNode2.empty(0);
            self.min_step = 1e-12;
            
        end
        
        function init(self, t0)
            
            self.time = t0;
            self.tlast = t0;
            
            % initialize state, quantized, history, and derivatives:
            
            for k = 1:self.nbranch
                self.branches(k).init(self.time);
            end  
            
            for k = 1:self.nnode
                self.nodes(k).init(self.time);
            end 
            
            % perform initial advance of all atoms:
            
            for k = 1:self.nbranch
                self.branches(k).update(self.time);
            end  
            
            for k = 1:self.nnode
                self.nodes(k).update(self.time);
            end 
            
        end
        
        function add_node(self, node)
            
            self.nnode = self.nnode + 1;
            self.nodes(self.nnode) = node;
            
        end
        
        function add_branch(self, branch)
            
            self.nbranch = self.nbranch + 1;
            self.branches(self.nbranch) = branch;
            
        end
               
        function advance(self)
            
            mintime = inf;
            
            for k = 1:self.nbranch
                mintime = min(self.branches(k).tnext, mintime);
                self.branches(k).trigger = 0;
            end
            
            for k = 1:self.nnode
                mintime = min(self.nodes(k).tnext, mintime);
                self.nodes(k).trigger = 0;
            end
                    
            self.time = min(mintime, self.tstop);
            
            force = self.time >= self.tstop;
            
            % perform next scheduled internal updates:

            for k = 1:self.nbranch
                if self.branches(k).tnext <= self.time || force
                   self.branches(k).update(self.time);
                end
            end  
            
            for k = 1:self.nnode
                if self.nodes(k).tnext <= self.time || force
                   self.nodes(k).update(self.time);
                end
            end 
            
            % now update externally triggered atoms:
            
            for k = 1:self.nbranch
                if self.branches(k).trigger
                   self.branches(k).update(self.time);
                end
            end  
            
            for k = 1:self.nnode
                if self.nodes(k).trigger
                   self.nodes(k).update(self.time);
                end
            end 
            
        end % function advance

    end % methods
    
end % classdef QssLimSystem