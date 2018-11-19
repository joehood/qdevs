

classdef LiqssLimSystem2 < handle

    properties
        
        dQcurrent;  % quantization resolution
        dQvoltage;  % quantization resolution
        eps_factor; % hysterisis factor (eps = eps_factor * dQ)      
        devices;    % device models     
        ndevice;    % number of lim nodes
        time0;      % initial time
        time;       % current simulation time
        tlast;      % time at last update
        tstop;      % simulation stop time
        min_step;   % minimum time step
        
    end  
  
    methods
  
        function self = LiqssLimSystem2(dQvoltage, dQcurrent, eps_factor)
            
            self.dQvoltage = dQvoltage;
            self.dQcurrent = dQcurrent;           
            self.eps_factor = eps_factor;
            
            self.ndevice = 0;
            self.devices = LiqssLimDevice.empty(0);
            self.min_step = 1e-12;
            
        end
        
        function init(self, t0)
            
            self.time0 = t0;
            self.time = t0;
            self.tlast = t0;
            
            % initialize all devices:
            
            for k = 1:self.ndevice
                self.devices(k).init(self.time);
            end  
            
            % perform initial advance of all atoms:
            
            self.update_all(); 
                 
        end
        
        function add_device(self, device)
            
            self.ndevice = self.ndevice + 1;
            self.devices(self.ndevice) = device;
            
        end
        
        function run(self, tstop)
            
            self.tstop = tstop;
            
            % force initial update for this run period:
            
            self.update_all();

            while (self.time < self.tstop)
                
                self.advance();
             
            end
            
            % force final update and save for this run period:
            
            self.update_all();
            
        end
        
        function update_all(self)
            
            for k = 1:self.ndevice
                self.devices(k).update_branches(self.time);
            end
            
            for k = 1:self.ndevice
                self.devices(k).update_nodes(self.time);
            end  
            
        end 
        
        function update_next(self)
            
            force = self.time >= self.tstop;
            
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
            
        end
        
        function update_triggered(self)
            
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
            
        end
               
        function advance(self)     
            
            % get min tnext from all atoms:
            
            mintime = self.get_mintime();
            
            % advance the simulation time to the min tnext (or tstop):

            self.time = min(mintime, self.tstop);
            
            % perform next scheduled internal updates:

            self.update_next();
            
            % now update externally triggered atoms:
            
            self.update_triggered();
            
        end % function advance
              
        function [mintime] = get_mintime(self)
            
            mintime = inf;
            
            for j = 1:self.ndevice
                
                for k = 1:self.devices(j).nbranch
                    mintime = min(self.devices(j).branches(k).tnext, mintime);
                    self.devices(j).branches(k).trigger = 0;
                end

                for k = 1:self.devices(j).nnodes
                    mintime = min(self.devices(j).nodes(k).tnext, mintime);
                    self.devices(j).nodes(k).trigger = 0;
                end
                
            end
            
        end % function get_mintime

    end % methods
    
end % classdef QssLimSystem