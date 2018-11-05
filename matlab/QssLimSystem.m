

classdef QssLimSystem < handle

    properties
        
        branches; % atomic qss branch models
        nodes;    % atomic qss node models
        nnode;    % number of lim nodes
        nbranch;  % number of lim branches
        time;     % current simulation time
        tstop;    % simulation stop time
        
    end  
  
    methods
  
        function obj = QssLimSystem()
            
            obj.nnode = 0;
            obj.nbranch = 0;
            obj.branches = QssLimBranch.empty(0);
            obj.nodes = QssLimNode.empty(0);
            
        end
        
        function init(obj, t0)
            
            obj.time = t0;
            
            % initialize state, quantized, history, and derivatives:
            
            for k = 1, obj.nbranch;
                obj.branches(k).init(obj.time);
            end  
            
            for k = 1, obj.nnode;
                obj.nodes(k).init(obj.time);
            end 
            
            % perform initial advance of all atoms:
            
            for k = 1, obj.nbranch;
                obj.branches(k).advance(obj.time);
            end  
            
            for k = 1, obj.nnode;
                obj.nodes(k).advance(obj.time);
            end 
            
        end
        
        function add_node(obj, node)
            
            obj.nnode = obj.nnode + 1;
            obj.nodes(obj.nnode) = node;
            
        end
        
        function add_branch(obj, branch)
            
            obj.nbranch = obj.nbranch + 1;
            obj.branches(obj.nbranch) = branch;
            
        end
        
        function run(obj, tstop)
            
            obj.tstop = tstop;
            
            while (obj.time > obj.tstop)
                
                obj.advance();
             
            end
            
        end
               
        function advance(obj)
            
            
            
        end

    end
    
end