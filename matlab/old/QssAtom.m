

classdef QssAtom < handle

    properties
        x      % internal state
        d      % derivative
        q      % external (quantized) state
        qlast  % ext state at last internal transition
        dQ     % quantum
        eps    % hysteresis width
        tlast  % time at last internal transition
        tnext  % next transition time (from ta or ext trigger)
        thist  % time history array
        qhist  % q history array
        k      % history arrayu index
    end  
  
    methods
  
        function obj = QssAtom(dQ, eps)
            obj.dQ = dQ;
            obj.eps = eps; 
        end
        
        function init(obj)
            obj.x = 0;
            obj.q = 0;
            obj.d = 0;
            obj.qlast = 0;
            obj.tlast = 0;
            obj.tnext = 0;
            obj.thist(1) = 0;
            obj.qhist(1) = 0;
            obj.k = 2;  
        end

        function [t, trigger] = update(obj, tstop)

            trigger = 0;
            
            % 1. delta_int function (update x)
            t = obj.tnext;
            dt = t - obj.tlast;
            obj.x = obj.x + obj.d * dt;  % integrate int state
            obj.tlast = t;
            
            % 2. lambda function (update q)
            if obj.x >= obj.q + obj.dQ - obj.eps
                obj.q = obj.q + obj.dQ;
            elseif obj.x <= obj.q - 0.5 * obj.dQ + obj.eps
                obj.q = obj.q - obj.dQ;
            end
            
            % 3. lambda function (tigger ext event if q changed):
            
            if obj.q ~= obj.qlast
            
                obj.qlast = obj.q;
                trigger = 1; % trigger ext event
                
                % save data:
                obj.thist(obj.k) = obj.tlast;
                obj.qhist(obj.k) = obj.q;
                obj.k = obj.k + 1;
                
            end

            % 3. ta function:
            
            if (obj.d > 0)
                obj.tnext = t + (obj.q + obj.dQ - obj.x) / obj.d;
            elseif (obj.d < 0)
                obj.tnext = t + (obj.q - 0.5 * obj.dQ - obj.x) / obj.d;
            else
                obj.tnext = inf;
            end
            
            obj.tnext = max(obj.tnext, 0);
            obj.tnext = min(obj.tnext, tstop);
            
        end
        
    end
    
end