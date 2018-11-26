
classdef QdlSystem < handle
    
    properties
        
        % scalar variables:
        
        nnode   % number of nodes
        nbranch % number of branchest 
        n       % nnode + nbranch (total number of states)
        npt     % size of output vectors
        dtmin   % minimum time step
        time    % current simualtion time
        tstop   % simulation stop time

        % liqss variables (dimension(n)):
        
        dq     % quantum level
        dqmin  % min quantum level
        dqmax  % max quantum level
        maxerr % nominal allowed error
        qlo    % lower boundary of quantized state
        qhi    % upper boundary of quantized state
        tlast  % last internal transition time of state
        tnext  % next predicted internal transition
        x      % internal states
        trig   % trigger flags for external transition
        iout   % output index
        source % const voltage value (node) or const current value (branch)
        
        % liqss variables with past values (n, 2):

        d % state derivatives
        q % quantized state

        % liqss output data arrays output (n, npt):
        
        qout % quantized states
        tout % output time points
        nupd % update count
        tupd % update count times
        iupd % update count array size

        % enhanced LIM model:
        
        A      % incidence matrix between node i and branch k (nnode, nbranch) 
        Cinv   % 1/capacitance at node i (nnode)
        G      % conductance at node i (nnode)
        H      % current injection at node i (nnode)
        B      % vccs gain matrix between node i and node p (nnode, nnode)
        S      % cccs gain matrix between node i and branch k (nnode, nbranch)
        Linv   % 1/inductance in branch k (nbranch)
        R      % series resistance in branch k (nbranch)
        E      % voltage at branch k from node i to node j (nbranch)
        T      % vcvs gain matrix between branch k and node i (nbranch, nnode)
        Z      % ccvs gain matrix between branch k and branch q (nbranch, nbranch)
        
        % switching:
        
        swatom % switching incidence matrix (n)
        swper  % switching period (n)
        swduty % switching duty (n)
        swon   % switching source value (n)
        swoff  % switching source value (n)

    end
    
    methods
        
        function self = QdlSystem(npt)
            
            self.npt     = npt;
            
            self.nnode   = 0;
            self.nbranch = 0;
              
            
            self.dtmin = 1e-12;  % hardcoded for now 

        end
        
        function build_lim(self)
            
            self.n = self.nnode + self.nbranch;
            
            self.A = zeros(self.nnode, self.nbranch);
            
            self.Cinv = zeros(self.nnode, 1);
            self.G    = zeros(self.nnode, 1);
            self.H    = zeros(self.nnode, 1);
            self.B    = zeros(self.nnode, self.nnode);
            self.S    = zeros(self.nnode, self.nbranch);

            self.Linv = zeros(self.nbranch, 1); 
            self.R    = zeros(self.nbranch, 1);
            self.E    = zeros(self.nbranch, 1);
            self.T    = zeros(self.nbranch, self.nnode);
            self.Z    = zeros(self.nbranch, self.nbranch);
            
            self.dq     = zeros(self.n, 1);
            self.dqmin  = zeros(self.n, 1);
            self.dqmax  = zeros(self.n, 1);
            self.maxerr = zeros(self.n, 1);
            
            self.source = zeros(self.n, 1);
            
            
        end
    
        function [idx] = add_node(self, node)
            
            self.nnode = self.nnode + 1;
            idx = self.nnode;
            
            self.Cinv(i) = 1/C;
            self.G(i) = G;
            self.H(i) = H;
            
            if p > 0
                self.B(i, p) = B;
            end
            
            if q > 0
                self.S(i, q) = S;
            end
            
            self.source(i) = source;
            
            self.dqmin(i) = dqmin;
            self.dqmax(i) = dqmax;
            self.maxerr(i) = maxerr;
            
            idx = self.nnode;
            
        end
        
        function [idx] = add_branch(self, k, i, j, L, R, E, p, T, q, Z, source, dqmin, dqmax, maxerr)
            
            % i(k)' = 1/L(k) * [e(k) + t(k,p) * v(p) + z(k,q) * i(q) + (v(i) - v(j)) - r(k) * i(k)]
            
            if i > 0
                self.A(i, k) = 1;
            end
            
            if j > 0
                self.A(j, k) = -1;
            end
            
            self.Linv(k) = 1/L;
            self.R(k) = R;
            self.E(k) = E;
            
            if p > 0
                self.T(k, p) = T;
            end
            
            if q > 0
                self.Z(k, q) = Z;
            end
            
            self.source(self.nnode+k) = source;
            
            self.dqmin(self.nnode+k) = dqmin;
            self.dqmax(self.nnode+k) = dqmax;
            self.maxerr(self.nnode+k) = maxerr;
            
            idx = self.nbranch;
            
        end
        
        function init(self)
            
            self.time = 0.0;
            
            self.qlo    = zeros(self.n, 1);
            self.qhi    = zeros(self.n, 1);     
            self.tlast  = zeros(self.n, 1);
            self.tnext  = zeros(self.n, 1);
            self.trig   = zeros(self.n, 1);
            self.x      = zeros(self.n, 1);
            
            self.d = zeros(self.n, 2);
            self.q = zeros(self.n, 2);           

            self.tnext(:) =  inf;
            self.dq(:)    =  self.dqmin(:);
            self.qhi(:)   =  self.dq(:); 
            self.qlo(:)   = -self.dq(:);

            % output data arrays:
            
            self.qout = zeros(self.n, self.npt);
            self.tout = zeros(self.n, self.npt);
            self.iout = zeros(self.n, 1);
            self.iout(:) = 1;
            
            % update counters:
            
            self.nupd = zeros(self.n, self.npt);
            self.tupd = zeros(self.n, self.npt);  
            self.iupd = zeros(self.n, 1);
            self.iupd(:) = 1;
            
            % init switching events:
            
            for iatom = 1:self.n
               if self.swatom(iatom)
                   self.source(iatom) = self.swon(iatom);
                   self.tnext(iatom) = self.swper(iatom) * self.swduty(iatom);
               end
            end
        end
        
        function runto(self, tstop)
            
            self.tstop = tstop;
            
            % force initial update and save for this run period:
            
            for iatom = 1:self.n
                self.update(iatom);
                self.save(iatom);
            end 
            
            % now advance the simulation until we reach tstop:
            
            while self.time < self.tstop
                self.advance();      
            end
            
            % force final update and save for this run period:
            
            self.time = self.tstop;
            
            for iatom = 1:self.n
                self.update(iatom);
                self.save(iatom);
            end 
        end
        
        function advance(self)
            
            % determine next time and advance time:
             
            self.time = min(min(self.tnext), self.tstop);
            self.time;

            % set force flag if we are at the simulation stop time:
            
            force = self.time >= self.tstop;
            
            % perform next scheduled internal updates:

            for iatom = 1:self.n
                if self.tnext(iatom) <= self.time || force
                   self.update(iatom);
                end
            end  
            
            % now update externally triggered atoms:
            
            for iatom = 1:self.n
                if self.trig(iatom)
                   self.update(iatom);
                end
            end  

        end
        
        function update(self, iatom)   
                      
            self.trig(iatom) = 0;  % reset trigger flag

            if self.source(iatom) > 0.0
                
                self.x(iatom) = self.source(iatom);
                
                if self.swatom(iatom)
                    tsw = mod(self.time, self.swper(iatom));
                    if tsw >= self.swper(iatom) * self.swduty(iatom)
                        self.source(iatom) = self.swoff(iatom);
                        self.tnext(iatom) = self.tlast(iatom) + self.swper(iatom);
                    else
                        self.source(iatom) = self.swon(iatom);
                        self.tnext(iatom) = self.tlast(iatom) + self.swper(iatom) * self.swduty(iatom);
                        self.tlast(iatom) = self.time;
                    end
                end
                
                self.q(iatom, 2) = self.x(iatom);
                self.d(iatom, :) = 0.0;
                    
            else
                self.dint(iatom);      % update internal state (delta_int function)
                self.quantize(iatom);  % quantize internal state
                self.d(iatom, 2) = self.f(iatom, self.q(iatom, 2));  % update derivative 
                self.ta(iatom);        % calculate new tnext
            end
            
            
            % trigger external update if quatized output changed:
            
            if self.q(iatom, 2) ~= self.q(iatom, 1) 
                self.q(iatom, 1) = self.q(iatom, 2);  % save last q
                self.save(iatom); % save output 
                self.dext(iatom); % set trigger flags
                self.update_dq(iatom);
            end
 
        end
        
        function [interp] = quantize(self, iatom)  
            
            % save old derivative:
            self.d(iatom, 1) = self.d(iatom, 2);
            
            interp = 0;
            change = 0;
            
            if self.x(iatom) >= self.qhi(iatom)
                self.q(iatom, 2) = self.qhi(iatom);
                self.qlo(iatom) = self.qlo(iatom) + self.dq(iatom);
                change = 1;
            elseif self.x(iatom) <= self.qlo(iatom)
                self.q(iatom, 2) = self.qlo(iatom);
                self.qlo(iatom) = self.qlo(iatom) - self.dq(iatom);
                change = 1;
            end
            
            self.qhi(iatom) = self.qlo(iatom) + 2 * self.dq(iatom);
            
            if change  % we've ventured out of qlo/self.qhi(iatom) bounds
            
                % calculate new derivative:
                self.d(iatom, 2) = self.f(iatom, self.q(iatom, 2));
                
                % if the derivative has changed signs, then we know 
                % we are in a potential oscillating situation, so
                % we will set the q such that the derivative=0:
                
                if self.d(iatom, 2) * self.d(iatom, 1) < 0  % if derivative has changed sign
                    flo = self.f(iatom, self.qlo(iatom)); 
                    fhi = self.f(iatom, self.qhi(iatom));
                    a = (2 * self.dq(iatom)) / (fhi - flo);
                    self.q(iatom, 2) = self.qhi(iatom) - a * fhi;
                    interp = 1;
                end
                
            end
            
        end
        
        function dint(self, iatom)
            
            self.x(iatom) = self.x(iatom) + self.d(iatom, 2) * (self.time - self.tlast(iatom));
            self.tlast(iatom) = self.time;
            
        end
        
        function ta(self, iatom)
            
            % estimate the time to the next self.qhi(iatom)/qlo crossing:
            
            if (self.d(iatom, 2) > 0)
                self.tnext(iatom) = self.time + (self.qhi(iatom) - self.x(iatom)) / self.d(iatom, 2);
                
            elseif (self.d(iatom, 2) < 0)
                self.tnext(iatom) = self.time + (self.qlo(iatom) - self.x(iatom)) / self.d(iatom, 2);
                
            else
                self.tnext(iatom) = inf;  % derivative == 0, so tnext is inf
                
            end
            
            % force delta t to be greater than 0:
            
            self.tnext(iatom) = max(self.tnext(iatom), self.tlast(iatom) + self.dtmin);
            
        end
        
        function save(self, iatom) 
            
            if self.time ~= self.tout(iatom, self.iout(iatom))
                self.iout(iatom) = self.iout(iatom) + 1;
                self.tout(iatom, self.iout(iatom)) = self.time;           
            end
            
            self.qout(iatom, self.iout(iatom)) = self.q(iatom, 2); 
            
            if self.time ~= self.tupd(self.iupd)
                self.iupd(iatom) = self.iupd(iatom) + 1;
                self.tupd(iatom, self.iupd(iatom)) = self.time;
            end
            
            self.nupd(iatom, self.iupd(iatom)) = self.nupd(iatom, self.iupd(iatom)) + 1;
             
        end
        
        function [d] = f(self, iatom, qval)
            
            if iatom <= self.nnode  % is node
                
                %  v(i)' = 1/c(i) * [h(i) + b(i,p) * v(p) + s(i,q) * i(q) + sum(ik) - g(i) * v(i)]
                
                isum = self.A(iatom, :) * self.q(self.nnode+1:end, 2);

                d = self.Cinv(iatom) * (self.H(iatom) + self.S(iatom,:) * self.x(self.nnode+1:end) ...
                    - qval * self.G(iatom) - isum);
                
            else  % is branch
                
                % i(k)' = 1/L(k) * [e(k) + t(k,p) * v(p) + z(k,q) * i(q) + (v(i) - v(j)) - r(k) * i(k)]
                
                ibranch = iatom-self.nnode;
                
                vij = self.A(:, ibranch)' * self.q(1:self.nnode, 2);
                
                d = self.Linv(ibranch) * (self.E(ibranch) + self.T(ibranch,:) * self.x(1:self.nnode) ...
                    - qval * self.R(ibranch) + vij);
                
            end
        
        end
        
        function dext(self, iatom)
            
           if iatom <= self.nnode
               
               %self.trig(self.nnode+1:end) = abs(self.A(iatom, :)) || abs(sign(self.S(iatom, :)));
               
               for ibranch = 1:self.nbranch
                   
                   self.trig(self.nnode+ibranch) = ...
                          self.A(iatom, ibranch) ~= 0 || ...
                          self.S(iatom, ibranch) ~= 0 || ...
                          self.T(ibranch, iatom) ~= 0;
                      
                   if self.T(ibranch, iatom) ~= 0
                      x=1;
                   end
               end
               
           else
               for inode = 1:self.nnode
                   
                   ibranch = iatom - self.nnode;
                   
                   self.trig(self.nnode) = ...
                          self.A(inode, ibranch) ~= 0 || ...
                          self.T(ibranch, inode) ~= 0 || ...
                          self.S(inode, ibranch) ~= 0;
                      
                  if self.S(inode, ibranch) ~= 0
                      x=1;
                  end
               end
           end

        end 
        
        function update_dq(self, iatom)
            
            if (self.dqmax(iatom) - self.dqmin(iatom)) < eps
                return
            end
            
            self.dq(iatom) = min(self.dqmax(iatom), max(self.dqmin(iatom), abs(self.maxerr(iatom) * self.q(iatom, 2)))); 
            
            self.qlo(iatom) = self.q(iatom) - self.dq(iatom);
            self.qhi(iatom) = self.q(iatom) + self.dq(iatom);
            
        end
            
        function [x, t] = run_ss(self, dt)
            
            a = zeros(self.n, self.n);
            b = zeros(self.n, self.n);
            
            for k=1:self.nnode
                
                a(k, k) = 1/C;
                
            end 
            
            t = 0:dt:self.tstop;
            npt = length(t);
            x = zeros(self.n, npt);
            
            for k=2:npt
                
            end
            
        end
            
    end

end