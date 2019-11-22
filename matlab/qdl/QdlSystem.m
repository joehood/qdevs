
classdef QdlSystem < handle  
    
    
    properties (Constant)
        
        StimNone    = 0;
        StimSource  = 1;
        StimSignal  = 2;
        
        SourceNone  = 0;
        SourceDC    = 1;
        SourceSINE  = 2;
        SourcePWM   = 3;
        
        SignalNone  = 0;
        SignalBGain = 1;
        SignalSGain = 2;
        SignalTGain = 3;
        SignalZGain = 4; 
        
        DefaultNpoints = 1e6;
        DefaultDtmin = 1e-12;
        
    end
   
    properties
        
        % scalar variables:
        
        nnode      % number of nodes
        nbranch    % number of branchest 
        nstim      % number of stimulus devices
        n          % nnode + nbranch (total number of lim atoms)
        ndevice    % nnode + nbranch + nstim (total number of qdl devices) 
        npt        % size of output vectors
        dtmin      % minimum time step
        time       % current simualtion time
        tstop      % simulation stop time
        def_dqmin  % default minimum dq
        def_dqmax  % default maximum dq
        def_dqerr  % default maximum error target
        
        % devs atom objects:
        
        nodes    % array of LIM nodes and ideal voltage sources (nnode)
        branches % array of LIM branches and ideal current sources (nbranch)
        stims    % array of stimuli (nstim)

        % liqss variables 
        
        dq     % quantum level (n)
        dqmin  % min quantum level (n)
        dqmax  % max quantum level (n)
        dqerr  % nominal allowed error (n)
        qlo    % lower boundary of quantized state (n)
        qhi    % upper boundary of quantized state (n)
        tlast  % last internal transition time of state (n)
        tnext  % next predicted internal transition (n)
        x      % internal states (n)
        x0     % initial internal states (n)
        trig   % trigger flags for external transition (n)
        iout   % output index for sparse output arrays (n)
        d      % state derivatives (n, 2)
        q      % quantized state (n, 2)
        M      % maps each atom to all receiver atoms for update broadcasting (n,n)

        % liqss output data arrays output :
        
        qout % quantized states (n, npt)
        tout % output time points (n, npt)
        nupd % update count (n, npt)
        tupd % update count times (n, npt)
        iupd % update count array size (n, npt)

        % enhanced LIM model:
        
        A  % incidence matrix between node i and branch k (nnode, nbranch)
        
        % data for latency nodes: 
        Cinv % 1/capacitance at node i (nnode)
        G    % conductance at node i (nnode)
        H    % current injection at node i (nnode)
        B    % vccs gain matrix between node i and node p (nnode, nnode)
        S    % cccs gain matrix between node i and branch k (nnode, nbranch)
        
        % data for latency branches:
        Linv % 1/inductance in branch k (nbranch)
        R    % series resistance in branch k (nbranch)
        E    % voltage at branch k from node i to node j (nbranch)
        T    % vcvs gain matrix between branch k and node i (nbranch, nnode)
        Z    % ccvs gain matrix between branch k and branch q (nbranch, nbranch)

        % stimulus params
        stim_type   % value from above stim type enum values (n)
        source_type % value from above source type enum values (n)
        signal_type % value from above signal type enum values (n)
        xdc         % offset value (n)
        xa          % sine amplitude (n)
        x1          % pwm first value (n)
        x2          % pwm second value (n)
        freq        % sine wave freqnecy or pwm switching freqency (Hz) (n)
        duty        % pwm duty cycle (pu) (n)
        phi         % sine wave phase (rad) (n)
        period      % cached 1/f (s) (n)
        
        % state space
        dtss % state space timestep
        tss  % state space current time
        xss  % state space state vector (for stepping method)
        Ass  % state space A matrix
        Bss  % state space B matix
        Uss  % state space U vector
        Apr  % state space descretixed A' matrix (backwards Euler)
        Bpr  % state space descretixed B' matrix (backwards Euler)
        
        timer

    end
    
    methods
        
        function self = QdlSystem(dqmin, dqmax, dqerr)
            
            self.def_dqmin = dqmin;
            self.def_dqmax = dqmax;
            self.def_dqerr = dqerr;
            self.npt      = self.DefaultNpoints;   
            self.nnode    = 0;
            self.nbranch  = 0;    
            self.nodes    = QdlNode.empty(0);
            self.branches = QdlBranch.empty(0);            
            self.dtmin    = self.DefaultDtmin;

        end
        
        function copyto(self, other)
            
            % scalar variables:
        
            other.nnode      = self.nnode;      
            other.nbranch    = self.nbranch;    
            other.n          = self.n;          
            other.npt        = self.npt;        
            other.dtmin      = self.dtmin;      
            other.time       = self.time;       
            other.tstop      = self.tstop;      
            other.def_dqmin  = self.def_dqmin;  
            other.def_dqmax  = self.def_dqmax;  
            other.def_dqerr  = self.def_dqerr;  
            
            % devs atom objects:
            
            other.nodes    = self.nodes;
            other.branches = self.branches;
            
            % liqss variables 
            
            other.dq     = self.dq;    
            other.dqmin  = self.dqmin;  
            other.dqmax  = self.dqmax;  
            other.dqerr  = self.dqerr;  
            other.qlo    = self.qlo;    
            other.qhi    = self.qhi;    
            other.tlast  = self.tlast;  
            other.tnext  = self.tnext;  
            other.x      = self.x;      
            other.x0     = self.x0;     
            other.trig   = self.trig;   
            other.iout   = self.iout;   
            other.d      = self.d;
            other.q      = self.q;
            other.M      = self.M;
            
            % enhanced LIM model:
            
            other.A    = self.A;
            
            % data for latency nodes: 
            other.Cinv = self.Cinv;
            other.G    = self.G;   
            other.H    = self.H;   
            other.B    = self.B;  
            other.S    = self.S;   
            
            % data for latency branches:
            other.Linv = self.Linv;
            other.R    = self.R;    
            other.E    = self.E;    
            other.T    = self.T;    
            other.Z    = self.Z;    
            
            % source stimulus params
            other.stim_type   = self.stim_type;
            other.source_type = self.source_type;
            other.signal_type = self.signal_type;
            other.xdc         = self.xdc;        
            other.xa          = self.xa;         
            other.x1          = self.x1;         
            other.x2          = self.x2;
            other.freq        = self.freq;       
            other.duty        = self.duty;       
            other.phi         = self.phi;        
            other.period      = self.period;     
            
            % state space
            other.dtss = self.dtss; 
            other.tss  = self.tss;  
            other.xss  = self.xss;  
            other.Ass  = self.Ass;  
            other.Bss  = self.Bss;  
            other.Uss  = self.Uss;  
            other.Apr  = self.Apr;  
            other.Bpr  = self.Bpr;  
            
        end
        
        function build_lim(self)
                        
            self.n = self.nnode + self.nbranch; % total number of atoms
            
            self.A = zeros(self.nnode, self.nbranch); % connectivity matrix
            
            self.M = zeros(self.n, self.n);  % trigger map

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
            
            self.x0    = zeros(self.n, 1);
            self.dq    = zeros(self.n, 1);
            self.dqmin = zeros(self.n, 1);
            self.dqmax = zeros(self.n, 1);
            self.dqerr = zeros(self.n, 1);
            
            self.xdc = zeros(self.n, 1);
            self.xa  = zeros(self.n, 1);
            self.x1  = zeros(self.n, 1);
            self.x2  = zeros(self.n, 1);

            self.stim_type   = zeros(self.n, 1);
            self.source_type = zeros(self.n, 1);
            self.signal_type = zeros(self.n, 1);
            self.freq        = zeros(self.n, 1);
            self.period      = zeros(self.n, 1);
            self.duty        = zeros(self.n, 1);
            self.phi         = zeros(self.n, 1);

            for inode = 1:self.nnode
                
                self.Cinv(inode) = 1 / self.nodes(inode).C;
                self.G(inode) = self.nodes(inode).G;
                self.H(inode) = self.nodes(inode).H; 
                
                self.x0(inode) = self.nodes(inode).v0;
                self.dqmin(inode) = self.nodes(inode).dqmin;
                self.dqmax(inode) = self.nodes(inode).dqmax;
                self.dqerr(inode) = self.nodes(inode).dqerr;
                
                self.stim_type(inode)   = self.nodes(inode).stim_type;
                self.source_type(inode) = self.nodes(inode).source_type;
                self.signal_type(inode) = self.nodes(inode).signal_type;
                self.xdc(inode)         = self.nodes(inode).vdc;
                self.xa(inode)          = self.nodes(inode).va;
                self.x1(inode)          = self.nodes(inode).v1;
                self.x2(inode)          = self.nodes(inode).v2;
                self.freq(inode)        = self.nodes(inode).freq;
                self.period(inode)      = 1/self.nodes(inode).freq;
                self.duty(inode)        = self.nodes(inode).duty;
                self.phi(inode)         = self.nodes(inode).phi;
                
                % add entries to B matrix:
                for ibnodes = 1:self.nodes(inode).nbnode
                    self.B(inode, self.nodes(inode).bnodes(ibnodes).index) ...
                        = self.nodes(inode).B(ibnodes);
                end
                
                % add entries to S matrix:
                for isbranch = 1:self.nodes(inode).nsbranch
                    self.S(inode, self.nodes(inode).sbranches(isbranch).bindex) ...
                        = self.nodes(inode).S(isbranch);
                end

            end 
            
            for ibranch = 1:self.nbranch
                
                iatom = self.nnode + ibranch;
                
                self.A(self.branches(ibranch).inode.index, ibranch) = 1; 
                self.A(self.branches(ibranch).jnode.index, ibranch) = -1; 
 
                self.Linv(ibranch) = 1 / self.branches(ibranch).L;
                self.R(ibranch) = self.branches(ibranch).R;
                self.E(ibranch) = self.branches(ibranch).E;
                
                self.x0(iatom) = self.branches(ibranch).i0;
                self.dqmin(iatom) = self.branches(ibranch).dqmin;
                self.dqmax(iatom) = self.branches(ibranch).dqmax;
                self.dqerr(iatom) = self.branches(ibranch).dqerr;
                
                self.stim_type(iatom)   = self.branches(ibranch).stim_type;
                self.source_type(iatom) = self.branches(ibranch).source_type;
                self.signal_type(iatom) = self.branches(ibranch).signal_type;
                self.xdc(iatom)         = self.branches(ibranch).idc;
                self.xa(iatom)          = self.branches(ibranch).ia;
                self.x1(iatom)          = self.branches(ibranch).i1;
                self.x2(iatom)          = self.branches(ibranch).i2;
                self.freq(iatom)        = self.branches(ibranch).freq;
                self.period(iatom)      = 1/self.branches(ibranch).freq;
                self.duty(iatom)        = self.branches(ibranch).duty;
                self.phi(iatom)         = self.branches(ibranch).phi;

                % add entries to T matrix:
                for itnodes = 1:self.branches(ibranch).ntnode
                    self.T(ibranch, self.branches(ibranch).tnodes(itnodes).index) ...
                        = self.branches(ibranch).T(itnodes);
                end
                
                % add entries to Z matrix:
                for izbranch = 1:self.branches(ibranch).nzbranch
                    self.Z(ibranch, self.branches(ibranch).zbranches(izbranch).bindex) ...
                        = self.branches(ibranch).Z(izbranch);
                end
                
            end
            
        end
        
        function build_map(self)
            
            for inode = 1:self.nnode
                
                % add connections from to B nodes:
                for ibnodes = 1:self.nodes(inode).nbnode
                    if self.nodes(inode).bnodes(ibnodes).source_type == self.StimNone
                        self.M(self.nodes(inode).bnodes(ibnodes).index, self.nodes(inode).index) = 1;
                    end
                end
                
                % add connections from to S branches:
                for isbranch = 1:self.nodes(inode).nsbranch
                    if self.nodes(inode).sbranches(isbranch).source_type == self.SourceNone
                        self.M(self.nodes(inode).sbranches(isbranch).index, self.nodes(inode).index) = 1;
                    end
                end
                
            end
            
            for ibranch = 1:self.nbranch
                
                iatom = self.nnode + ibranch;
                
                if self.branches(ibranch).inode.source_type == self.SourceNone
                    self.M(iatom, self.branches(ibranch).inode.index) = 1;
                end
                
                if self.branches(ibranch).jnode.source_type == self.SourceNone
                    self.M(iatom, self.branches(ibranch).jnode.index) = 1;
                end
                
                if self.branches(ibranch).source_type == self.SourceNone
                    self.M(self.branches(ibranch).inode.index, iatom) = 1;
                    self.M(self.branches(ibranch).jnode.index, iatom) = 1;
                end
                
                % add connections from T nodes:
                for itnodes = 1:self.branches(ibranch).ntnode
                    if self.branches(ibranch).tnodes(itnodes).source_type == self.SourceNone
                        self.M(self.branches(ibranch).tnodes(itnodes).index, self.branches(ibranch).index) = 1;
                    end
                end
                
                % add connections from Z branches:
                for izbranch = 1:self.branches(ibranch).nzbranch
                    if self.branches(ibranch).zbranches(izbranch).source_type == self.SourceNone
                        self.M(self.branches(ibranch).zbranches(izbranch).index, self.branches(ibranch).index) = 1;
                    end
                end
                
            end
            
        end
    
        function add_node(self, node)

            if node.dqmin <= 0
                node.dqmin = self.def_dqmin;
            end
            
            if node.dqmax <= 0
                node.dqmax = self.def_dqmax;
            end
            
            if node.dqerr <= 0
                node.dqerr = self.def_dqerr;
            end
            
            self.nnode = self.nnode + 1;
            node.index = self.nnode;
            self.nodes(self.nnode) = node;
            
            % increment branch indeces as these are offest by the nodes:
            
            for kbranch = 1:self.nbranch
                self.branches(kbranch).index = self.branches(kbranch).index + 1;
            end

        end
        
        function add_branch(self, branch)
            
            if branch.dqmin <= 0
                branch.dqmin = self.def_dqmin;
            end
            
            if branch.dqmax <= 0
                branch.dqmax = self.def_dqmax;
            end
            
            if branch.dqerr <= 0
                branch.dqerr = self.def_dqerr;
            end
            
            self.nbranch = self.nbranch + 1;
            branch.bindex = self.nbranch; % not offset by nnode
            branch.index = self.nbranch + self.nnode;  % offest by nnode
            self.branches(self.nbranch) = branch;
            
        end
        
        function init(self)
            
            self.build_lim();
            
            self.build_map();
            
            self.build_qdl();

            self.time = 0.0;
  
        end
        
        function build_qdl(self)
            
            % dimension liqss arrays:
            
            self.qlo    = zeros(self.n, 1);
            self.qhi    = zeros(self.n, 1);     
            self.tlast  = zeros(self.n, 1);
            self.tnext  = zeros(self.n, 1);
            self.trig   = zeros(self.n, 1);
            self.x      = zeros(self.n, 1);
            self.d      = zeros(self.n, 2);
            self.q      = zeros(self.n, 2);  

            % initialize liqss arrays:
            
            self.tnext(:) =  inf;
            self.dq(:)    =  self.dqmin(:);
            self.qhi(:)   =  self.dq(:); 
            self.qlo(:)   = -self.dq(:);
            self.x(:)     =  self.x0(:);
            self.q(:, 1)  =  self.x0(:);
            self.q(:, 2)  =  self.x0(:);

            % output data arrays:
            
            self.qout = zeros(self.n, self.npt);
            self.tout = zeros(self.n, self.npt);
            self.iout = zeros(self.n, 1);
            self.iout(:) = 1;
            
            % set up counters:
            
            self.nupd = zeros(self.n, self.npt);
            self.tupd = zeros(self.n, self.npt);  
            self.iupd = zeros(self.n, 1);
            self.iupd(:) = 1;
            
        end 
        
        function runto(self, tstop)
            
            self.tstop = tstop;
            
            % force initial update and save for this run period:
            
            for iatom = 1:self.n
                self.internal_update(iatom);
                self.save(iatom);
            end 
            
            % now update intially externally triggered atoms:
            
            while any(self.trig)
                for iatom = 1:self.n
                   if self.trig(iatom)
                       self.internal_update(iatom);
                   end
                end 
            end
            
            % now advance the simulation until we reach tstop:
            
            %self.timer = timer();
            %self.timer.ExecutionMode = 'fixedRate';
            %self.timer.Period = 1;
            
            %self.timer.TimerFcn = {@(~,~,obj, x) disp(['sim time: ' ...
            %    num2str(round(obj.time)) ' of ' x]), self, num2str(tstop)};
            
            %start(self.timer);
            
            tdisp = 0;
            tint = 10;
            
            while self.time < self.tstop
                self.advance();   
                if self.time - tdisp > tint
                    disp(['sim time = ' num2str(round(self.time))]);
                    tdisp = self.time;
                end
            end
            
            %stop(self.timer);
            
            % force final update and save for this run period:
            
            self.time = self.tstop;
            
            for iatom = 1:self.n
                self.internal_update(iatom);
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
                   self.internal_update(iatom);
                end
            end  
            
            % now update externally triggered atoms:
            
            while any(self.trig)
                for iatom = 1:self.n
                    if self.trig(iatom)
                       self.internal_update(iatom);
                    end
                end 
            end

        end
        
        function internal_update(self, iatom)   
            
            self.trig(iatom) = 0;  % reset trigger flag
            
            self.dint(iatom);      % update internal state (delta_int function)
            self.quantize(iatom);  % quantize internal state
            self.d(iatom, 2) = self.f(iatom, self.q(iatom, 2));  % update derivative 
            self.ta(iatom);        % calculate new tnext

            % trigger external update if quatized output changed:
            
            if self.q(iatom, 2) ~= self.q(iatom, 1)
                self.save(iatom); % save output 
                self.q(iatom, 1) = self.q(iatom, 2);  % save last q
                self.trigger(iatom); % set trigger flags
                self.update_dq(iatom);
            end
 
        end
        
        function external_update(self, iatom)
           
            self.trig(iatom) = 0;  % reset trigger flag
            
            self.dext(iatom);      % update internal state
            self.quantize(iatom);  % quantize internal state
            self.d(iatom, 2) = self.f(iatom, self.q(iatom, 2));  % update derivative 
            self.ta(iatom);        % calculate new tnext

            % trigger external update if quatized output changed:
            
            if self.q(iatom, 2) ~= self.q(iatom, 1)
                self.save(iatom); % save output 
                self.q(iatom, 1) = self.q(iatom, 2);  % save last q
                self.trigger(iatom); % set trigger flags replace with tnext = time + dtmin for all trigger items?
                self.update_dq(iatom);
            end
            
        end
        
        function dint(self, iatom)
            
            if self.source_type(iatom) == self.SourceNone
                
                self.x(iatom) = self.x(iatom) + self.d(iatom, 2) * (self.time - self.tlast(iatom));
            
            elseif self.source_type(iatom) == self.SourceDC
                
                self.x(iatom) = self.xdc(iatom);
                
            elseif self.source_type(iatom) == self.SourcePWM
                
                if self.duty(iatom) == 0
                    self.x(iatom) = self.x1(iatom);
                    
                elseif self.duty(iatom) == 1
                    self.x(iatom) = self.x2(iatom);

                else

                    w = mod(self.time, self.period(iatom));

                    if self.isnear(w, self.period(iatom))
                        self.x(iatom) = self.x1(iatom);
                        
                    elseif self.isnear(w, self.period(iatom) * self.duty(iatom))
                        self.x(iatom) = self.x2(iatom);
                        
                    elseif w < self.period(iatom) * self.duty(iatom)
                        self.x(iatom) = self.x1(iatom);
                        
                    else
                        self.x(iatom) = self.x2(iatom);
                        
                    end

                end
            
            elseif self.source_type(iatom) == self.SourceSINE
                
                self.x(iatom) = self.xdc(iatom) + self.xa(iatom) * ...
                    sin(2*pi*self.freq(iatom)*self.time + self.phi(iatom));
                
            end
            
            self.tlast(iatom) = self.time;
            
        end
        
        function interp = quantize(self, iatom)  
            
            if self.source_type(iatom) == self.SourceNone
            
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
                
            else
                
                self.q(iatom, 2) = self.x(iatom);  % all sources fix q=x
                
            end
            
        end
        
        function ta(self, iatom)
            
            if self.source_type(iatom) == self.SourceNone
                
                if (self.d(iatom, 2) > 0)
                    self.tnext(iatom) = self.time + (self.qhi(iatom) - self.x(iatom)) / self.d(iatom, 2);
                elseif (self.d(iatom, 2) < 0)
                    self.tnext(iatom) = self.time + (self.qlo(iatom) - self.x(iatom)) / self.d(iatom, 2);
                else
                    self.tnext(iatom) = inf;  % derivative == 0, so tnext is inf
                end
             
            elseif self.source_type(iatom) == self.SourceDC
                
                self.tnext(iatom) = inf;
                
            elseif self.source_type(iatom) == self.SourcePWM
                
                if self.duty(iatom) == 0 || self.duty(iatom) == 1
                    self.tnext(iatom) = inf;
                    
                else

                    w = mod(self.time, self.period(iatom));

                    if self.isnear(w, self.period(iatom))
                        self.tnext(iatom) = self.time + self.period(iatom) * self.duty(iatom);
                    elseif self.isnear(w, self.period(iatom) * self.duty(iatom))
                        self.tnext(iatom) = self.time + self.period(iatom) - w; 
                    elseif w < self.period(iatom) * self.duty(iatom)
                        self.tnext(iatom) = self.time + self.period(iatom) * self.duty(iatom) - w;
                    else
                        self.tnext(iatom) = self.time + self.period(iatom) - w;
                    end

                end
            
            elseif self.source_type(iatom) == self.SourceSINE
                
                [q, self.tnext(iatom)] = self.update_sine(self.xdc(iatom), ...
                    self.xa(iatom), self.freq(iatom), self.phi(iatom), ...
                    self.time, self.dq(iatom));
                
            end
            
            % force delta t to be greater than 0:
            
            self.tnext(iatom) = max(self.tnext(iatom), self.tlast(iatom) + self.dtmin);
            
        end
        
        function save(self, iatom) 
            
            if 0  % enable for zoh recording
                self.iout(iatom) = self.iout(iatom) + 1;
                self.tout(iatom, self.iout(iatom)) = self.time - self.dtmin/2;
                self.qout(iatom, self.iout(iatom)) = self.q(iatom, 1);
            end
            
            %if ~self.isnear(self.time, self.tout(iatom, self.iout(iatom)))
            %if abs(self.time - self.tout(iatom, self.iout(iatom))) >= self.dtmin
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
        
        function d = f(self, iatom, qval)
            
            d = 0;
            
            if self.source_type(iatom) == self.SourceNone  % only for latent atoms

                if iatom <= self.nnode  % node

                    %  v(i)' = 1/c(i) * [h(i) + b(i,p) * v(p) + s(i,q) * i(q) + sum(ik) - g(i) * v(i)]

                    isum = self.A(iatom, :) * self.q(self.nnode+1:end, 2);

                    % d = 1/C * [ H + B * qnode + S * qbranch - q * G - isum ] 
                    
                    d = self.Cinv(iatom) * (self.H(iatom) + self.B(iatom,:) * self.q(1:self.nnode, 2) ...
                        + self.S(iatom,:) * self.q(self.nnode+1:end, 2) - qval * self.G(iatom) - isum);

                else  % branch

                    % i(k)' = 1/L(k) * [e(k) + t(k,p) * v(p) + z(k,q) * i(q) + (v(i) - v(j)) - r(k) * i(k)]

                    ibranch = iatom-self.nnode;

                    vij = self.A(:, ibranch)' * self.q(1:self.nnode, 2);
                    
                    % d = 1/L * [ E + T * qnode + Z * qbranch - q * R - vij ] 

                    d = self.Linv(ibranch) * (self.E(ibranch) + self.T(ibranch,:) * self.q(1:self.nnode, 2) ...
                        + self.Z(ibranch,:) * self.q(self.nnode+1:end, 2) - qval * self.R(ibranch) + vij);

                end
                
            elseif self.source_type(iatom) == self.SourceSINE
                
                d1 = 2 * pi * self.freq(iatom) * self.xa(iatom) * ...
                    cos(2 * pi * self.freq(iatom) * self.time + self.phi(iatom));
                
                d2 = -4 * pi * pi * self.freq(iatom) * self.freq(iatom) * self.xa(iatom) * ...
                    sin(2 * pi * self.freq(iatom) * self.time + self.phi(iatom));
                
                if abs(d1) > abs(d2)
                    d = d1;
                else
                    d = d2;
                end
                
            end
        
        end
        
        function dext(self, iatom)
            
            if self.source_type(iatom) == self.SourceNone
                
                self.x(iatom) = self.x(iatom) + f(self, iatom, self.x(iatom)) * (self.time - self.tlast(iatom));
            
            end
            
            self.tlast(iatom) = self.time;

        end 
        
        function trigger(self, iatom)
            
            for jatom = 1:self.n
                if self.M(iatom, jatom)
                    self.trig(jatom) = 1;
                end
            end
            
        end
        
        function update_dq(self, iatom)
            
            if (self.dqmax(iatom) - self.dqmin(iatom)) < eps
                return
            end
            
            self.dq(iatom) = min(self.dqmax(iatom), max(self.dqmin(iatom), ...
                                 abs(self.dqerr(iatom) * self.q(iatom, 2)))); 
            
            self.qlo(iatom) = self.q(iatom) - self.dq(iatom);
            self.qhi(iatom) = self.q(iatom) + self.dq(iatom);
            
        end
            
        function build_ss(self)
            
            n = self.n;
            nn = self.nnode;
            nb = self.nbranch;

            Ann = diag(self.Cinv) * (self.B - diag(self.G));
            Anb = diag(self.Cinv) * (self.S - self.A);
            Abn = diag(self.Linv) * (self.T + self.A');
            Abb = diag(self.Linv) * (self.Z - diag(self.R));
            
            self.Ass = [Ann, Anb; Abn, Abb];
            
            self.Bss = diag(cat(1, self.Cinv, self.Linv));
            
            self.Uss = [self.H; self.E];
          
        end
        
        function [t, x] = run_ss_to(self, dt, tstop)
            
            self.build_ss();
            
            self.dtss = dt;
            
            t = self.time : dt : tstop;
            npt = length(t);
            
            x = zeros(self.n, npt);
            
            x(:, 1) = self.x;
            
            self.Apr = inv(eye(self.n) - dt * self.Ass);
            self.Bpr = self.Apr * self.Bss * dt;
            
            for k = 2:npt
                
                x(:,k) = self.Apr * x(:,k-1) + self.Bpr * self.Uss;
                
            end
            
        end
        
        function init_ss(self, dt)
            
            self.build_ss();
            
            self.dtss = dt;
            
            if self.tss == 0
                self.xss = self.x0(:);
            end
            
            self.Apr = inv(eye(self.n) - dt * self.Ass);
            self.Bpr = self.Apr * self.Bss * dt;
        end
        
        function step_ss(self)
            
            self.tss = self.tss + self.dtss;
            self.xss = self.Apr * self.xss + self.Bpr * self.Uss;
            
        end
        
        function plot(self, atom, dots, lines, upd, cumm_upd, bins, xlbl, tss, xss, ymax)
            
            k = atom.index;
            
            ss = length(tss) > 1;
            
            if upd
                yyaxis left
            end
            
            if ss
                plot(tss, xss, 'c--', ...
                    'DisplayName','ss'); hold on;
            end
            
            if dots
                plot(self.tout(k,1:1:self.iout(k)), ...
                    self.qout(k,1:1:self.iout(k)), 'k.', ...
                    'DisplayName','qss'); hold on;
            end
            
            if lines
                plot(self.tout(k,1:self.iout(k)), ...
                    self.qout(k,1:self.iout(k)), 'b-', ...
                    'DisplayName', 'qss'); hold on;
            end

            ylabel('atom state');
            
            if upd
                yyaxis right
                
                updstr = strcat('updates per: ', num2str(self.time/bins), ' s');
                
                if cumm_upd
                    plot(self.tupd(k,1:1:self.iupd(k)), cumsum(self.nupd(k,1:1:self.iupd(k)))/2, 'r-',...
                        'DisplayName', 'updates');
                else
                    histogram(self.tupd(k,1:1:self.iupd(k)), bins, 'EdgeColor', 'none', 'FaceAlpha', 0.3, ...
                        'DisplayName', updstr);
                end
                title(atom.name);
                ylabel('updates');
                if ymax > 0
                    ylim([0, ymax]);
                end
            end
            
            xlim([-1.0, self.time]);
            
            loc = 'southeast';

            %leg = legend();
            %leg.Location = loc;
            
            if xlbl
               xlabel('t (s)');
            end
            
        end
         
    end
    
    methods(Static)
        
        function result = isnear(a, b)
            
            result = abs(a-b) < 1e4*eps(min(abs(a), abs(b)));
            
        end
        
        function [q, tnext] = update_sine(x0, xa, f, phi, t, dq)

            T = 1/f;                 % period
            w = mod(t, T);           % cycle time
            t0 = t - w;              % cycle base time
            omega = 2*pi*f;          % angular velocity
            theta = omega*w + phi;   % wrapped angular position
            x = xa * sin(2*pi*f*t);  % magnitude shifted by x0

            % determine next transition time. Saturate at +/- xa:
            
            if theta < pi/2        % quadrant I
                tnext = t0 + (asin(min(1, (x + dq)/xa))) / omega;
            elseif theta < pi      % quadrant II
                tnext = t0 + T/2 - (asin(max(0, (x - dq)/xa))) / omega;
            elseif theta < 3*pi/2  % quadrant III
                tnext = t0 + T/2 - (asin(max(-1, (x - dq)/xa))) / omega;
            else                   % quadrant IV
                tnext = t0 + T + (asin(min(0, (x + dq)/xa))) / omega;
            end
            
            % update state for current time:
            q = x0 + xa * sin(2*pi*f*t + phi);

        end
        
    end

end