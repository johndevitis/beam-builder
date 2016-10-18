classdef beam_builder < handle
%% classdef beam_builder
% 
% 
% 
% author: john devitis
% create date: 18-Oct-2016 14:06:42 (refactor)
%   * originally created 7/25/13, and heavily modified 082015-112015

%% object properties
	properties
        L % total beam length [in]
        nn % number of dof (# of inner nodes excluding end bc's)
        mchk = 1 % 1 = lumped mass, else = continuous
        bc % boundary condition node indices
        k % element stiffness matrix
        m % element mass matrix
        K % global stiffness matrix
        M % global mass matrix      
	end

%% dependent properties
	properties (Dependent)
        nel % number of beam elements
        l % element length [in]
        x % 
    end

%% private properties
	properties (Access = private)
	end

%% dynamic methods
	methods
	%% constructor
		function self = beam_builder()
        end
        
        
	%% dependent methods
        function x = get.x(self)
            % xdimension
            x = linspace(0,self.L,self.nn);
        end
    
        function nel = get.nel(self)
		%% number of beam elements
			nel = self.nn-1;
        end

		function l = get.l(self)
		%% element length [in]
			l = self.L/self.nel;
        end
       
    end
%% static methods
	methods (Static)
        bh = example2;
        bh = example5;
        [Kr,Mr] = beam_condense(K,M,rm,mchk);
	end

%% protected methods
	methods (Access = protected)
	end

end
