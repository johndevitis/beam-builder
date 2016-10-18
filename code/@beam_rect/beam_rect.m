classdef beam_rect < handle
%% classdef beam_rect
% 
% 
% 
% author: john devitis
% create date: 18-Oct-2016 14:55:17

%% object properties
	properties
        b = 1;          % section width [in]
        h = 12;         % section height [in]
        E = 29e6;       % modulus [psi] 
        ro = .29        % density [lb/in^3]
        grav = 386.4    % gravity [in/sec^2]
	end

%% dependent properties
	properties (Dependent)
        A % section area [in^2]
        I % moment of inertia (vertical) [in^4]
        mbar % self-weight [lbm/in] (m=f/a)
	end

%% private properties
	properties (Access = private)
	end

%% dynamic methods
	methods
	%% constructor
		function self = beam_rect()
		end

	%% ordinary methods


	%% dependent methods

		function A = get.A(self)
		%% area [in^2]
			A = self.b * self.h;
		end

		function I = get.I(self)
		%% moment of inertia (vertical) [in^4]
			I = self.b*self.h^3/12;
        end
        
        function mbar = get.mbar(self)
        %% self-weight [lbm/in] (m=f/a)
            mbar = self.A*self.ro; % lbf/in
            mbar = mbar/self.grav; % lbm/in
        end

	end

%% static methods
	methods (Static)
	end

%% protected methods
	methods (Access = protected)
	end

end
