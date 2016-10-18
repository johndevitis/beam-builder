function show(self)
%% disp
% 
% 
% 
% author: john devitis
% create date: 18-Oct-2016 17:03:11
	
	
%% Summary
% System State Information: 
fprintf('Mass: \n');
disp(self.M); 

fprintf('Stiffness: \n');
disp(self.K); 

fprintf('Natural Frequencies [Hz]: \n');
disp(self.F); 

fprintf('Damped natural frequencies [rad/sec]: \n');
disp(self.Wn); 

fprintf('Modal Mass: \n');
disp(self.Mr);

fprintf('Modal Stiffness: \n');
disp(self.Kr);

% fprintf('Complex Roots:\n')
% fprintf('Pole %d\t Damping Factor: %6.3f\t Positive Pole: %6.3f\n',...
%     [1:length(root); real(root)'; imag(root)'])
	
	
end
