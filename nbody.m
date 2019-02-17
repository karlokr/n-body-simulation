function [t, r, v, m, E, T, V] = nbody(tmax, level, r0, v0, m0, tracefreq)
%  nbody Solves the n-body problem using an O(deltat^2) FDA
%  
%  Input arguments
%
%      tmax:      (real scalar) Final solution time.
%      level:     (integer scalar) Discretization level.
%      r0:        (nb*3 matrix) Initial positions of the n bodies in
%	          Cartesian coordinates, nb is the number of bodies
%      v0:        (nb*3 matrix) Initial velocities of the n bodies in
%		  Cartesian coordinates, nb is the number of bodies
%      m0:	  (1*nb array) Masses of the n bodies, nb is the number
%	          of bodies
%      tracefreq: (optional integer scalar) Frequency of tracing output, 
%                 0 disables tracing.
%
%  Output arguments
%
%      t:      (real vector) Vector of length nt = 2^level + 1 containing
%              discrete times (time mesh).
%      r:      (nb*3*nt matrix) Matrix of dimension nb*3*nt, where nt is the
%	       total number of timesteps and nb is the number of bodies, 
%	       containing computed values of the positions of the n bodies
%	       at discrete times t(n). 
%      v:      (nb*3*nt matrix) Matrix of dimension nb*3*nt, where nt is the
%	       total number of timesteps and nb is the number of bodies,
%	       containing computed values of the positions of the n bodies 
%	       at discrete times t(n).
%      m:      (real vector) Vector of length nb, where nb is the number of
%	       bodies, containing the masses of the n bodies
%      E:      (real vector) Vector of length nt = 2^level + 1 containing
%	       the total energy at discrete times t(n).
%      T:      (real vector) Vector of length nt = 2^level + 1 containing 
%	       the total kinetic energy at discrete times t(n)
%      V:      (real vector) Vector of length nt = 2^level + 1 containing 
%	       the toal potential energy at discrete times t(n).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	% Tracing control: if 6th arg is supplied base tracing on that input,
	% otherwise use local defaults.
	if nargin > 5
		if tracefreq == 0
			trace = 0;
		else
			trace = 1;
		end
	else
		trace = 1;
		tracefreq = 100;
	end

	if trace
		fprintf('In nbody: Argument dump follows\n');
		tmax, level, r0, v0, m0
	end

	% Number of bodies
	nb = length(m0);

	% Number of time steps; create t, r, v,  m, E, K and T matrices
	nt = 2^level + 1;
	t = linspace(0.0, tmax, nt);
	r = zeros(nb, 3, nt);
	v = zeros(nb, 3, nt);
	m = zeros(nb, 1);
	E = zeros(1, nt);
	T = zeros(1, nt);
	V = zeros(1, nt);

	% The time step
	deltat = t(2) - t(1);

	% Initialize the masses of the particles
	m = m0;

	% Initialize the first value of the particles velocities
	v(:,:,1) = v0;

	% Initialize the first two values of the particles positions
	r(:,:,1) = r0;
	r(:,:,2) = r0 + deltat * v0 + 0.5 * deltat^2 * nbodyaccn(r0, m);

	if trace 
		fprintf('deltat=%g \n', deltat);
		fprintf('r(:,:,1)=\n');
		fprintf([repmat('%f\t', 1, size(r(:,:,1), 2)) '\n'], r(:,:,1)');
		fprintf('r(:,:,2)=\n');
		fprintf([repmat('%f\t', 1, size(r(:,:,2), 2)) '\n'], r(:,:,2)');
	end

	%------------------------
	% Perform the simulation
	%------------------------

	for n = 2:nt-1
		% Compute position of the particles at later times
		r(:,:,n+1) = 2 * r(:,:,n) - r(:,:,n-1) + deltat^2 * nbodyaccn(r(:,:,n), m);

		% Compute the velocities of the particles
		v(:,:,n) = ( r(:,:,n+1) - r(:,:,n-1) ) / (2 * deltat);
	end

	% Determine the value of the velocity at the final time step
	v(:,:,nt) = 2 * v(:,:,nt-1) - v(:,:,nt-2);

	% Compute the total kinetic energy at all times
	for s=1:length(t)
		for i=1:nb
			T(s) = T(s) + 0.5 * m(i) * ( (v(i,1,s))^2 + (v(i,2,s))^2 + (v(i,3,s))^2 );
		end
	end

	% Compute the total potential energy at all times
	for l=1:length(t)
		for i=1:nb
                if i+1 <= nb
                        for j=i+1:nb
                                V(l) = V(l) - ( (m(i) * m(j)) / ( (r(j,1,l) - r(i,1,l))^2 + (r(j,2,l) - r(i,2,l))^2 + (r(j,3,l) - r(i,3,l))^2 )^(1/2) );
                        end
                end
		end
        end

	% Compute the total energy at all times
	for l=1:length(t)
		E(l) = T(l) + V(l);
	end

end
