function [a] = nbodyaccn(r, m)
	% The number of particles
	n = length(m);

	% Initiate the acceleration matrix
	A = zeros(size(r));

	% Compute the acceleration matrix
	for i = 1:n
		l = 1:n;
		l(i) = [];
		for j = l
			direction = r(j,:) - r(i,:);
			magr = sum(direction.^2)^(1/2);
			coupling = m(j) / magr^3;
			A(i,:) = coupling * direction(:);
		end
	end

	a = A;

end
