function nbodyout(fname, t, r, m, rgb)%
% nbodyout Output N-body data in xfpp3d format
%
% nbodyout writes the values computed by a typical N-body code to the file 
% fname in the format required by xfpp3d.
%
% Input arguments
%
%   fname (character string): Name of output file.
%
%   t (real 1D array, length nt): Contains output times.  At the current
%     time, nt must be >= 2 (see note below).
%
%   r (real 3D array, dimensions np x 3 x nt): (x,y,z) coordinates of
%     all particles at all time steps. np is the number of particles. 
%     nt is the number of time steps stored. 
%     
%     (r(i,1,n), r(i,2,n), r(i,3,n)) are the x, y and z coordinates, 
%     respectively, of the i-th particle at the n-th time step.
% 
%     NOTE: In MATLAB/octave it is apparently impossible to define 
%     a 3D array with a third dimension of length 1: for example,
%     although zeros will happily appear to create one if you type 
%     'zeros(2,2,1)', the result will actually be a 2D array with 
%     size 2 x 2.  Although nt = 1 could be implemented as a special 
%     case here, I have not yet done so.  Thus the restriction to 
%     nt >= 2.
%
%   m (OPTIONAL, real 1D array, dimension np): Defines masses of 
%     particles.  Particles will be visualized by xfpp3d using markers
%     with (linear) sizes that are PROPORTIONAL to m(i).  Thus, you 
%     may want to rescale m relative to the true masses should 
%     max(m) / min(m) be large, so that the particle sizes do not 
%     span too large a range.
%
%     If this argument is not supplied, then m(i) = 1 will be used for 
%     all particles.
%
%   rgb (OPTIONAL, real 2D array, dimension np x 3): Defines RGB color
%     that will be used for each particle.  (rgb(i,1), rgb(i,2), rgb(i,3))
%     defines the RGB triple, (red, green, blue), for particle i.  
%     Each RGB value must be in the range 0.0 to 1.0.  
%
%     NOTE: xfpp3d must be invoked with the -c option if RGB information 
%     is included in the input file.
%
%     If this argument is not supplied then RGB information will NOT be 
%     included in the file, and xfpp3d should not be called with the 
%     -c argument.
%
% Output arguments
%
%   None
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   % BEGIN argument processing / checking ...
   err = 0;
   np = 0;
   nt = 0;

   % Need at least three arguments ...
   if nargin < 3 
      fprintf('nbodyout: Must be called with at least 3 arguments.\n');
      return;
   end
 
   % Check that fname is a string ...
   if ~ ischar(fname)
      fprintf('nbodyout: First argument [fname] must be a string.\n');
      err = 1;
   end

   % Check that t is a vector; if it is, set nt to its length,
   % and ensure that nt >= 2 ...
   if ~ isvector(t) 
      fprintf('nbodyout: Second argument [t] must be a vector.\n');
      err = 1;
   else
      nt = length(t);
      if nt < 2 
         fprintf('nbodyout: Number of time steps [length(t)=%d] must be >= 2.\n', nt);
         err = 1;
      end
   end

   % Check dimensions of r ...
   szr = size(r);
   szszr = size(size(r));
   rankr = szszr(2);
   if rankr ~= 3
      fprintf('nbodyout: Third argument [r] must be a 3D array.\n');
      err = 1;
   else
      if szr(2) ~= 3
         fprintf('nbodyout: Second dimension of third argument [r] must be 3.\n');
         err = 1;
      else
         % Determine number of particles from dimensions of r and 
         % ensure that length of 3rd dimension is nt ...
         np = szr(1);
         if szr(3) ~= nt 
            fprintf('nbodyout: Third dimension of third argument [r] ');
            fprintf('should be %d but is %d.\n', nt, szr(3));
            err = 1;
         end
      end
   end

   % If array of masses has been supplied, ensure that it is a vector
   % of length np; if it hasn't define it to be a vector of that length
   % with all elements = 1 ...
   if np
      if nargin > 3
         if ~ isvector(m)
            fprintf('nbodyout: Fourth argument [m] must be a 1d array.\n');
            err = 1;
         else
            if length(m) ~= np
               fprintf('nbodyout: Length of fourth argument [m] is %d,', ...
                  length(m));
               fprintf(' but should be %d.\n', np);
               err = 1;
            end
         end
      else
         m = ones(1,np);
      end
   end

   % If rgb array has been supplied, ensure that it is np x 3 and 
   % that all elements are in the range 0.0 - 1.0.  Set flag rgbout
   % to indicate whether or not RGB info is to be written to the 
   % file ...
   if np & nargin > 4
      szrgb = size(rgb);
      if ~ (szrgb(1) == np & szrgb(2) == 3)
            fprintf('nbodyout: Fifth argument [rgb] must be a ');
            fprintf('%d x %d array, but is %d x %d.\n', np, 3, ...
               szrgb(1), szrgb(2));
            err = 1;
      else
         if min(min(rgb)) < 0.0 | max(max(rgb)) > 1.0
            fprintf('nbodyout: All elements of fifth argument [rgb] must ');
            fprintf('be in range 0.0 to 1.0.\n');
            err = 1;
         end
      end
      rgbout = 1;
   else
      rgbout = 0;
   end

   if err 
      fprintf('nbodyout: One or more errors detected in input arguments.\n');
      fprintf('nbodyout: Aborting without writing to file.\n');
      return
   end

   % ... END argument processing / checking

   % Open file for write ...
   fid = fopen(fname, 'w');
   if fid < 0 
      fprintf('nbodyout: Error opening ''%s'' for write.\n', fname);
      return;
   end

   % Write number of particles ...
   fprintf(fid, '%d\n', np);

   % Write particle masses, and, if supplied, RGB values ...
   if rgbout
      for i = 1 : np
         fprintf(fid, '%g %g %g %g\n', m(i), rgb(i,1), rgb(i,2), rgb(i,3));
      end
   else
      for i = 1 : np
         fprintf(fid, '%g\n', m(i));
      end
   end

   % For each time step ...
   for n = 1 : nt
      % Write time ...
      fprintf(fid, '%g\n', t(n));
      % Write (x,y,z) for all particles at time step n ...
      for j = 1 : np
         fprintf(fid, '%g %g %g\n', r(j,:,n));
      end
   end

   % Close file and return ...
   fclose(fid);
end
