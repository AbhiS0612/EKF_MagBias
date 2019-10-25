function rotmatrix = rotxyz(rph)
% ROTXYZ  Compute a rotation matrix about the XYZ-axes.
%   R = ROTXYZ(RPH) returns [3x3] rotation matrix R.  Note: RPH
%   is a 3-vector of Euler angles [roll,pitch,heading] measured in
%   radians.  RPH measures orientation of coordinate frame 2
%   relative to coordinate frame 1.  Multiplication by rotation
%   matrix R rotates a vector in coordinate frame 2 into coordinate
%   frame 1.
%
%-----------------------------------------------------------------
%    History:
%    Date            Who         What
%    -----------     -------     -----------------------------
%    17 March 2001   LLW         Created and Written
%    12-02-2002      rme         Added matlab help text and changed
%                                rotation matrix to be consistent with Fossen.
%    07-20-2003      rme         Fixed a typo in rotmatrix code comment
%    08-14-2003      rme         Explicitly wrote out rotation matrix to
%                                speed up the calculation.
  
% rotmatrix = yaw_matrix' * pitch_matrix' * roll_matrix'
% rotmatrix = rotz(rph(3))' * roty(rph(2))' * rotx(rph(1))';

% note: the 3 principle rotations above are equivalent to:
cr = cos(rph(1)); sr = sin(rph(1));
cp = cos(rph(2)); sp = sin(rph(2));
ch = cos(rph(3)); sh = sin(rph(3));

rotmatrix = [ ch*cp,   (-sh*cr + ch*sp*sr),   ( sh*sr + ch*sp*cr); ...
	      sh*cp,   ( ch*cr + sh*sp*sr),   (-ch*sr + sh*sp*cr); ...
	      -sp,          cp*sr,                  cp*cr        ];

