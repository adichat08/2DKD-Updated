%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 2DKD -- Two-Dimensional Krawtchouk Descriptors
% 
% Copyright (C) 2019, Julian S DeVille, Daisuke Kihara, Atilla Sit, Purdue 
% University, and Eastern Kentucky University.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% All of codes were implemented by Julian S DeVille and Atilla Sit, and checked 
% and run by Daisuke Kihara
% 
% Octave Release: 5.1.0 (or MATLAB Release: 9.6.0 (R2019a))
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% If you use these programs, please cite the following references:
% 
% [1] AUTHORS: Julian S DeVille, Daisuke Kihara, Atilla Sit
% TITLE: 2DKD: A Toolkit for Content Based Local Image Search
% JOURNAL: Source Code for Biology and Medicine, BMC, submitted, 2019.
% 
% [2] AUTHORS: Atilla Sit, Daisuke Kihara
% TITLE: Comparison of Image Patches Using Local Moment Invariants
% JOURNAL: IEEE Transactions on Image Processing, 23(5), 2369-2379, 2014.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% To download the current version of this code, please visit the website:
% <https://github.com/kiharalab/2DKD>
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% computeDesc.m -- This file is part of 2DKD.
% 
% This script computes the two-dimensional Krawtchouk descriptors (2DKD) of 
% order up to 3. The inputs are an image density function f(x,y) and the 
% point-of-interest location (xp,yp). The output vector V contains the six 
% invariant local descriptors.
% 
% See the references [1] and [2] for details.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function V = compDesc(I_f, I_xf, I_yf, I_fW, istart, jstart, S, const)

  % Determine ending indices
  iend = istart + S - 1;
  jend = jstart + S - 1;
  
  % 🔹 Compute moments using the integral images (O(1))
  M00 = integralSum(I_f, istart, iend, jstart, jend);
  M10 = integralSum(I_xf, istart, iend, jstart, jend);
  M01 = integralSum(I_yf, istart, iend, jstart, jend);
  
  % Compute center of mass for the patch
  xtilde = M10 / M00;
  ytilde = M01 / M00;
  
  % 🔹 Compute second-order central moments using integral images
  [X, Y] = meshgrid(0:S-1, 0:S-1);
  x_coords = X - xtilde;
  y_coords = Y - ytilde;
  
  % Extract ftilde from the weighted integral image
  ftilde = integralSum(I_fW, istart, iend, jstart, jend);
  
  % Compute second-order central moments
  mu20 = sum(sum((x_coords.^2) .* ftilde));
  mu02 = sum(sum((y_coords.^2) .* ftilde));
  mu11 = sum(sum((x_coords .* y_coords) .* ftilde));
  
  % 🔹 Compute the rotation-invariant moment-based angle
  mu20MinuSmu02 = mu20 - mu02;
  if mu11 == 0
      if mu20MinuSmu02 == 0
          theta = 0;
      else
          theta = (mu20MinuSmu02 > 0) * 0 + (mu20MinuSmu02 < 0) * (-pi/2);
      end
  else
      ksi = 2 * mu11 / mu20MinuSmu02;
      theta = 0.5 * atan(ksi);
      if mu11 < 0 && mu20MinuSmu02 < 0
          theta = theta - pi/2;
      elseif mu11 > 0 && mu20MinuSmu02 < 0
          theta = theta + pi/2;
      end
  end
  
  % 🔹 Compute Qtilde using precomputed constants
  lambdatilde = zeros(4,4);
  for i = 0:3
      for j = 0:3
          if i + j <= 3
              lambdatilde(i+1,j+1) = sum(sum((X.^i) .* (Y.^j) .* ftilde));
          end
      end
  end
  
  % Normalize lambdatilde
  lambdatilde = lambdatilde / M00;
  
  % Compute Qtilde
  Qtilde = zeros(4,4);
  for n = 0:3
      for m = 0:3
          if n+m <= 3
              Qtilde(n+1,m+1) = const.a(1:n+1,n+1)' * lambdatilde(1:n+1,1:m+1) * const.a(1:m+1,m+1);
          end
      end
  end
  
  % Normalize Qtilde
  Qtilde = Qtilde ./ sqrt(const.rho' * const.rho);
  
  % Construct final descriptor vector (using non-redundant descriptors)
  V = [ Qtilde(3,1) Qtilde(1,3) Qtilde(2,3) Qtilde(3,2) Qtilde(4,1) Qtilde(1,4) ];
end
