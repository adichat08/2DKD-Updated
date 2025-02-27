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
% dbSearch.m -- This file is part of 2DKD.
% 
% This script is responsible for searching the output of 'dbIndex' for 
% descriptors similar to the ones corresponding the query. A query image is 
% supplied as input, then 'compDesc' is run on the query, producing descriptors 
% for it, and then the matrix from 'dbIndex' is sorted by Euclidean distance of 
% descriptors to the new ones obtained, giving a ranked list of the most similar
% regions to the query from all subimages in the database.
% 
% See the references [1] and [2] for details.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function top_matches = dbSearch(queryPath, DB, dbPath, imageList, S, frameSize, xp, yp)

  % 🔹 Read the query image and retrieve integral images
  [fQuery, ~, I_f, I_xf, I_yf, I_fW] = readImage(queryPath, prepStep(S));

  % 🔹 Compute the query descriptor using integral images
  Vquery = compDescDP(I_f, I_xf, I_yf, I_fW, xp, yp, S, prepStep(S));

  % Number of entries in the database
  dbSize = size(DB, 1);

  % 🔹 Compute Euclidean distances efficiently
  EucDist = sum((repmat(Vquery, dbSize, 1) - DB(:, 4:9)).^2, 2);

  % Append Euclidean distances to the database for sorting
  DB = [DB EucDist];

  % Sort the database by Euclidean distances
  DB = sortrows(DB, 10);

  % 🔹 Remove redundant matches in the top 2000 results
  k = min(dbSize, 2000);
  DB = DB(1:k, :);
  dbSize = size(DB, 1);

  Ind = ones(1, k);
  for i = 1:k-1
      if Ind(i) == 1
          % Faster way to check for nearby duplicate regions
          T = abs(DB(i+1:k, 1:3) - DB(i, 1:3)) <= [0 10 10];
          T = sum(T, 2);
          [T, ~] = find(T == 3);
          Ind(T + i) = 0;
      end
  end
  DB(Ind == 0, :) = [];
  dbSize = size(DB, 1);

  % 🔹 Select top 15 matches for final output
  k = min(dbSize, 15);
  top_matches = [repmat('  ', [k 1]), num2str((1:k)'), repmat('    ', [k 1]), ...
                 imageList(DB(1:k, 1), :), repmat('  ', [k 1]), num2str(DB(1:k, 2:3))];

  % 🔹 Visualization (unchanged)
  X = 0:frameSize-1;
  Y = X;
  [X, Y] = ndgrid(X, Y);

  [~, fileName, ext] = fileparts(queryPath);
  s = get(0, 'ScreenSize');
  figure('Position', [10 10 s(3)-100 s(4)-100]);

  subplot(3,6,1);
  [fQueryLocal,~,~] = squareCrop(fQuery, xp, yp, frameSize);
  pl = pcolor(X, -Y, fQueryLocal);
  axis equal tight;
  colormap gray;
  set(pl, 'edgecolor', 'none');
  axis off;
  title({'Query', [fileName ext], ['(x,y) = (' num2str(xp) ',' num2str(yp) ')']}, ...
        'Interpreter', 'none', 'Fontsize', 8);

  % 🔹 Display Top 15 Matches
  for i = 1:k
      [~, fileName, ext] = fileparts(imageList(DB(i, 1), :));
      [fMatch, ~, I_f, I_xf, I_yf, I_fW] = readImage([dbPath fileName ext(1:4)], prepStep(S));
      xMatch = DB(i, 2);
      yMatch = DB(i, 3);
      [fMatchLocal,~,~] = squareCrop(fMatch, xMatch, yMatch, frameSize);

      subplotIdx = i + 1 + (i > 5) + (i > 10);
      subplot(3, 6, subplotIdx);
      pl = pcolor(X, -Y, fMatchLocal);
      axis equal tight;
      colormap gray;
      set(pl, 'edgecolor', 'none');
      axis off;
      title({['Top ' num2str(i)], [fileName ext], ['(x,y) = (' num2str(xMatch) ',' num2str(yMatch) ')']}, ...
            'Interpreter', 'none', 'Fontsize', 8);
  end
end
