function sumVal = integralSum(I, x1, x2, y1, y2)
  % Computes the sum of pixel values in the region (x1:x2, y1:y2) using an integral image I.
  % I: Integral image (cumulative sum)
  % (x1, y1): Top-left coordinate of the rectangle
  % (x2, y2): Bottom-right coordinate of the rectangle
  % Returns: sumVal, the sum over the specified rectangle.

  % Ensure coordinates are within the image bounds
  [H, W] = size(I);
  
  % Adjust out-of-bounds indices to avoid indexing errors
  x1 = max(x1, 1);
  y1 = max(y1, 1);
  x2 = min(x2, H);
  y2 = min(y2, W);
  
  % Compute the sum using the summed-area table formula
  sumVal = I(x2, y2);
  
  if x1 > 1
      sumVal = sumVal - I(x1 - 1, y2);
  end
  
  if y1 > 1
      sumVal = sumVal - I(x2, y1 - 1);
  end
  
  if x1 > 1 && y1 > 1
      sumVal = sumVal + I(x1 - 1, y1 - 1);
  end
end
