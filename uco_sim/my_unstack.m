function A = my_unstack(stacked, rows)

% 2019-05-22 LLW returns matrix with specified number of rows from stacked vec
%                does not check to verify conformability of arguments

cols = length(stacked)/rows;
  
for i = 0:(cols-1)
  A(:,i+1) = stacked((i*rows)+1:(i*rows)+rows,1);
end

  

