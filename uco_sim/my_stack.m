function stacked = my_stack(A)
% 2019-05-22 LLW returns stack of argument array

[rows, cols] = size(A);

for i = 0:(cols-1)
  stacked((i*rows)+1:(i*rows)+rows,1) = A(:,i+1);
end

  

