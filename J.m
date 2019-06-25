%Skew Symmetric Operator
function skew=J(e) 
   skew=[[0 -e(3,1) e(2,1)];[e(3,1) 0 -e(1,1)];[-e(2,1) e(1,1) 0]];
end