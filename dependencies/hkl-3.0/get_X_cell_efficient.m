



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function Xloci1 = get_X_cell_efficient(hull,Xs,ind_Xs,d_Xs,p,n);

temps1 = [];
for j=1:p
    temps1(j) = d_Xs{j}(hull(j));
end

Xloci1 = ones(n,prod(temps1));

for j=1:prod(temps1)
    temp = ind2subv(temps1,j);
    for k=1:p
        Xloci1(:,j) = Xloci1(:,j) .* Xs(:,ind_Xs{k}{hull(k)}(temp(k)));
    end
end