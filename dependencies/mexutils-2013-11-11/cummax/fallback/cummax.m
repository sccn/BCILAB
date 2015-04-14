function v = cummax(v)
m = v(1);
for k = 2:numel(v)
    if v(k) <= m
        v(k) = m;
    else
        m = v(k);
    end
end
