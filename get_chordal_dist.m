function r = get_chordal_dist(x, y, z)

n = length(x);
r = zeros(n);
for j = 1:n
    for i = 1:j
        s = [x(i); y(i); z(i)];
        t = [x(j); y(j); z(j)];
        r(i ,j) = norm(s-t);
    end
end

end