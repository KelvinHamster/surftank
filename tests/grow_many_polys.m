rows = 5;
cols = 5;

tiledlayout(rows,cols)

for NPOLY=1:(rows*cols)
    p = test_grow_random_poly();
    nexttile
    plot([p(:,1);p(1,1)],[p(:,2);p(1,2)])
    hold on
    scatter(p(:,1),p(:,2))

end