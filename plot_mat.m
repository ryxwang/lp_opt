function [] = plot_mat(A, bd, start_ind, col)
bd1 = bd-start_ind+1;
imagesc(A)
colorbar()
hold on
for i=1:length(bd)
line([bd1(i,1), bd1(i,2)], [bd1(i,1), bd1(i,1)], 'Color', col, 'Linewidth',5)
line([bd1(i,2), bd1(i,2)], [bd1(i,1), bd1(i,2)], 'Color', col, 'Linewidth',5)
end
hold off
