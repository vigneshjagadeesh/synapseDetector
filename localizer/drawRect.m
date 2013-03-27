function drawRect(x,y,w,h, inpText, colVector)
hold on;
rectangle('Position', [x y w h], 'LineWidth', 5, 'EdgeColor', colVector );
text(x, y, inpText, 'HorizontalAlignment','left', 'BackgroundColor',[.7 .9 .7]);
hold off;