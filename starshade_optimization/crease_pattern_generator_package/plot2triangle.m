function plot2triangle(points, fig)
   n = size(points, 2);
   
   figure(fig);
   hold on
   
   for i = 1:(n - 2)
       %plot(points(1, :), points(2, :), '.')
       
       line([points(1, i) points(1, i + 1)],...
            [points(2, i) points(2, i + 1)])
        
       line([points(1, i + 1) points(1, i + 2)],...
            [points(2, i + 1) points(2, i + 2)])
        
       line([points(1, i + 2) points(1, i)],...
            [points(2, i + 2) points(2, i)])
   end
   
   hold off   
   axis equal
end