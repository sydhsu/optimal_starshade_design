function plot3triangle(points, fig)
   n = size(points, 2);
   
   figure(fig);
   hold on
   
   for i = 1:(n-2)
       plot3(points(1, :), points(2, :), points(3, :), '.')
       
       line([points(1, i) points(1, i + 1)],...
            [points(2, i) points(2, i + 1)],...
            [points(3, i) points(3, i + 1)])
        
       line([points(1, i + 1) points(1, i + 2)],...
            [points(2, i + 1) points(2, i + 2)],...
            [points(3, i + 1) points(3, i + 2)])
        
       line([points(1, i + 2) points(1, i)],...
            [points(2, i + 2) points(2, i)],...
            [points(3, i + 2) points(3, i)])
   end
   
   hold off
   
   axis equal
end