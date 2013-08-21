function start_display(ring,fig)
  
 set(fig,'XData',(1:length(ring)),'YData',ring); axis([1,length(ring),-0.3,0.3])
 drawnow
 
return
