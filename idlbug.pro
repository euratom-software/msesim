set_plot,'PS'
!p.font=0
device,/Color, File='bug.eps', /encapsulated

plot, [0,1.0],[0,1.249e4],       $
        xtitle='xtitle', $
        ytitle='ytitle',         $
        title ='graphs title'

device,/close

end
