# ANIMATE BALLS for MD codes (disks)

# usage example: "gnuplot -e 'nt=2000; LX=10; LY=10' --persist animate2_ball.gnuplot"

set t x11
#set t png
set xrange[-LX*0.5:LX*0.5]
set yrange[-LY*0.5:LY*0.5]
set size ratio LY/LX

cd '/home/fvega/Datos/pyMD/'

do for[count=0:nt]{
    pause 0.1
    frame = 'xy'.count
    if (count<1000) {frame='xy0'.count}
    if (count<100) {frame='xy00'.count}
    if (count<10) {frame='xy000'.count}
    archivo = frame.'.dat'
#    image = frame.'.png'
#    set o image
    plot archivo u 1:2:(1) w circles title archivo
}
