

# Configuración del estilo de la línea y puntos
set style line 1 lc rgb "black" lw 2 pt 7 ps 1.5 # Puntos negros y gruesos
set style line 2 lc rgb "red" lw 1 pt 7 ps 1.5   # Puntos rojos y gruesos

# Configuración del título
set title "Pulso gaussiano"

# Configuración de los ejes
set xlabel "t"
set ylabel "Energia
set grid

# Cargar los datos de la tabla
# Asumiendo que los datos están en el archivo "datos.txt" y están separados por espacios
# La primera columna es X, la segunda columna es Y, la tercera columna es el error de Y
plot "Energia.txt" using 1:2 with line ls 1 lc rgb "black" title " "

