 #!/bin/bash

fractalParam="-v -w 1920 -h 1080 -t 5000 -r -1.992561950824231 -i 0"
# Maximum zoom should be 2.5*10^9

# fractalParam="-v -w 1920 -h 1080 -t 5000 -j -A -1.9999 -B 0 -c 4"
# Maximum zoom should be 1.5*10^19

maxZoom=$(echo $(bc <<< "2.5*10^9") | awk '{print int($1+0.5)}')
base="1.1"
exponent=1
zoomLevel=1
while [ $(bc <<< "${zoomLevel} < ${maxZoom}") -gt 0 ]
do
    echo "--------------------------------------------"
    echo Zoom level: $(bc <<< "${base}^${exponent}")
    ./liten_fraktaler $fractalParam -s $(bc <<< "${base}^${exponent}") -o $exponent
    zoomLevel=$(echo $(bc <<< "${base}^${exponent}") | awk '{print int($1+0.5)}')
    exponent=$((exponent+1))
    #echo $var | awk '{print int($1+0.5)}'
done

