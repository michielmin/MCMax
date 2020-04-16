#!/bin/bash
# remove extra linebrakes in the visbility vs. wavelength files (VISIBILITY)
# expects 5 files, H, K, L, M and N band
# will work for the following setup:
# HKLM: 6 wavelengths grid, starting and ending with wavelength as used here
# N 15: 15 wavelengths grid, starting and ending with wavelength as used here
# for any other setup, the wavelength values used in the calculation need to replace the wavelength values provided here

H=$1
K=$2
L=$3
M=$4
N=$5

if [ "$#" -ne 5 ]; then
    echo "please provide the five filenames H K L M N, ordered"
fi

if [[ -n "$H" ]]; then

    sed -i ':a;N;$!ba;s/\n/ /g' $H
    sed -i 's/#/\n#/g' $H
    sed -i 's/1.5000000/\n  1.5000000/g' $H
    sed -i 's/1.5557059/\n  1.5557059/g' $H
    sed -i 's/1.6134806/\n  1.6134806/g' $H
    sed -i 's/1.6734009/\n  1.6734009/g' $H
    sed -i 's/1.7355465/\n  1.7355465/g' $H
    sed -i 's/1.8000000/\n  1.8000000/g' $H
    sed -i 1d $H

else   
    echo "sorry, something went wrong reading the H-band file"  
fi

if [[ -n "$K" ]]; then

    sed -i ':a;N;$!ba;s/\n/ /g' $K
    sed -i 's/#/\n#/g' $K
    sed -i 's/2.0000000/\n  2.0000000/g' $K
    sed -i 's/2.0742745/\n  2.0742745/g' $K
    sed -i 's/2.1513075/\n  2.1513075/g' $K
    sed -i 's/2.2312012/\n  2.2312012/g' $K
    sed -i 's/2.3140620/\n  2.3140620/g' $K
    sed -i 's/2.4000000/\n  2.4000000/g' $K
    sed -i 1d $K

else   
    echo "sorry, something went wrong reading the K-band file"  
fi
    
if [[ -n "$L" ]]; then

    sed -i ':a;N;$!ba;s/\n/ /g' $L
    sed -i 's/#/\n#/g' $L
    sed -i 's/3.0000000/\n  3.0000000/g' $L
    sed -i 's/3.1776715/\n  3.1776715/g' $L
    sed -i 's/3.3658654/\n  3.3658654/g' $L
    sed -i 's/3.5652049/\n  3.5652049/g' $L
    sed -i 's/3.7763500/\n  3.7763500/g' $L
    sed -i 's/4.0000000/\n  4.0000000/g' $L
    sed -i 1d $L

else   
    echo "sorry, something went wrong reading the L-band file"     
fi    

if [[ -n "$M" ]]; then

    sed -i ':a;N;$!ba;s/\n/ /g' $M
    sed -i 's/#/\n#/g' $M
    sed -i 's/4.0000000/\n  4.0000000/g' $M
    sed -i 's/4.1825582/\n  4.1825582/g' $M
    sed -i 's/4.3734482/\n  4.3734482/g' $M
    sed -i 's/4.5730505/\n  4.5730505/g' $M
    sed -i 's/4.7817624/\n  4.7817624/g' $M
    sed -i 's/5.0000000/\n  5.0000000/g' $M
    sed -i 1d $M

else   
    echo "sorry, something went wrong reading the M-band file"     
fi 

if [[ -n "$N" ]]; then

    sed -i ':a;N;$!ba;s/\n/ /g' $N
    sed -i 's/#/\n#/g' $N
    sed -i 's/8.0000000/\n  8.0000000/g' $N
    sed -i 's/8.2822996/\n  8.2822996/g' $N
    sed -i 's/8.5745610/\n  8.5745610/g' $N
    sed -i 's/8.8771355/\n  8.8771355/g' $N
    sed -i 's/9.1903871/\n  9.1903871/g' $N
    sed -i 's/9.5146925/\n  9.5146925/g' $N
    sed -i 's/9.8504419/\n  9.8504419/g' $N
    sed -i 's/10.198039/\n  10.198039/g' $N
    sed -i 's/10.557901/\n  10.557901/g' $N
    sed -i 's/10.930463/\n  10.930463/g' $N
    sed -i 's/11.316171/\n  11.316171/g' $N
    sed -i 's/11.715490/\n  11.715490/g' $N
    sed -i 's/12.128900/\n  12.128900/g' $N
    sed -i 's/12.556898/\n  12.556898/g' $N
    sed -i 's/13.000000/\n  13.000000/g' $N
    sed -i 1d $N

else   
    echo "sorry, something went wrong reading the N-band file"    
fi

echo "all files can now be read with python"
