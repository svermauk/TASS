# To compute cummulative average of system temperature.
awk '{print NR,$1,(p+=$1)/NR}'  OFS='\t' sys_temp.dat >out.dat
