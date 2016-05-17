#!/bin/csh -f

set thedir=./v_i
set files = `ls $thedir` 
set outputdir = nothing 
echo $files
foreach file ($files:q) 
        set prog = `echo $thedir\/$file` 
        echo "running for $prog" 
         
        python simulateNDR.py $prog 65536
        python simulateNDR.py $prog 1048576
        python simulateNDR.py $prog 10485760
        python simulateNDR.py $prog 104857600
end 

