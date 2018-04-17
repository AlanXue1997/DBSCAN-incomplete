import os
import sys

for i in range(10):
    os.system('DBSCAN-incomplete '+ sys.argv[1] +'%-miss-' + str(i) + '.txt ' + sys.argv[2])
    
os.system('log.txt');