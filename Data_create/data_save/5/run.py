import os
import sys

for i in range(10):
    os.system('DBSCAN-incomplete '+ sys.argv[1] +'%-miss-' + str(i) + '.txt 0.047')
    
os.system('log.txt');