import os
import sys

f = open('log.txt', 'w')
f.close()
for j in range(5,100,5):
    for i in range(10):
        os.system('iDBSCAN '+ str(j) +'%-miss-' + str(i) + '.txt ' + sys.argv[1])
    f = open('log.txt','a')
    f.write('-1 -1\n')
    f.close()
    
#os.system('log.txt');
os.system('average < log.txt > output2.txt')