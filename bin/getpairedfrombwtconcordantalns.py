#!/usr/bin/env python
import csv
import sys

f = csv.reader(sys.stdin, dialect="excel-tab")
of = csv.writer(sys.stdout, dialect="excel-tab")
last_read = None
c = 0
str = []
for line in f :
    #take care of the header
    if(line[0][0] == "@") :
        of.writerow(line)
        continue

    if (last_read == None) :
        last_read = line
        c = 0
        str.append(line)
        continue
        
    if (last_read[0] != line[0]) :
        if (c == 1):
            of.writerow(str[0])
            of.writerow(str[1])
        str = []
        c = 0
        str.append(line)
    else :
        c += 1
        str.append(line)

    last_read = line

if (c == 1):
    of.writerow(str[0])
    of.writerow(str[1])

#     if(last_read == None) : 
#         last_read = line
#     else :
#         if(last_read[0] == line[0]) :
#             of.writerow(last_read)
#             of.writerow(line)
#             last_read = None
#         else :
#             last_read = line
