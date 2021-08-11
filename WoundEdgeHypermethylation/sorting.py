#!/usr/bin/env python

#Copyright (C) <2021>  <The Ohio State University>       

#This program is free software: you can redistribute it and/or modify                              
#it under the terms of the GNU General Public License as published by 
#the Free Software Foundation, either version 3 of the License, or    
#(at your option) any later version.                                                                                       
#This program is distributed in the hope that it will be useful, 
#but WITHOUT ANY WARRANTY; without even the implied warranty of           
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       
#GNU General Public License for more details.                                                                             
#You should have received a copy of the GNU General Public License 
#along with this program.  If not, see <https://www.gnu.org/licenses/>.

#This code was written by Logan Walker

import sys

file_in = open(sys.argv[1], 'r')
file_out = open(sys.argv[2], 'w')

print >>file_out, file_in.readline().strip() #Header

out = []

for line in file_in:
    line = line.strip()
    cols = line.split("\t")

    if float(cols[1]) > 10.0:
        out.append( cols )

out = sorted(out, key=lambda x:x[1])

for line in out:
    print >>file_out, line
