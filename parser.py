import re
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))
my_file = os.path.join(dir_path, 'rrms.txt')
my_fileout = os.path.join(dir_path, 'output.txt')

out=""

infile = open(my_file,'r')

cols=40
rarray= []
stringOut="table(1)%min_a=1 \n table(1)%z=["

i=1
j=1
k=1

for x in infile:
    x=re.sub("RRMSVEC\("," ",x)
    x=re.sub("D 00"," ",x)
    x=re.sub("\n"," ",x)
    arrx=re.split(",|\)|=",x)

    z=arrx[0].replace(" ","")
    zval=int(z)
    a=arrx[1].replace(" ","")
    aval=int(a)
    r=arrx[3].replace(" ","")

    while k !=zval and j<=cols:
        stringOut+='0.0000'+','
        i+=1
        j+=1     
    
    if j<=cols:
        if i==aval:
            stringOut+=r+','
            i+=1
            j+=1
        else:
            while i != aval and j<=cols: 
                stringOut+='0.0000'+','
                i+=1
                j+=1
            if j>cols:
                i=aval
                j=1
                stringOut = stringOut[:-1]
                stringOut+="] \n table("+z+")%min_a="+a+" \n  table("+z+")%z=["
                k+=1
                while i != aval and j<=cols: 
                    stringOut+='0.0000'+','
                    i+=1
                    j+=1
                if i==aval:
                    stringOut+=r+','
                    i+=1
                    j+=1
            else:
                stringOut+=r+','
                i+=1
                j+=1
    else:
        i=aval
        j=1
        stringOut = stringOut[:-1]
        stringOut+="] \n table("+z+")%min_a="+a+" \n table("+z+")%z=["
        k+=1
        while i != aval and j<=cols: 
            stringOut+='0.0000'+','
            i+=1
            j+=1
        if i==aval:
            stringOut+=r+','
            i+=1
            j+=1
while j<=cols: 
    stringOut+='0.0000'+','
    j+=1
stringOut = stringOut[:-1]
stringOut+="]"

outfile = open(my_fileout, "w")
outfile.write(stringOut)
outfile.close