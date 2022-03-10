from cao_pylib import operation

def schedule(x):
    res=[x]
    for i in range(35):
        tmp=[]
        tmp.append(x[1])
        tmp.append(x[0])
        tmp.append(x[4])
        tmp.append(x[2])
        tmp.append(x[3])
        tmp.append(x[5])
        for i in range(6,16):
            tmp.append(x[i])
        tmp=operation.roate_left(tmp, 6)
        x=tmp
        res.append(x)
    return res
res=schedule('1000000000000000')

def print_x(x):
    counter=1
    for i in x:
        tmp=''
        for j in range(len(i)):
            if j%4==0:
                tmp+=' '
            tmp+=i[j]
        print(str(counter),tmp)
        counter+=1
print_x(res)