from cao_pylib import operation

'''
s3=[0,0,0,1,0,1,1,1]
s2=[0,0,0,1]
t=[[0 for i in range(2)] for i in range(8)]

for x1 in range(8):
    for x2 in range(8):
        t[x1^x2][s3[x1]^s3[x2]]+=1
for i in t:
    print(i)

t2=[[0 for i in range(2)] for i in range(4)]
for x1 in range(4):
    for x2 in range(4):
        t2[x1^x2][s3[x1]^s3[x2]]+=1
for i in t2:
    print(i)

'''
for diff in range(16):
    for x1 in range(16):
        x2=diff^x1
        for key in range(16):
            diff1=((x1+key)%16)^((x2+key)%16)
            if diff!=diff1:
                print(diff,x1,x2,key,diff1)