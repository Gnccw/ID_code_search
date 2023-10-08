from cao_pylib import sbox_model
import math
import copy

"""
cnf="(2'+4')(0+2')(1+3')(0+4')(0'+2+4)(0'+1)(1'+3+4)"
cnf_k="(1+2')(1'+2)(0'+2)"
cnf1="(0')(1+2')(1'+2)"
cnf2="(0'+2')(1+2')(0'+1)(0+1'+2)"
sbox_model.gen_ine_from_cnf(cnf2,'ine//num1.txt')
"""
n=38
x=[207 for i in range(n)]
rounds=31
num=n
block_size=256
res=0
for i in x:
    res+=2**(i+1)
a=math.log2(res)-math.log2(rounds)
print(a)
print(math.log2((2**(block_size-num*math.log2(math.e)))))
b=2**a+2**(block_size-num*math.log2(math.e))
print(math.log2(b))

counter=0
s='000000000000000000000000000000000000000000000000000000000000000'
for i in s:
    if i=='0':
        counter+=1
print(counter)
'''
def all_poss_bin_list(a):
    res=[[]]
    res1=[]
    res2=[]
    num=1
    for i in a:
        if i==2:
            num*=2
            tmp=copy.deepcopy(res)
            res=[]
            for j in range(len(tmp)):
                res.append(tmp[j]+[0])
                res.append(tmp[j]+[1])
                res.append(tmp[j]+[2])
        else:
            for j in range(len(res)):
                res[j]=res[j]+[i]
    for i in res:
        if 2 not in i:
            res1.append(i)
        else:
            res2.append(i)
    return res1,res2


s='00000000000000000000000010000000\n00000000000000000100000010000001\n00000000000000000100000010000000\n00000000000000000000000010000001\n00000000000000000000001010000000\n00000000000000000000001000000000\n00000000000000000000001000000001\n00000000000000000000001000000101\n00000000000000000000001000000100'
s=s.split('\n')

res=[]
for i in range(len(s)):
    tmp=[]
    tmp1=[]
    for j in range(len(s[i])):
        tmp1.append(eval(s[i][j]))
        if s[i][j]=='1':
            tmp.append(2)
        else:
            tmp.append(0)
    res.append(tmp)
    s[i]=tmp1

restmp=copy.deepcopy(res)
for i in res:
    res1,res2=all_poss_bin_list(i)
    restmp+=res2


ss=''
for i in restmp:
    tmpres1=[]
    res1,res2=all_poss_bin_list(i)
    for j in res1:
        if 1 in j:
            tmpres1.append(j)
    
    flag=1
    for j in tmpres1:
        if j not in s:
            flag=0
    if flag==1:
        tmp=''
        for k in i:
            tmp+=str(k)
        ss+=("'"+tmp+"'")
        ss+=','
print(ss)

for i in range(3):
    for j in range(4):
        if j<2:
            print(i,j)
        else:
            break
'''