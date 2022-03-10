from cao_pylib import sbox_model
from cao_pylib import sbox
from cao_pylib import operation
from cao_pylib import cipher_sbox
import copy
import numpy as np
import math
from gurobipy import *

_3xor=[[0,0,0,0],[0,0,1,1],[0,0,2,2],[0,1,0,1],[0,1,1,0],[0,1,2,2],[1,0,0,1],[1,0,1,0],[1,0,2,2],[1,1,0,0],
       [1,1,1,1],[1,1,2,2],[0,2,0,2],[0,2,1,2],[0,2,2,2],[1,2,0,2],[1,2,1,2],[1,2,2,2],[2,1,0,2],[2,1,1,2],
       [2,1,2,2],[2,0,0,2],[2,0,1,2],[2,0,2,2],[2,2,0,2],[2,2,1,2],[2,2,2,2]]
_2and=[[0,0,0],[0,1,2],[0,2,2],[1,0,2],[1,1,2],[1,2,2],[2,0,2],[2,1,2],[2,2,2]]
_xor=[[0,0,0],[0,1,1],[0,2,2],[1,0,1],[1,1,0],[1,2,2],[2,1,2],[2,0,2],[2,2,2]]
_word_xor=[[0,0,0],[0,1,1],[0,2,2],[1,0,1],[1,1,2],[1,2,2],[2,1,2],[2,0,2],[2,2,2]]
_word_3xor=[[0,0,0,0],[0,0,1,1],[0,0,2,2],[0,1,0,1],[0,1,1,2],[0,1,2,2],[1,0,0,1],[1,0,1,2],[1,0,2,2],[1,1,0,2],
       [1,1,1,2],[1,1,2,2],[0,2,0,2],[0,2,1,2],[0,2,2,2],[1,2,0,2],[1,2,1,2],[1,2,2,2],[2,1,0,2],[2,1,1,2],
       [2,1,2,2],[2,0,0,2],[2,0,1,2],[2,0,2,2],[2,2,0,2],[2,2,1,2],[2,2,2,2]]


#将一组列表写入到文件中
def write_list_to_file(a,filename):
    s=''
    for i in range(len(a)):
        for v in a[i]:
            s+=str(v)+','
        s=s[:len(s)-1]+'\n'
    fp=open(filename,'w')
    fp.write(s[:len(s)-1])
    fp.close()


#从文件中读取一组列表
def read_list_from_file(filename):
    fp=open(filename,'r')
    s=fp.read()
    fp.close()

    res=[]
    s=s.split('\n')
    for i in s:
        i=i.split(',')
        tmp=[]
        for j in i:
            tmp.append(eval(j))
        res.append(tmp)
    return res

#向量a,b的点乘
def mul(a,b):
    assert len(a)==len(b),'mul err'
    res=0
    for i in range(len(a)):
        res+=a[i]*b[i]
    return res

#将list转换位一行连续的字符
def list2str(x):
    s=''
    for i in range(len(x)):
        s+=str(x[i])
    return s


def str2list(s):
    res=[]
    for i in range(len(s)):
        res.append(eval(s[i]))
    return res


#将形如[2,1]转变位[1,1,0,1]
def triList2binList(a):
    res=[]
    for i in range(len(a)):
        if a[i]==0:
            res+=[0,0]
        if a[i]==1:
            res+=[0,1]
        if a[i]==2:
            res+=[1,1]
    return res


def test_ine(inefile,posspoint):
    """
    测试不等式模型是否正确
    """
    imposspoint=[]
    width=int(len(posspoint[0])/2)
    for i in range(3**width):
        tmppoint=int2tri(width,i)
        if triList2binList(tmppoint) not in posspoint:
            imposspoint.append(triList2binList(tmppoint))
    
    m=Model()
    x=m.addVars([i for i in range(width*2)],vtype=GRB.BINARY,name='x')
    sbox_model.gen_model_from_ine(m,x,inefile)

    for point in posspoint:
        for i in range(width*2):
            m.addConstr(x[i]==point[i],name='c'+str(i))
        m.Params.OutputFlag=0
        m.optimize()
        if m.status==GRB.Status.INFEASIBLE:
            print('poss->imposs err:',point)
        for i in range(width*2):
            m.remove(m.getConstrByName('c'+str(i)))
        m.update()
    
    for point in imposspoint:
        for i in range(width*2):
            m.addConstr(x[i]==point[i],name='c'+str(i))
        m.Params.OutputFlag=0
        m.optimize()
        if m.status!=GRB.Status.INFEASIBLE:
            print('impposs->poss err:',point)
        for i in range(width*2):
            m.remove(m.getConstrByName('c'+str(i)))
        m.update()
    print('---right---')


def int2tri(width,n):
    #将10进制数转为3进制数
    assert 3**width-1>=n,"int2tri:number too big"
    res=[0 for i in range(width)]

    for i in range(width-1,-1,-1):
        res[i]=n%3
        n=int(n/3)
        if n<3:
            res[i-1]=n
            break
    return res


def dtri2dbin(a):
    #将包含不定变量的列表转换为所有可能取值的集合
    width=len(a)
    num=0
    index=[]
    res=[]
    
    for i in range(width):
        if a[i]==2:
            num+=1
            index.append(i)
    
    for i in range(2**num):
        tmp=copy.deepcopy(a)
        ibin=operation.num2binlist(num,i)
        for j in range(num):
            tmp[index[j]]=ibin[j]
        res.append(tmp)
    return res


def undisturb_sbox(s):
    n=sbox.get_input_size_sbox(s)
    ddt=sbox.gen_ddt(s)
    res=[]
    disturb=set() #该输入下不包含任何不动差分
    
    for in_diff in range(3**n):
        #print(in_diff)
        #print(int2tri(n,in_diff))
        inset_bin=dtri2dbin(int2tri(n,in_diff))
        inset_int=set()
        for i in inset_bin:
            inset_int=inset_int.union([operation.binlist2num(i)])
        
        tmp=inset_int.union(disturb)
        if len(tmp)!=len(disturb)+len(inset_int):
            disturb=tmp
            break
        else:
            outbin=[0 for i in range(n)]
            outnum=0
            for indiff_int in inset_int:
                for outdiff in range(2**n):
                    if ddt[indiff_int][outdiff]!=0:
                        outdiff_bin=operation.num2binlist(n,outdiff)
                        for i in range(n):
                            outbin[i]+=outdiff_bin[i]
                        outnum+=1
        for i in range(n):
            if outbin[i]!=0 and outbin[i]!=outnum:
                outbin[i]='?'
            if outbin[i]==outnum:
                outbin[i]=1
            if outbin[i]=='?':
                outbin[i]=2
        res.append(int2tri(n,in_diff)+outbin)
    return res


def undisturb_sbox_invert(s):
    invert_s=[0 for i in range(len(s))]
    for i in range(len(s)):
        invert_s[s[i]]=i
    return undisturb_sbox(invert_s)
        

def write_patterns_csv(patterns,filename):
    """
    将所有的patterns写入到csv文件中，变量名用1，2，3，4,....
    """
    point_num=len(patterns[0])
    input_num=int(point_num/2)
    s=''
    for i in range(point_num):
        s+=str(i)+','
        if i==input_num-1:
            s+=','
    s=s[:len(s)-1]+'\n'
    
    outputset=dict()
    for i in range(len(patterns)):
        output=''
        tmp=''
        for j in range(point_num):
            tmp+=str(patterns[i][j])+','
            if j==input_num-1:
                tmp+=','
            if j>input_num-1:
                output+=str(patterns[i][j])
        tmp=tmp[:len(tmp)-1]+'\n'
        s+=tmp
        outputset[output]=0

    fp=open(filename,'w')
    fp.write(s[:len(s)-1])
    fp.close()
    return outputset


def trans_undistrub_bin(undisturb):
    res=[]
    for i in undisturb:
        res.append(trilist2binlist(i))
    return res


def gen_csv_undisturb_bit_under8bit(s,sbox_name):
    input_size=sbox.get_input_size_sbox(s)
    undisturbbit=undisturb_sbox_invert(s)

    for i in range(input_size):
        tmp=[]
        for j in range(len(undisturbbit)):
            tmp1=[]
            for k in range(input_size):
                tmp1.append(undisturbbit[j][k])
            tmp1.append(undisturbbit[j][input_size+i])
            tmp.append(tmp1)
            tmp2=[0 for k in range(len(undisturbbit))]
        print(tmp)
        for k in range(len(undisturbbit)):
            tmp2[k]=trilist2binlist(tmp[k])
        sbox.write_patterns_csv(tmp2,sbox_name+str(i)+'.csv')    


def _int2x(x,a,n):
    """
    进制转换，将10进制转换为x进制，位数为n,高位在前
    """
    assert x<a**n,"number too large!"
    res=[0 for i in range(n)]

    counter=n-1
    while x>=a:
        res[counter]=x%a
        x=int(x/a)
        counter-=1
    res[counter]=x
    return res


def xor2csv(n):
    #n比特不变差分比特的异或
    res=[]
    resbin=[]
    for i in range(3**n):
        inputs=_int2x(i,3,n)
        if 2 in inputs:
            inputs.append(2)
        else:
            output=0
            for i in range(n):
                output^=inputs[i]
            inputs.append(output)
        res.append(inputs)
    for i in res:
        resbin.append(trilist2binlist(i))
    sbox.write_patterns_csv(resbin,"csv//"+'_'+str(n)+'xor.csv')


def _dia_array(a):
    b=copy.deepcopy(a)
    if len(b)==1:
        return b
    else:
        if b[0][0]!=1:
            for row in range(1,len(b)):
                if b[row][0]==1:
                    tmp=a[row]
                    b[row]=b[0]
                    b[0]=tmp
                    break
        for row in range(1,len(b)):
            if b[row][0]==1:
                for col in range(len(a[0])):
                    b[row][col]^=b[0][col]
        c=[[0 for i in range(len(a[0])-1)] for i in range(len(a)-1)]
        for row in range(1,len(b)):
            for col in range(1,len(b[0])):
                c[row-1][col-1]=b[row][col]
        d=_dia_array(c)
        for row in range(1,len(b)):
            for col in range(1,len(b[0])):
                b[row][col]=d[row-1][col-1]
        return b


def inv_array(a):
    """
    生成异或操作的可逆矩阵
    """
    b=copy.deepcopy(a)

    num=len(b)
    for i in range(num):
        tmpb=[0 for i in range(num)]
        tmpb[i]=1
        b[i]+=tmpb
    
    b=_dia_array(b)
    res=[[0 for i in range(len(a))] for i in range(len(a))]
    for i in range(len(a)-1,-1,-1):
        for j in range(i):
            if b[j][i]!=0:
                for col in range(len(b[0])):
                    b[j][col]^=b[i][col]
    return b
    for row in range(len(b)):
        for col in range(len(b),len(b[0])):
            res[row][col-len(b)]=b[row][col]
    return res


def print_table(a):
    for i in a:
        print(i)


def tri_pset2bin_pset(a):
    b=[]
    for i in a:
        b.append(triList2binList(i))
    return b


'''
for i in range(8):
    res=undisturb_sbox_invert(cipher_sbox.lblock_sbox[i])
    posspoint=[]
    for j in range(len(res)):
        posspoint.append(triList2binList(res[j]))
    sbox.write_patterns_csv(posspoint,'lblock/lblock_sbox'+str(i)+'_invert.csv')
'''

def undisturb_sbox_linear(s):
    n=sbox.get_input_size_sbox(s)
    ddt=sbox.gen_ldt(s)
    res=[]
    disturb=set() #该输入下不包含任何不动差分
    
    for in_diff in range(3**n):
        #print(in_diff)
        #print(int2tri(n,in_diff))
        inset_bin=dtri2dbin(int2tri(n,in_diff))
        inset_int=set()
        for i in inset_bin:
            inset_int=inset_int.union([operation.binlist2num(i)])
        
        tmp=inset_int.union(disturb)
        if len(tmp)!=len(disturb)+len(inset_int):
            disturb=tmp
            break
        else:
            outbin=[0 for i in range(n)]
            outnum=0
            for indiff_int in inset_int:
                for outdiff in range(2**n):
                    if ddt[indiff_int][outdiff]!=0:
                        outdiff_bin=operation.num2binlist(n,outdiff)
                        for i in range(n):
                            outbin[i]+=outdiff_bin[i]
                        outnum+=1
        for i in range(n):
            if outbin[i]!=0 and outbin[i]!=outnum:
                outbin[i]='?'
            if outbin[i]==outnum:
                outbin[i]=1
            if outbin[i]=='?':
                outbin[i]=2
        res.append(int2tri(n,in_diff)+outbin)
    return res


def undisturb_sbox_invert_linear(s):
    invert_s=[0 for i in range(len(s))]
    for i in range(len(s)):
        invert_s[s[i]]=i
    return undisturb_sbox_linear(invert_s)

p=[[0,1,3,6,7,4,5,2],
   [0,1,7,4,3,6,5,2],
   [0,3,1,6,7,5,4,2],
   [0,7,3,5,1,4,6,2]]

csv=[]
res=undisturb_sbox_invert(cipher_sbox.skinny_sbox_4bit)
'''
for i in res:
    if 0 in i[4:] or 1 in i[4:]:
        print(i)

for i in res:
    csv.append(triList2binList(i))
#sbox.write_patterns_csv(csv,'print/print_sbox_inv.csv')
'''


ssbox=cipher_sbox.t_twine_pi
print(ssbox)
s=undisturb_sbox_invert(ssbox)
for i in s:
    print(i)