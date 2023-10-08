import math
from cao_pylib import operation
import copy

def print_x(m,var_name,num,length):
    """
    var:变量名，字符串的形式
    num:变量的数量
    length:一般是block_size或half_block
    """
    '''以length比特为一行输出某个变量的值,该值存入列表中'''

    tmp=[0 for i in range(num)]
    for v in m.getVars():
        if v.varName[:len(var_name)]==var_name and v.varName[len(var_name)]=='[':
            tmp[eval(v.varName[len(var_name)+1:len(v.varName)-1])]=int(v.x)
    res=[]
    for r in range(int((num/length)/2)):
        s=''
        for i in range(length):
            if i%4==0:
                s+=' '
            s+=str(tmp[(r*length+i)*2]+tmp[(r*length+i)*2+1])
        res.append(s)
    return res

def print_k(m,var_name,num,length):
    tmp=[i for i in range(num)]
    for v in m.getVars():
        if v.varName[:len(var_name)]==var_name and v.varName[len(var_name)]=='[':
            tmp[eval(v.varName[len(var_name)+1:len(v.varName)-1])]=int(v.x)
    
    res=[]
    for r in range(int(num/length)):
        s=''
        for i in range(length):
            if i%4==0:
                s+=' '
            s+=str(tmp[r*length+i])
        res.append(s)
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

def undisturb_sbox_linear(ldt):
    """
    根据线性分布表，生成线性undisturbed bits的分布表
    """

    input_n=int(math.log2(len(ldt)))
    output_n=int(math.log2(len(ldt[0])))
    ddt=ldt
    res=[]
    disturb=set() #该输入下不包含任何不动差分
    
    for in_diff in range(3**input_n):
        #print(in_diff)
        #print(int2tri(n,in_diff))
        inset_bin=dtri2dbin(int2tri(input_n,in_diff))
        inset_int=set()
        for i in inset_bin:
            inset_int=inset_int.union([operation.binlist2num(i)])
        
        tmp=inset_int.union(disturb)
        if len(tmp)!=len(disturb)+len(inset_int):
            disturb=tmp
            break
        else:
            outbin=[0 for i in range(output_n)]
            outnum=0
            for indiff_int in inset_int:
                for outdiff in range(2**output_n):
                    if ddt[indiff_int][outdiff]!=0:
                        outdiff_bin=operation.num2binlist(output_n,outdiff)
                        for i in range(output_n):
                            outbin[i]+=outdiff_bin[i]
                        outnum+=1
        for i in range(output_n):
            if outbin[i]!=0 and outbin[i]!=outnum:
                outbin[i]='?'
            if outbin[i]==outnum:
                outbin[i]=1
            if outbin[i]=='?':
                outbin[i]=2
        res.append(int2tri(input_n,in_diff)+outbin)
    return res


def undisturb_sbox_diff(ddt):
    """
    根据差分分布表，生成差分的undisturbed bits的分布表
    """
    input_n=int(math.log2(len(ddt)))
    output_n=int(math.log2(len(ddt[0])))
    res=[]
    disturb=set() #该输入下不包含任何不动差分
    
    for in_diff in range(3**input_n):
        #print(in_diff)
        #print(int2tri(n,in_diff))
        inset_bin=dtri2dbin(int2tri(input_n,in_diff))
        inset_int=set()
        for i in inset_bin:
            inset_int=inset_int.union([operation.binlist2num(i)])
        
        tmp=inset_int.union(disturb)
        if len(tmp)!=len(disturb)+len(inset_int):
            disturb=tmp
            break
        else:
            outbin=[0 for i in range(output_n)]
            outnum=0
            for indiff_int in inset_int:
                for outdiff in range(2**output_n):
                    if ddt[indiff_int][outdiff]!=0:
                        outdiff_bin=operation.num2binlist(output_n,outdiff)
                        for i in range(output_n):
                            outbin[i]+=outdiff_bin[i]
                        outnum+=1
        for i in range(output_n):
            if outbin[i]!=0 and outbin[i]!=outnum:
                outbin[i]='?'
            if outbin[i]==outnum:
                outbin[i]=1
            if outbin[i]=='?':
                outbin[i]=2
        res.append(int2tri(input_n,in_diff)+outbin)
    return res

def tri2bin_list(x):
    """
    将含有未知变量的列表转换为0/1的点
    """
    res=[]
    for i in range(len(x)):
        if x[i]==0:
            res+=[0,0]
        if x[i]==1:
            res+=[0,1]
        if x[i]==2:
            res+=[1,1]
    return res