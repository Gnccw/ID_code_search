from cao_pylib import sbox
from cao_pylib import operation
import copy


def _mul(x,y):
    """
    多项式的乘法，下标是set，常量是int [{0,1},1]表示：x0x1+1,[{0,1},{1,2},1]:表示x0x1+x1x2+1
    """
    if len(x)==0:
        return y
    if len(y)==0:
        return x
    res=[]
    for i in range(len(x)):
        for j in range(len(y)):
            if type(x[i])==set and type(y[j])==set:
                res.append(x[i].union(y[j]))
            elif type(x[i])==set and y[j]==1:
                res.append(x[i])
            elif type(y[j])==set and x[i]==1:
                res.append(y[j])
            else:
                res.append(1)
    return res


def _merge(x):
    """
    消去多项式中多余的式子
    """
    x=copy.deepcopy(x)
    ss={}
    res=[]
    for i in range(len(x)):
        if type(x[i])==set:
            s=''
            x[i]=list(x[i])
            x[i].sort()
            for j in x[i]:
                s+=str(j)+','
            ss[s]=0
            x[i]=s
        else:
            ss[x[i]]=0
    for i in range(len(x)):
        ss[x[i]]+=1
    for key in ss.keys():
        if type(key)==str:
            if ss[key]%2!=0:
                tmp=set()
                key=key.split(',')
                key.remove('')
                for k in key:
                    tmp.add(int(k))
                res.append(tmp)
        else:
            if ss[key]%2!=0:
                res.append(1)
    return res


def _gen_dnf(points):
    res=''
    for point in points:
        for i in range(len(point)):
            if point[i]==0:
                res+=str(i)+"' "
            else:
                res+=str(i)+' '
        res+='+'
    return res


def _dnf2anf(s):
    cnf_list=[]
    clauses=s.split('+')
    for clause in clauses:
        tmp=[]
        res=[]
        clause=clause.split(' ')
        while '' in clause:
            clause.remove('')
        for var in clause:
            if var[len(var)-1]=="'":
                tmp.append([set([int(var[:len(var)-1])]),1])
            else:
                tmp.append([set([int(var)])])
        for i in range(len(tmp)):
            res=_mul(res,tmp[i])
        cnf_list+=res
    return _merge(cnf_list)


def _print_anf(x):
    ss=''
    for i in range(len(x)):
        if type(x[i])==set:
            for var in x[i]:
                ss+='x'+str(var)
            ss+=' + '
    if 1 not in x:
        ss=ss[:len(ss)-3]
    else:
        ss+='1'
    return ss

def gen_anf(ciphersbox):
    input_size=sbox.get_input_size_sbox(ciphersbox)
    output_size=sbox.get_output_size_sbox(ciphersbox)
    res=''
    posspoint=[[] for i in range(output_size)]
    for i in range(2**input_size):
        xbin=operation.num2binlist(input_size,i)
        ybin=operation.num2binlist(output_size,ciphersbox[i])
        for j in range(output_size):
            if ybin[j]==1:
                posspoint[j].append(xbin)
    for i in range(output_size):
        dnf=_gen_dnf(posspoint[i])
        anf=_dnf2anf(dnf)
        res+='y'+str(i)+' = '+_print_anf(anf)+'\n'
    return res