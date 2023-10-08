
"""
更新时间：2020-11-7

生成差分分布表的某一个点的父集
包括以下函数：
1.gen_non_choice(n)
2.gen_index_begin(n)
3.clause2index(clause,n,index_begin,non_choice)
4.index2clause(pos,n,index_begin,non_choice)
5.gen_father_set(x,n)
6.gen_father_index(x,n)
"""

import operation
import scipy.special as sp
import copy


def gen_non_choice(n):
    '''
    n: 表示总的变量的个数
    non_choice: 非选择，non_choice[1]=[1,2,4,8] 表示有一个没有选择
    '''
    non_choice=[[0]]
    for i in range(1,n):
        non_choice.append(operation.gen_seq_weight(n,i))
    return non_choice


def gen_index_begin(n):
    #n表示总的变量的个数
    index=[0,0]
    for i in range(1,n+1):
        index.append(int(sp.comb(n,i)*2**i)+index[len(index)-1])
    return index


def clause2index(clause,n,index_begin,non_choice):
    """
    将clause转变成下标，clause的表示（-1，1，0，0）：~a+b,n=4
    n表示总的变量的个数
    """
    row=0   #表示已选择的个数
    value_list=[]
    non_choice_value_list=[]
    for i in range(len(clause)):
        if clause[i]!=0:
            row+=1
            value_list.append(int((clause[i]+1)/2))
            non_choice_value_list.append(0)
        else:
            non_choice_value_list.append(1)
    non_choice_value=operation.binlist2num(non_choice_value_list)
    value=operation.binlist2num(value_list)
    for i in range(len(non_choice[n-row])):
        if non_choice[n-row][i]==non_choice_value:
            pos=i
            break
    return index_begin[row]+(2**row)*pos+value


def index2clause(pos,n,index_begin,non_choice):
    """
    将pos转变成clause，clause的表示（-1，1，0，0）：~a+b,n=4
    n表示总的变量的个数
    row:取的变量的个数
    """
    res=[0 for i in range(n)]
    for i in range(len(index_begin)-1,-1,-1):
        if pos>=index_begin[i]:
            row=i
            break
    pos=pos-index_begin[row]
    #print('row',row)
    non_choice_col=pos//2**row
    #print('non_choice',non_choice_col)
    choice_value=pos%(2**row)
    non_choice_value=non_choice[n-row][non_choice_col]
    non_choice_bin=operation.num2binlist(n,non_choice_value)
    choice_bin=operation.num2binlist(row,choice_value)
    for i in range(row):
        if choice_bin[i]==0:
            choice_bin[i]=-1
    for i in range(n):
        non_choice_bin[i]^=1
    for i in range(row):
        for j in range(n):
            if non_choice_bin[j]==1:
                res[j]=choice_bin[i]
                non_choice_bin[j]=0
                break
    return res


def gen_father_set(x,n):
    """
    生成逻辑表达式代表的值的父集，输出的集合是逻辑表达式
    """
    father_set=[]
    xbin=operation.num2binlist(n,x)
    for i in range(n):
        if xbin[i]==0:
            xbin[i]=-1
    for i in range(2**n):
        y=operation.num2binlist(n,i)
        y=operation.list_mul(xbin,y)
        flag=0
        for j in y:
            if j!=0:
                flag=1
                break
        if flag==1:
            father_set.append(y)
    return father_set


def gen_father_index(x,n):
    """
    生成x的父集的所有的变量的下标
    Args:
        x:the position of truthtable
        n:number of var
    """
    index_begin=gen_index_begin(n)
    non_choice=gen_non_choice(n)
    res=[]
    for x in gen_father_set(x,n):
        res.append(clause2index(x,n,index_begin,non_choice))
    return res