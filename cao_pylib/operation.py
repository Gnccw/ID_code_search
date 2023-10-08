# -*- coding: utf-8 -*-

"""
更新时间：2020-11-7

包括以下函数：
1.roate_right(a,n)
2.roate_left(a,n)
3.num2binlist(width,num)
4.binlist2num(binlist)
5.gen_seq_weight(width,weight)
6.list_mul(x,y)
"""
def roate_right(a,n):
    """
    将列表a中的元素循环右移n位，输出新的列表
    """
    tmp=[]
    width=len(a)
    for i in range(width):
        tmp.append(a[(i+width-n)%width])
    return tmp


def roate_left(a,n):
    """
    将列表中的元素循环左移n位，输出新的列表
    """
    tmp=[]
    width=len(a)
    for i in range(width):
        tmp.append(a[(i+width+n)%width])
    return tmp


def num2binlist(width,num):
    '''
    生成num的width-bit的列表[x0,x1,x2,x3],x0为最高位
    '''
    tmp=[]
    num_str='0'*(width+2-len(bin(num)))+bin(num)[2:]
    for i in range(width):
        tmp.append(eval(num_str[i]))
    return tmp

def binlist2num(binlist):
    '''
    将形为[x0,x1,x2,x3]列表变为数字，x0为最高位（权重最大）
    '''
    width=len(binlist)
    res=0
    for i in range(width):
        res+=binlist[i]*(2**(width-i-1))
    return res


def gen_seq_weight(width,weight):
    '''
    生成宽度为width,汉明重量为weight的所有的数
    '''
    seq=[]
    for i in range(2**width):
        res=0
        ibin=num2binlist(width,i)
        for j in ibin:
            if j==1:
                res+=1
        if res==weight:
            seq.append(i)
    return seq


def list_mul(x,y):
    '''
    两个列表x,y的对应位置的乘积，值放在一个新的列表中
    '''
    res=[]
    for i in range(len(x)):
        res.append(x[i]*y[i])
    return res

def list2str(x):
    """
    将list中的元素变为字符，并用逗号隔开
    """
    s=''
    for i in x:
        s+=str(i)+','
    return s[:len(s)-1]

def str2list(s):
    """
    将用逗号隔开的str变为list
    """
    res=[]
    string=s.split(',')
    for var in string:
        res.append(eval(var))
    return res


def write_d2list(table,filename):
    """
    将二维数组写入到文件中，每行存一个一维数组，元素用逗号隔开
    """
    s=''
    for row in table:
        tmp=''
        for element in row:
            tmp+=str(element)+','
        s+=tmp[:len(tmp)-1]
        s+='\n'
    fp=open(filename,'w')
    fp.write(s[:len(s)-1])
    fp.close()


def read_d2list(filename):
    """
    读取一个文件中的二维数组
    """
    fp=open(filename,'r')
    s=fp.read()
    fp.close()

    res=[]
    s=s.split('\n')
    for row in s:
        tmp=[]
        row=row.split(',')
        for element in row:
            tmp.append(eval(element))
        res.append(tmp)
    return res