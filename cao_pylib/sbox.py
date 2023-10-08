#coding utf-8
import math
from cao_pylib import operation
import copy

def get_input_size_sbox(sbox):
    return int(math.log(len(sbox),2))


def get_output_size_sbox(sbox):
    sbox_set=set(sbox)
    return int(math.log(len(sbox_set),2))


def gen_ddt(sbox):
    input_size=get_input_size_sbox(sbox)
    output_size=get_output_size_sbox(sbox)
    ddt=[[0 for i in range(2**output_size)] for j in range(2**input_size)]
    for input_1 in range(2**input_size):
        for input_2 in range(2**input_size):
            output_1=sbox[input_1]
            output_2=sbox[input_2]
            ddt[input_1^input_2][output_1^output_2]+=1
    return ddt


def gen_ldt(sbox):
    input_size=get_input_size_sbox(sbox)
    output_size=get_output_size_sbox(sbox)
    ldt=[[0 for i in range(2**output_size)] for j in range(2**input_size)]

    for input_mask in range(2**input_size):
        for  output_mask in range(2**output_size):
            num0=0
            for input1 in range(2**input_size):
                strbin=bin(input_mask&input1)[2:]+bin(sbox[input1]&output_mask)[2:]
                sum_value=0
                for i in range(len(strbin)):
                    sum_value^=int(strbin[i])
                if sum_value==0:
                    num0+=1
            ldt[input_mask][output_mask]=abs(num0-int(2**input_size/2))
    return ldt


def poss_outdiff_of_fixed_inputdiff(sbox,input_diff):
    """output all possible out difference of sbox in a fixed input difference

    Args: ...

    return: a two dimensional list contain all possible output patterns. eg:[[0,0],[0,1]]  
    """
    tmp=[]
    output_size=get_output_size_sbox(sbox)
    ddt=gen_ddt(sbox)

    for i in range(2**output_size):
        if ddt[input_diff][i]!=0:
            tmp1=[]
            istr='0'*(output_size-(len(bin(i)[2:])))+bin(i)[2:]
            for j in range(output_size):
                tmp1.append(eval(istr[j]))
            tmp.append(tmp1)
    return tmp


def all_poss_diff_patterns(sbox):
    """
    return: a two dimensional list contain all possible [input,output] patterns.
    """
    temp = []
    input_size = get_input_size_sbox(sbox)
    for input_diff in range(2 ** input_size):
        possible_output_diff = poss_outdiff_of_fixed_inputdiff(sbox,input_diff)
        input_diff_pattern = [0 for i in range(input_size + 2 - len(bin(input_diff)))] + list(map(int, bin(input_diff)[2:]))
        if possible_output_diff !=[]:
            for output_diff_pattern in possible_output_diff:
                point = input_diff_pattern + output_diff_pattern
                temp.append(point)
    return temp


def all_imposs_diff_patterns(sbox):
    """
    return: a two dimensional list contain all impossible [input,output] patterns.
    """
    temp = []
    input_size = get_input_size_sbox(sbox)
    output_size = get_output_size_sbox(sbox)
    all_possible_diff_pattern = all_poss_diff_patterns(sbox)
    for input_diff in range(2 ** input_size):
        for output_diff in range(2 ** output_size):
            point = [0 for i in range(input_size + 2 - len(bin(input_diff)))] + list(map(int, bin(input_diff)[2:])) + \
                    [0 for i in range(output_size + 2 - len(bin(output_diff)))] + list(map(int, bin(output_diff)[2:]))
            if point not in all_possible_diff_pattern:
                temp.append(point)
    return temp    


def poss_outmask_of_fixed_inputmask(sbox,input_mask):
    """output all possible out mask of sbox in a fixed input mask

    Args: ...

    return: a two dimensional list contain all possible output patterns. eg:[[0,0],[0,1]]  
    """
    tmp=[]
    output_size=get_output_size_sbox(sbox)
    ldt=gen_ldt(sbox)

    for i in range(2**output_size):
        if ldt[input_mask][i]!=0:
            tmp1=[]
            istr='0'*(output_size-(len(bin(i)[2:])))+bin(i)[2:]
            for j in range(output_size):
                tmp1.append(eval(istr[j]))
            tmp.append(tmp1)
    return tmp


def all_poss_line_patterns(sbox):
    """
    return: a two dimensional list contain all possible [input,output] patterns.
    """
    temp = []
    input_size = get_input_size_sbox(sbox)
    for input_mask in range(2 ** input_size):
        possible_output_mask = poss_outmask_of_fixed_inputmask(sbox,input_mask)
        input_mask_pattern = [0 for i in range(input_size + 2 - len(bin(input_mask)))] + list(map(int, bin(input_mask)[2:]))
        if possible_output_mask !=[]:
            for output_mask_pattern in possible_output_mask:
                point = input_mask_pattern + output_mask_pattern
                temp.append(point)
    return temp


def all_imposs_line_patterns(sbox):
    """
    return: a two dimensional list contain all impossible [input,output] patterns.
    """
    temp = []
    input_size = get_input_size_sbox(sbox)
    output_size = get_output_size_sbox(sbox)
    all_possible_line_pattern = all_poss_line_patterns(sbox)
    for input_mask in range(2 ** input_size):
        for output_mask in range(2 ** output_size):
            point = [0 for i in range(input_size + 2 - len(bin(input_mask)))] + list(map(int, bin(input_mask)[2:])) + \
                    [0 for i in range(output_size + 2 - len(bin(output_mask)))] + list(map(int, bin(output_mask)[2:]))
            if point not in all_possible_line_pattern:
                temp.append(point)
    return temp


def trans_table(tab):
    """
    将分布表从概率变为差分重量形式
    """
    pro_set=set()
    for i in range(len(tab)):
        pro_set|=set(tab[i])
    
    flag=True
    for pro in pro_set:
        if pro!=0:
            if math.log(pro,2)-int(math.log(pro,2))!=0:
                flag=False
    assert flag,'some probability of table is not 2^(-n)'

    table=copy.deepcopy(tab)
    input_num=len(table)
    for i in range(len(table)):
        for j in range(len(table[0])):
            table[i][j]/=input_num
            if table[i][j]==0:
                table[i][j]='x'
            else:
                table[i][j]=int(math.log(table[i][j],2))*(-1)
    return table


def all_poss_diff_patterns_table_withPro(tab):
    """
    table中的值是差分重量
    差分重量和p的比特的关系：p0*max+...+pn
    """
    table=copy.deepcopy(tab)
    table=trans_table(table)
    poss_patterns=[]
    input_num=len(table)
    input_size=int(math.log(input_num,2))
    output_num=len(table[0])
    output_size=int(math.log(output_num,2))

    max_weight=0
    for i in range(input_num):
        for j in range(output_num):
            if table[i][j]!='x':
                max_weight=max(table[i][j],max_weight)
    p_num=int(math.log(max_weight,2))+1

    for i in range(input_num):
        for j in range(input_num):
            if table[i][j]!='x':
                poss_patterns.append(operation.num2binlist(input_size,i)+operation.num2binlist(output_size,j)+operation.num2binlist(p_num,table[i][j]))
    return poss_patterns


def all_imposs_diff_patterns_table_withPro(tab):
    table=copy.deepcopy(tab)
    all_poss_patterns=all_poss_diff_patterns_table_withPro(table)
    width=len(all_poss_patterns[0])
    all_imposs_patterns=[]

    for i in range(2**width):
        ibin_list=operation.num2binlist(width,i)
        if ibin_list not in all_poss_patterns:
            all_imposs_patterns.append(ibin_list)
    return all_imposs_patterns


def all_poss_diff_patterns_table(tab):
    # generate a two dimensional list cotains all [input,output] list,no probability
    table=copy.deepcopy(tab)
    input_num=len(table)
    input_size=int(math.log(input_num,2))
    output_num=len(table[0])
    output_size=int(math.log(output_num,2))

    poss_pattern=[]
    for i in range(len(table)):
        for j in range(len(table[0])):
            if table[i][j]!=0:
                poss_pattern.append(operation.num2binlist(input_size,i)+operation.num2binlist(output_size,j))
    return poss_pattern


def all_imposs_diff_patterns_table(tab):
    # generate a two dimensional list cotains all impossible [input,output] list,no probability
    table=copy.deepcopy(tab)
    all_poss_patterns=all_poss_diff_patterns_table(table)
    width=len(all_poss_patterns[0])
    all_imposs_patterns=[]

    for i in range(2**width):
        ibin_list=operation.num2binlist(width,i)
        if ibin_list not in all_poss_patterns:
            all_imposs_patterns.append(ibin_list)
    return all_imposs_patterns


def write_patterns_csv(patterns,filename):
    """
    将所有的patterns写入到csv文件中，变量名用1，2，3，4,....
    """
    point_num=len(patterns[0])
    s=''
    for i in range(point_num):
        s+=str(i)+','
    s+=',F\n'

    for i in range(len(patterns)):
        tmp=''
        for j in range(point_num):
            tmp+=str(patterns[i][j])+','
        tmp+=',1\n'
        s+=tmp

    fp=open(filename,'w')
    fp.write(s[:len(s)-1])
    fp.close()


def write_poss_pattern(poss_pattern,filename):
    """将差分\线性模式写入到文件中

    Args:
        poss_pattern: 差分\线性模式
        filename: 要写入的文件
    """
    fp=open(filename,'w+')
    string='['
    for pattern in poss_pattern:
        string+='['
        for var in pattern:
            string+=str(var)+','
        string=string[:len(string)-1]+'],'
    string=string[:len(string)-1]+']'
    fp.write(string)
    fp.close()


def ddt2truthtable(ddt):
    """
    transform ddt to truthtable,the probability add to output,eg:output:1101,p=3->1101 011
    """
    ddt=copy.deepcopy(ddt)
    ddt=trans_table(ddt)
    input_width=int(math.log2(len(ddt)))

    dist=set()
    for i in range(len(ddt)):
        for j in range(len(ddt[0])):
            if ddt[i][j]!='x':
                dist.add(ddt[i][j])
    p_width=int(math.log2(max(dist)))+1
    output_num=2**(input_width+p_width)
    truthtable=[[0 for i in range(output_num)] for i in range(2**input_width)]

    for i in range(len(ddt)):
        for j in range(len(ddt[0])):
            if ddt[i][j]!='x':
                pbin=operation.num2binlist(p_width,ddt[i][j])
                outbin=operation.num2binlist(input_width,j)+pbin
                output=operation.binlist2num(outbin)
                truthtable[i][output]=1
    return truthtable


def truthtable_reverse(tab):
    table=copy.deepcopy(tab)
    for i in range(len(table)):
        for j in range(len(table[0])):
            table[i][j]=table[i][j]^1
    return table


def all_poss_line_patterns_ldt_with_pro_cor(ldt):
    """
    根据线性分布表返回带概率的高维线性的点
    """
    input_size=int(math.log2(len(ldt)))
    output_size=int(math.log2(len(ldt[0])))
    res=copy.deepcopy(ldt)
    res1=[]

    p=[]
    for input_mask in range(2**input_size):
        for output_mask in range(2**output_size):
            pro=ldt[input_mask][output_mask]/(2**output_size)
            if pro!=0:
                wight=-math.log2(pro)
                if int(wight)-wight!=0:
                    print('err, the weight of linear is not a integer!')
                    return False
                else:
                    res[input_mask][output_mask]=int(wight)-1
                    if int(wight)-1 not in p:
                        p.append(int(wight)-1)
            else:
                res[input_mask][output_mask]=-1
    width=math.ceil(math.log2(len(p)))

    for input_mask in range(2**input_size):
        for output_mask in range(2**output_size):
            if res[input_mask][output_mask]>=0:
                res1.append(operation.num2binlist(input_size,input_mask)+operation.num2binlist(output_size,output_mask)+operation.num2binlist(width,res[input_mask][output_mask]))
    return res1
    
    



def get_p(in_diff,out_diff,truthtable):
    """
    return a binary expression of probability ,LSB
    """
    res=[]
    in_num=len(truthtable)

    pnum=int(len(truthtable[0])/len(truthtable))
    for i in range(pnum):
        res.append(truthtable[in_diff][out_diff*pnum+i])
    return res


def truthtable2ddt(table):
    input_num=len(table)
    ddt=[[0 for i in range(input_num)] for i in range(input_num)]
    for in_diff in range(input_num):
        for out_diff in range(input_num):
            p=get_p(in_diff,out_diff,table)
            for i in range(len(p)):
                if p[i]==0:
                    ddt[in_diff][out_diff]=input_num//2**(i)
    return ddt

def get_invert_sbox(s):
    invert_s=[0 for i in range(len(s))]
    for i in range(len(s)):
        invert_s[s[i]]=i
    return invert_s
    

