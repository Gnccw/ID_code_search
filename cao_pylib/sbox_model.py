"""
更新时间：2020-11-9

11-9：增加gen_ine_from_cnf

"""
import copy
import numpy
from gurobipy import *
from cao_pylib import ine
from cao_pylib import operation
from pyeda.inter import *

def gen_ine_from_cnf(cnf,best_ine):
    """
    从cnf生成不等式，并存入best_ine文件,文字的命名：0,1，2，3...
    """
    vars=set()
    all_clause=cnf[1:len(cnf)-1].split(')(')
    for clause in all_clause:
        liter=clause.split('+')
        for var in liter:
            if var[len(var)-1]=="'":
                vars=vars.union([eval(var[:len(var)-1])])
            else:
                vars=vars.union([eval(var)])
    ine_all=[]
    for clause in all_clause:
        liter=clause.split('+')
        const=-1
        tmp=[0 for i in range(len(vars)+2)]
        for index in vars:
            if str(index) in liter:
                tmp[index]=1
            if str(index)+"'" in liter:
                const+=1
                tmp[index]=-1
        tmp[len(vars)]=const
        ine_all.append(tmp)

    #最优不等式写入文件
    fp=open(best_ine,'w+')
    res=''
    for line in ine_all:
        s=ine.list2ine(line)
        s+='\n'
        res+=s
    fp.write(res[:len(res)-1])
    fp.close()


def min_ine_of_sbox(cipher_name,file_from_sage,imposs_patterns):
    """
    file_from_sage:由可能差分模式得到的不等式
    得到形如"a_0 * x_0 + a_1 * x_1 + ... + a_{n-1} * x_{n-1} + c >= 0"，用[...,c,0]来表示
    """
    output_filename='best_ine_of_'+cipher_name+'.txt'
    
    #读sage生成的不等式
    ine_all = []
    with open(file_from_sage, 'r') as f:
        for line in f:
            line_str = copy.deepcopy((line.replace('(', '').replace(') x + ', ', ').replace(') x - ', ', -').replace(
                ' >=', ',').replace('\n', '').replace(' ', '')))
            line_int = copy.deepcopy(list(map(int, line_str.split(','))))
            ine_all.append(line_int)
    all_imposs_patterns=imposs_patterns
    best_cutoff_ine=ine.best_cutoff_ine(ine_all,all_imposs_patterns)

    #最优不等式写入文件
    fp=open(output_filename,'w+')
    res=''
    for line in best_cutoff_ine:
        s=ine.list2ine(line)
        s+='\n'
        res+=s
    fp.write(res[:len(res)-1])
    fp.close()

def gen_model_from_ine(m,x,ine):
    """
    generate a MILP model from the inequalities in the file "ine"
    """

    #read inequalities from file
    ine_all = []
    with open(ine, 'r') as f:
        for line in f:
            line_str = copy.deepcopy((line.replace('(', '').replace(') x + ', ', ').replace(') x  - ', ', -').replace(
                ' >=', ',').replace('\n', '').replace(' ', '')))
            line_int = copy.deepcopy(list(map(int, line_str.split(','))))
            ine_all.append(line_int)
    for ine in ine_all:
        constr=LinExpr()
        for i in range(len(ine)-2):
            constr+=ine[i]*x[i]
        constr+=(ine[len(ine)-2])
        m.addConstr(constr>=ine[len(ine)-1])


def test_ine_of_sbox(ine,posspatterns,imposspatterns):
    """
    ine:存储不等式的文件：-1, 1, 0, 0：-x0+x1+0>=0
    测试生成的不等式是否正确
    """
    flag=0
    for point in posspatterns:
        m=Model()
        var_num=len(posspatterns)
        x=m.addVars([i for i in range(var_num)],vtype=GRB.BINARY,name='x')
        gen_model_from_ine(m,x,ine)
        for i in range(len(point)):
            m.addConstr(x[i]==point[i])
        
        m.Params.OutputFlag=0
        m.optimize()

        if m.status==GRB.Status.INFEASIBLE:
            flag=1
            print('error:poss is imposs',point) 

    for point in imposspatterns:
        m=Model()
        var_num=len(imposspatterns)
        x=m.addVars([i for i in range(var_num)],vtype=GRB.BINARY,name='x')
        gen_model_from_ine(m,x,ine)
        for i in range(len(point)):
            m.addConstr(x[i]==point[i])
        
        m.Params.OutputFlag=0
        m.optimize()
        
        if m.status!=GRB.Status.INFEASIBLE:
            flag=1
            print('error:imposs is poss:',point)
    if flag==0:
        print('------------ine is right!----------')

def num_of_equalities(inefile):
    """
    """
    ine_all=[]
    resss=[]
    with open(inefile, 'r') as f:
        for line in f:
            line_str = copy.deepcopy((line.replace('(', '').replace(') x + ', ', ').replace(') x - ', ', -').replace(
                ' >=', ',').replace('\n', '').replace(' ', '')))
            line_int = copy.deepcopy(list(map(int, line_str.split(','))))
            ine_all.append(line_int)
    var_num=len(ine_all[0])-2

    for ine in ine_all:
        res=0
        coeffei=numpy.array(ine[:len(ine)-2])
        for i in range(2**var_num):
            ibin=numpy.array(operation.num2binlist(var_num,i))
            if numpy.dot(coeffei,ibin)+ine[len(ine)-2]==0:
                res+=1
        resss.append([ine,res])
    return resss


def _expr2cnf(s,num_of_bits):
    s=s.replace('Or', '')
    s=s.split('A')
    n=num_of_bits
    res=''
    clauses=[]

    for ss in s:
        if 'd' in ss:
            clauses.append(ss)
        elif 'x' in ss:
            clauses.append(ss)
    
    for cla in clauses:
        liter=[]
        cla=cla.split(',')
        for i in range(len(cla)):
            if 'x' in cla[i]:
                liter.append(cla[i])
        cnf_cla='('
        for i in range(len(liter)):
            tmp=''
            for j in range(len(liter[i])):
                if liter[i][j].isdigit():
                    tmp+=liter[i][j]
            tmp=str(num_of_bits-1-eval(tmp))
            if '~' not in liter[i]:
                tmp+="'"
            cnf_cla+=tmp+'+'
        cnf_cla=cnf_cla[:len(cnf_cla)-1]+')'
        res+=cnf_cla
    return res


def gen_cnf_from_point(table,ine_file):
    """
    利用pyeda模块下的espresso算法得到可行point点的cnf表达,然后将表达式转换为不等式并存入ine_file
    """
    tts=''
    ones=[]
    for point in table:
        ones.append(operation.binlist2num(point))
    for i in range(2**len(table[0])):
        if i in ones:
            tts+='0'
        else:
            tts+='1'
    X=ttvars('x',len(table[0]))
    f=truthtable(X,tts)
    fm=espresso_tts(f)[0]
    cnf=_expr2cnf(str(fm), len(table[0]))
    gen_ine_from_cnf(cnf, ine_file)
