'''
from knot.knot import Knot
from knot.knot import Inv_knot
from knot.knot import get_value
from ascon.ascon import Ascon
from ascon.ascon import Ascon_inv
from Lea.lea import Lea
from Lea.lea import Lea_inv

from simon.simon import Simon
from simon.simon import Simon_inv
from rectangle.rectangle import Rectangle
from rectangle.rectangle import Rectangle_inv
'''
#from present.present import present_inv
#from present.present import present
from simon.simon import Simon
from simon.simon import Simon_inv
from gurobipy import *
from simeck.simeck import Simeck
from simeck.simeck import Simeck_inv
from hight.hight import Hight
from hight.hight import Hight_inv
import copy
import os

#rounds存在，返回rounds的矛盾
def longest_ID(cla0,cla1,block_size,rounds=''):
    ciphername=cla0.lower()
    cla0=globals() [cla0]
    cla1=globals() [cla1]
    longestID=0
    constradict=[]

    M0_res=[]
    M1_res=[]

    if not os.path.exists(ciphername+'//result//'+ciphername+'_'+str(block_size)+'_forward.txt'):
        flag=1
        for r in range(1,100):
            if flag==1:
                flag=0
                for index in range(block_size):
                    for value in range(2):
                        print(r,index,value)
                        M0=cla0(r,block_size,index,value)
                        M0.gen_constr()
                        M0.m.Params.OutputFlag=0
                        #M0.m.Params.MIPFocus=1
                        M0.m.optimize()
                        if M0.m.Status!=GRB.status.INFEASIBLE:
                            M0_res.append([r,index,value])
                            flag=1
                            print('y')
        s=''
        for var in M0_res:
            s+=str(var[0])+','
            s+=str(var[1])+','
            s+=str(var[2])+'\n'
        fp=open(ciphername+'//result//'+ciphername+'_'+str(block_size)+'_forward.txt','w')
        fp.write(s)
        fp.close()
    else:
        fp=open(ciphername+'//result//'+ciphername+'_'+str(block_size)+'_forward.txt','r')
        s=fp.read()
        fp.close()
        s=s.split('\n')[:len(s.split())]
        for var in s:
            var=var.split(',')
            tmp=[]
            for i in range(3):
                tmp.append(eval(var[i]))
            M0_res.append(tmp)


    if not os.path.exists(ciphername+'//result//'+ciphername+'_'+str(block_size)+'_backward.txt'):
        flag=1
        for r in range(1,100):
            if flag==1:
                flag=0
                for index in range(block_size):
                    for value in range(2):
                        print(r,index,value)
                        M1=cla1(r,block_size,index,value)
                        M1.gen_constr()
                        M1.m.Params.OutputFlag=0
                        #M1.m.Params.MIPFocus=1
                        M1.m.optimize()
                        if M1.m.Status!=GRB.status.INFEASIBLE:
                            M1_res.append([r,index,value])
                            flag=1
                            print('y')
        s=''
        for var in M1_res:
            s+=str(var[0])+','
            s+=str(var[1])+','
            s+=str(var[2])+'\n'
        fp=open(ciphername+'//result//'+ciphername+'_'+str(block_size)+'_backward.txt','w')
        fp.write(s)
        fp.close()
    else:
        fp=open(ciphername+'//result//'+ciphername+'_'+str(block_size)+'_backward.txt','r')
        s=fp.read()
        fp.close()
        s=s.split('\n')[:len(s.split())]
        for var in s:
            var=var.split(',')
            tmp=[]
            for i in range(3):
                tmp.append(eval(var[i]))
            M1_res.append(tmp)


    if rounds:
        for v in M0_res:
            if [rounds-v[0],v[1],v[2]^1] in M1_res:
                constradict.append(v)
        s=''
        for var in constradict:
            s+=str(var[0])+','
            s+=str(var[1])+','
            s+=str(var[2])+'\n'
        fp=open(ciphername+'//result//'+ciphername+'_'+str(block_size)+'_r'+str(rounds)+'_constradict.txt','w')
        fp.write(s)
        fp.close()
    else:
        dictM0={}
        dictM1={}
        for v in M0_res:
            dictM0[str(v[1])+'_'+str(v[2])]=v[0]
        for v in M1_res:
            dictM1[str(v[1])+'_'+str(v[2]^1)]=v[0]
        for key in dictM0.keys():
            if key in dictM1.keys():
                if dictM1[key]+dictM0[key]>longestID:
                    longestID=dictM1[key]+dictM0[key]
        return longestID
    

def var2value(v):
    res=[]
    for i in range(int(len(v)/2)):
        if [v[i*2],v[i*2+1]]==[0,0]:
            res.append(0)
        if [v[i*2],v[i*2+1]]==[0,1]:
            res.append(1)
        if [v[i*2],v[i*2+1]]==[1,1]:
            res.append(2)
    return res


#返回2*block_size个值的列表
def get_inputdiff(cla):
    res=[]
    x=[0 for i in range(len(cla.x))]
    for v in cla.m.getVars():
        if v.varName[0]=='x' and v.varName[1]=='[':
            x[eval(v.varName[2:len(v.varName)-1])]=int(v.x)
  
    for i in range(cla.block_size*2):
        res.append(x[i])
    return res



#返回2*block_size个值的列表
def get_output_diff(cla):
    res=[]
    x=[0 for i in range(len(cla.x))]
    for v in cla.m.getVars():
        if v.varName[0]=='s' and v.varName[1]=='o':
            x[eval(v.varName[5:len(v.varName)-1])]=int(v.x)
    for i in range(cla.block_size*2):
        res.append(x[i+(cla.r-1)*cla.block_size*2])
    return res


#给2*block_size个输入var添加约束
def add_input_constr(cla,inputdiff):
    constr=LinExpr()
    const=1
    for i in range(len(inputdiff)):
        if inputdiff[i]==1:
            constr-=cla.x[i]
            const-=1
        else:
            constr+=cla.x[i]
    cla.m.addConstr(constr>=const)


#给2*block_size个输入var赋值
def add_input_value(cla,inputdiff):
    for i in range(len(inputdiff)):
        cla.m.addConstr(cla.x[i]==inputdiff[i])


def list2str(inputlist):
    s=''
    for i in range(int(len(inputlist)/2)):
        if [inputlist[i*2],inputlist[i*2+1]]==[0,0]:
            s+='0'
        if [inputlist[i*2],inputlist[i*2+1]]==[0,1]:
            s+='1'
        if [inputlist[i*2],inputlist[i*2+1]]==[1,1]:
            s+='2'
    return s


def var2value(v):
    res=[]
    for i in range(int(len(v)/2)):
        if [v[i*2],v[i*2+1]]==[0,0]:
            res.append(0)
        if [v[i*2],v[i*2+1]]==[0,1]:
            res.append(1)
        if [v[i*2],v[i*2+1]]==[1,1]:
            res.append(2)
    return res


def value2var(v):
    res=[]
    for i in v:
        if i==0:
            res+=[0,0]
        if i==1:
            res+=[0,1]
        if i==2:
            res+=[1,1]
    return res


def all_poss_bin_list(a):
    res=[[]]
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
    return res


def num_value(x):
    res=0
    x=var2value(x)
    for i in range(len(x)):
        if x[i]==2:
            res+=1
    if 1 in x:
        return 2**res
    else:
        return 2**res-1


#inputdiff是2*block_size个值的列表
def add_input_constr_with2(cla,inputdiff):
    num=0
    tmp=var2value(inputdiff)
    inputdiffs=all_poss_bin_list(tmp)
    for i in inputdiffs:
        a=value2var(i)
        add_input_constr(cla,a)


def printf(s):
    res=''
    for i in range(int(len(s)/2)):
        if [s[i*2],s[i*2+1]]==[0,0]:
            res+='0'
        if [s[i*2],s[i*2+1]]==[0,1]:
            res+='1'
        if [s[i*2],s[i*2+1]]==[1,1]:
            res+='2'
    print(res)

#输出给定输入下的r轮的所有不可能输出差分
def all_outputs_given_input(inputs,cla0,cla1,constradict,block_size,rounds,rounds_M):
    print(inputs)
    outputnum=0
    all_outputs_res=[]
    all_outputs=[]
    middle_bits=[]
    for r in range(rounds_M[0],rounds_M[1]+1):
        M0=cla0(r,block_size,1,1)
        M0.gen_constr()
        M0.m.update()
        M0.m.remove(M0.m.getConstrByName('output0'))
        M0.m.remove(M0.m.getConstrByName('output1'))
        M0.m.update()
        add_input_value(M0,inputs)
        #M0.m.Params.OutputFlag=0
        M0.m.optimize()
        if M0.m.Status!=GRB.status.INFEASIBLE:
            outputdiff=get_output_diff(M0)
            outputdiff_value=var2value(outputdiff)
        for i in range(block_size):
            if [r,i,outputdiff_value[i]] in constradict:
                middle_bits.append([r,i,outputdiff_value[i]])
    print(middle_bits)
    
    flag=1
    for v in middle_bits:
        if flag==1:
            print(v)
            M1=cla1(rounds-v[0],block_size,v[1],v[2]^1)
            #print([rounds-v[0],block_size,v[1],v[2]^1])
            M1.gen_constr()
            M1.m.Params.OutputFlag=0

            if len(all_outputs)!=0:
                for i in range(len(all_outputs)):
                    add_input_constr(M1,all_outputs[i])
            
            M1.m.optimize()
            while M1.m.Status!=GRB.status.INFEASIBLE:
                if int(M1.m.getObjective().getValue())>=15:
                    input_diff_M1=get_inputdiff(M1)
                    all_outputs_res.append(input_diff_M1)
                    outputnum+=num_value(input_diff_M1)
                    flag=0
                    break 
                input_diff_M1=get_inputdiff(M1)
                all_outputs_res.append(input_diff_M1)
                tmp=var2value(input_diff_M1)
                tmp=all_poss_bin_list(tmp)
                
                sss=''
                for i in range(64):
                    sss+=str(tmp[0][i])
                print(sss)
                print(len(all_outputs))
                for i in range(len(tmp)):
                    all_outputs.append(value2var(tmp[i]))
                add_input_constr_with2(M1,input_diff_M1)
                outputnum+=num_value(input_diff_M1)

                M1.m.update()
                M1.m.optimize()
            
    return [all_outputs_res,outputnum]
    

def search_all_ID(cla0,cla1,block_size,rounds):
    ciphername=cla0
    #------------------read constradict from file--------------------------------
    constradict_file=ciphername.lower()+'//result//'+ciphername.lower()+'_'+str(block_size)+'_r'+str(rounds)+'_constradict.txt'
    fp=open(constradict_file,'r')
    s=fp.read()
    fp.close()

    constradict=[]
    s=s.split('\n')[:len(s.split())]
    for i in range(len(s)):
        tmp=s[len(s)-1-i].split(',')
        for j in range(3):
            tmp[j]=eval(tmp[j])
        print(tmp)
        constradict.append(tmp)
    #------------------read constradict from file--------------------------------

    cla0=globals() [cla0]
    cla1=globals() [cla1]

    rounds_M0=[100,0]
    for i in range(len(constradict)):
        if constradict[i][0]<rounds_M0[0]:
            rounds_M0[0]=constradict[i][0]
        if constradict[i][0]>rounds_M0[1]:
            rounds_M0[1]=constradict[i][0]
    
    all_inputs=[]
    input_output_num=[]
    for v in constradict:
        print('--------------',v)
        if v[0]<=165:
            print(v)
            #--------------construct M0 and add all_inputs constraint-------------
            M0=cla0(v[0],block_size,v[1],v[2])
            if len(all_inputs)!=0:
                for inputs in all_inputs:
                    add_input_constr(M0,inputs)
            M0.gen_constr()
            M0.m.Params.OutputFlag=0
            M0.m.optimize()
            #--------------construct M0 and add all_inputs constraint-------------
            while M0.m.Status!=GRB.status.INFEASIBLE:
                s='inputdiff:'
                inputdiff_M0=get_inputdiff(M0)
                all_inputs.append(inputdiff_M0)

                s+=list2str(inputdiff_M0)+'\n'
                tmp=all_outputs_given_input(inputdiff_M0,cla0,cla1,constradict,block_size,rounds,rounds_M0)
                for tmpvar in tmp[0]:
                    s+='          '+list2str(tmpvar)+'\n'
                input_output_num.append(tmp[1])
                add_input_constr(M0,inputdiff_M0)
                M0.m.update()
                M0.m.optimize()
                file_name=ciphername.lower()+'//result//'+ciphername.lower()+'_'+str(block_size)+'_all_imposs_diff_r'+str(rounds)+'.txt'
                fp=open(file_name,'a')
                fp.write(s+str(input_output_num)+'\n')
                fp.close()

'''
test=Knot(8,384,288,0)
test.gen_constr()
test.m.optimize()
get_value(test,384)
test=Inv_knot(8,384,288,1)
test.gen_constr()
test.m.optimize()
get_value(test,384)
'''
#search_all_ID('Gift','Gift_inv',64,7)
#print(longest_ID('Hight','Hight_inv',64,13))
