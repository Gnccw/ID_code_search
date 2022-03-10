
from typing import List
from gurobipy import *
from cao_pylib import operation
from cao_pylib import sbox_model

class Twine():
    def __init__(self,rounds,block_size,index,value):
        self.p=[5,0,1,4,7,12,3,8,13,6,9,2,15,10,11,14]
        self.rounds=rounds
        self.block_size=16
        self.index=index
        self.value=operation.num2binlist(3,value)

        self.m=Model()
        self.x=self.m.addVars([i for i in range(self.block_size*(self.rounds+1)*3)],vtype=GRB.BINARY,name='x')
        self.fout=self.m.addVars([i for i in range((int(self.block_size*1.5))*self.rounds)],vtype=GRB.BINARY,name='fout')

        input_diff=LinExpr()
        for i in range(48):
            input_diff+=self.x[i]
        self.m.addConstr(input_diff>=1)

    def gen_constr(self):
        for r in range(self.rounds):
            x0=[]
            x1=[]
            fout=[]

            for i in range(16):
                x0.append([self.x[(r*16+i)*3],self.x[(r*16+i)*3+1],self.x[(r*16+i)*3+2]])
                x1.append([self.x[((r+1)*16+i)*3],self.x[((r+1)*16+i)*3+1],self.x[(r*16+i)*3+2]])
            for i in range(8):
                fout.append([self.fout[(r*8+i)*3],self.fout[(r*8+i)*3+1],self.fout[(r*8+i)*3+2]])
            self.m.update()
            
            for i in range(8):
                sbox_model.gen_model_from_ine(self.m,x0[i*2]+fout[i],'twine//ine//F_ine.txt')
                sbox_model.gen_model_from_ine(self.m,fout[i]+x0[i*2+1]+x1[self.p[i*2+1]],'twine//ine//xor_ine.txt')
                self.m.addConstr(x0[i*2][0]==x1[self.p[i*2]][0])
                self.m.addConstr(x0[i*2][1]==x1[self.p[i*2]][1])
                self.m.addConstr(x0[i*2][2]==x1[self.p[i*2]][2])
            
        index=self.index+self.rounds*16
        self.m.addConstr(self.x[index*3]==self.value[0])
        self.m.addConstr(self.x[index*3+1]==self.value[1])
        self.m.addConstr(self.x[index*3+2]==self.value[2])


class Twine_inv():
    def __init__(self,rounds,block_size,index,value):
        self.p=[5,0,1,4,7,12,3,8,13,6,9,2,15,10,11,14]
        self.rounds=rounds
        self.block_size=16
        self.index=index
        self.value=operation.num2binlist(3,value)

        self.m=Model()
        self.x=self.m.addVars([i for i in range(self.block_size*(self.rounds+1)*3)],vtype=GRB.BINARY,name='x')
        self.fout=self.m.addVars([i for i in range((int(self.block_size*1.5))*self.rounds)],vtype=GRB.BINARY,name='fout')

        input_diff=LinExpr()
        for i in range(48):
            input_diff+=self.x[i]
        self.m.addConstr(input_diff>=1)
    
    def gen_constr(self):
        for r in range(self.rounds):
            x0=[]
            x1=[]
            fout=[]

            for i in range(16):
                x0.append([self.x[(r*16+i)*3],self.x[(r*16+i)*3+1],self.x[(r*16+i)*3+2]])
                x1.append([self.x[((r+1)*16+i)*3],self.x[((r+1)*16+i)*3+1],self.x[(r*16+i)*3+2]])
            for i in range(8):
                fout.append([self.fout[(r*8+i)*3],self.fout[(r*8+i)*3+1],self.fout[(r*8+i)*3+2]])
            self.m.update()

            for i in range(8):
                sbox_model.gen_model_from_ine(self.m,x1[i*2]+fout[i],'twine//ine//F_ine.txt')
                sbox_model.gen_model_from_ine(self.m,fout[i]+x0[self.p[i*2+1]]+x1[i*2+1],'twine//ine//xor_ine.txt')
                self.m.addConstr(x1[i*2][0]==x0[self.p[i*2]][0])
                self.m.addConstr(x1[i*2][1]==x0[self.p[i*2]][1])
                self.m.addConstr(x1[i*2][2]==x0[self.p[i*2]][2])
            
        index=self.index+self.rounds*16
        self.m.addConstr(self.x[index*3]==self.value[0])
        self.m.addConstr(self.x[index*3+1]==self.value[1])
        self.m.addConstr(self.x[index*3+2]==self.value[2])

def get_var(test):
    x=[0 for i in range((test.rounds+1)*48)]
    f=[0 for i in range(test.rounds*24)]
    for v in test.m.getVars():
        if v.varName[0]=='x':
            x[eval(v.varName[2:len(v.varName)-1])]=int(v.x)
    for v in test.m.getVars():
        if v.varName[0]=='f':
            f[eval(v.varName[5:len(v.varName)-1])]=int(v.x)
    
    ss=[]
    for i in range(test.rounds+1):
        tmp=''
        for j in range(16):
            if j%2==0:
                tmp+=' '
            tmp+=str(operation.binlist2num([x[(i*16+j)*3],x[(i*16+j)*3+1],x[(i*16+j)*3+2]]))
        ss.append(tmp)
    
    for i in range(test.rounds):
        tmp=''
        for j in range(8):
            tmp+=' '
            tmp+=str(operation.binlist2num([f[(i*8+j)*3],f[(i*8+j)*3+1],f[(i*8+j)*3+2]]))
        ss[i]+=tmp
    
    for i in ss:
        print(i)


for i in range(16):
    test=Twine_inv(7,16,i,0)
    test.gen_constr()
    test.m.Params.OutputFlag=0
    test.m.optimize()

    if test.m.Status!=GRB.status.INFEASIBLE:
        print('y',i)
        #get_var(test)
        #break
