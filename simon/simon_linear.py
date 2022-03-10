from gurobipy import *
from cao_pylib import operation
from cao_pylib import sbox_model
import time

def var2value(x):
    res=[]
    varnum=int(len(x)/2)
    for i in range(varnum):
        if [x[i*2],x[i*2+1]]==[0,0]:
            res.append(0)
        if [x[i*2],x[i*2+1]]==[0,1]:
            res.append(1)
        if [x[i*2],x[i*2+1]]==[1,1]:
            res.append(2)
    return res

class Simon_L():
    def __init__(self,rounds,block_size,index,value):
        self.r=rounds
        self.block_size=block_size
        self.index=index
        if value==0:
            self.value=[0,0]
        if value==1:
            self.value=[0,1]
        
        self.m=Model()
        self.x=self.m.addVars([i for i in range((self.r+1)*self.block_size*2)],vtype=GRB.BINARY,name='x')
        self.a=self.m.addVars([i for i in range((self.r)*self.block_size)],vtype=GRB.BINARY,name='a')
        self.b=self.m.addVars([i for i in range((self.r)*self.block_size)],vtype=GRB.BINARY,name='b')
        self.c=self.m.addVars([i for i in range((self.r)*self.block_size)],vtype=GRB.BINARY,name='c')
        self.aout=self.m.addVars([i for i in range((self.r)*self.block_size)],vtype=GRB.BINARY,name='aout')

        constr=LinExpr()
        for i in range(self.block_size*2):
            constr+=self.x[i]
        self.m.addConstr(constr>=1)
        #self.m.addConstr(constr<=1)
    
    def gen_constr(self):
        for r in range(self.r):
            index=self.block_size*r
            x0=[]
            x1=[]
            a=[]
            b=[]
            c=[]
            aout=[]
            for i in range(self.block_size):
                x0.append([self.x[(index+i)*2],self.x[(index+i)*2+1]])
                x1.append([self.x[(index+i+self.block_size)*2],self.x[(index+i+self.block_size)*2+1]])
            for i in range(int(self.block_size/2)):
                a.append([self.a[index+i*2],self.a[index+i*2+1]])
                b.append([self.b[index+i*2],self.b[index+i*2+1]])
                c.append([self.c[index+i*2],self.c[index+i*2+1]])
                aout.append([self.aout[index+i*2],self.aout[index+i*2+1]])
            
            for i in range(int(self.block_size/2)):
                sbox_model.gen_model_from_ine(self.m,x0[i]+a[i]+b[i]+c[i]+x1[i+int(self.block_size/2)],'ine/_4xor_ine.txt')

            a=operation.roate_left(a,1)
            b=operation.roate_left(b,8)
            c=operation.roate_left(c,2)
            for i in range(int(self.block_size/2)):
                sbox_model.gen_model_from_ine(self.m,a[i]+b[i]+aout[i],'simon/br_ine.txt')
            
            
            for i in range(int(self.block_size/2)):
                self.m.addConstr(aout[i][0]==x0[i+int(self.block_size/2)][0])
                self.m.addConstr(aout[i][1]==x0[i+int(self.block_size/2)][1])
            for i in range(int(self.block_size/2)):
                self.m.addConstr(c[i][0]==x1[i][0])
                self.m.addConstr(c[i][1]==x1[i][1])
            for i in range(int(self.block_size/2)):
                self.m.addConstr(c[i][0]==aout[i][0])
                self.m.addConstr(c[i][1]==aout[i][1])
        
        index=self.index+self.r*self.block_size
        self.m.addConstr(self.x[index*2]==self.value[0])
        self.m.addConstr(self.x[index*2+1]==self.value[1])

def get_value(simon):
    xnum=len(simon.x)
    x=[0 for i in range(xnum)]

    for v in simon.m.getVars():
        if v.varName[0]=='x':
            x[eval(v.varName[2:len(v.varName)-1])]=int(v.x)
    
    x=var2value(x)

    for r in range(simon.r+1):
        s=''
        for i in range(simon.block_size):
            s+=str(x[r*simon.block_size+i])
            if i==int(simon.block_size/2)-1:
                s+=' '
        print(s)
    
    
    
for i in range(0,64,16):
    print(i)
    test=Simon_L(6,64,i,1)
    test.gen_constr()
    test.m.optimize()
    test.m.Params.OutputFlag=0
    if test.m.Status!=GRB.status.INFEASIBLE:
        get_value(test)