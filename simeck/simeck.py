from gurobipy import *
from cao_pylib import operation
from cao_pylib import sbox_model

class Simeck():
    def __init__(self,rounds,blocksize,index,value):
        self.r=rounds
        self.block_size=blocksize
        self.index=index
        if value==1:
            self.value=[0,1]
        if value==0:
            self.value=[0,0]
        
        self.m=Model()
        self.x=self.m.addVars([i for i in range(self.block_size*(self.r+1)*2)],vtype=GRB.BINARY,name='x')
        self.aout=self.m.addVars([i for i in range(self.r*self.block_size)],vtype=GRB.BINARY,name='aout')

        
        #输入差分不为0
        constr=LinExpr()
        for i in range(self.block_size*2):
            constr+=self.x[i]
        self.m.addConstr(constr>=1)

    def gen_constr(self):
        for r in range(self.r):
            index_start=r*self.block_size*2
            aindex_start=r*self.block_size
            x0=[]
            x1=[]
            x2=[]
            x3=[]
            x4=[]
            x5=[]
            a=[]
            for i in range(int(self.block_size/2)):
                tmp=[self.x[index_start+2*i],self.x[index_start+2*i+1]]
                tmp1=[self.x[index_start+self.block_size+2*i],self.x[index_start+self.block_size+2*i+1]]
                tmp2=[self.x[index_start+self.block_size*2+2*i],self.x[index_start+self.block_size*2+2*i+1]]
                tmp3=[self.aout[aindex_start+i*2],self.aout[aindex_start+i*2+1]]
                x0.append(tmp)
                x1.append(tmp1)
                x5.append(tmp2)
                a.append(tmp3)
            x2=operation.roate_left(x0,0)
            x3=operation.roate_left(x0,5)
            x4=operation.roate_left(x0,1)
            self.m.update()
            
            for i in range(int(self.block_size/2)):
                sbox_model.gen_model_from_ine(self.m,x2[i]+x3[i]+a[i],'ine//2and_ine.txt')
                sbox_model.gen_model_from_ine(self.m,a[i]+x1[i]+x4[i]+x5[i],'ine//3xor_ine.txt')
            
            for i in range(int(self.block_size/2)):
                self.m.addConstr(self.x[index_start+i*2]==self.x[index_start+3*self.block_size+2*i])
                self.m.addConstr(self.x[index_start+i*2+1]==self.x[index_start+3*self.block_size+2*i+1])
    
        self.m.addConstr(self.x[self.r*self.block_size*2+self.index*2]==self.value[0],name='output0')
        self.m.addConstr(self.x[self.r*self.block_size*2+self.index*2+1]==self.value[1],name='output1')
        
class Simeck_inv():
    def __init__(self,rounds,blocksize,index,value):
        self.r=rounds
        self.block_size=blocksize
        self.index=index
        if value==1:
            self.value=[0,1]
        if value==0:
            self.value=[0,0]
        
        self.m=Model()
        self.x=self.m.addVars([i for i in range(self.block_size*(self.r+1)*2)],vtype=GRB.BINARY,name='x')
        self.aout=self.m.addVars([i for i in range(self.r*self.block_size)],vtype=GRB.BINARY,name='aout')
        self.d=self.m.addVars([i for i in range(self.block_size)],vtype=GRB.BINARY,name='d')

        #输入差分不为0
        constr=LinExpr()
        for i in range(self.block_size*2):
            constr+=self.x[i]
        self.m.addConstr(constr>=1)
        
         #输入差分为2的个数
        constr=LinExpr()
        for i in range(self.block_size):
            constr+=self.d[i]
        self.m.setObjective(constr,GRB.MAXIMIZE)
        
        for i in range(self.block_size):
            self.m.addConstr(self.x[i*2]>=self.d[i])
            self.m.addConstr(self.x[i*2+1]>=self.d[i])
            self.m.addConstr(self.x[i*2]+self.x[i*2+1]-self.d[i]<=1)
    def gen_constr(self):
        for r in range(self.r):
            index_start=r*self.block_size*2
            aindex_start=r*self.block_size
            x0=[]
            x1=[]
            x2=[]
            x3=[]
            x4=[]
            x5=[]
            a=[]
            for i in range(int(self.block_size/2)):
                tmp=[self.x[index_start+2*i],self.x[index_start+2*i+1]]
                tmp1=[self.x[index_start+self.block_size+2*i],self.x[index_start+self.block_size+2*i+1]]
                tmp2=[self.x[index_start+self.block_size*2+2*i],self.x[index_start+self.block_size*2+2*i+1]]
                tmp3=[self.aout[aindex_start+i*2],self.aout[aindex_start+i*2+1]]
                x0.append(tmp)
                x1.append(tmp1)
                x5.append(tmp2)
                a.append(tmp3)
            x2=operation.roate_left(x0,0)
            x3=operation.roate_left(x0,5)
            x4=operation.roate_left(x0,1)
            self.m.update()
            
            for i in range(int(self.block_size/2)):
                sbox_model.gen_model_from_ine(self.m,x2[i]+x3[i]+a[i],'ine//2and_ine.txt')
                sbox_model.gen_model_from_ine(self.m,a[i]+x1[i]+x4[i]+x5[i],'ine//3xor_ine.txt')
            
            for i in range(int(self.block_size/2)):
                self.m.addConstr(self.x[index_start+i*2]==self.x[index_start+3*self.block_size+2*i])
                self.m.addConstr(self.x[index_start+i*2+1]==self.x[index_start+3*self.block_size+2*i+1])
        
        
        self.index=(self.index+int(self.block_size/2))%self.block_size
        self.m.addConstr(self.x[self.r*self.block_size*2+self.index*2]==self.value[0],name='output0')
        self.m.addConstr(self.x[self.r*self.block_size*2+self.index*2+1]==self.value[1],name='output1')


def get_value(simon,block_size):
    xnum=len(simon.x)
    anum=len(simon.aout)
    x=[0 for i in range(xnum)]
    aout=[0 for i in range(anum)]
    for v in simon.m.getVars():
        if v.varName[0]=='x':
            x[eval(v.varName[2:len(v.varName)-1])]=int(v.x)
        if v.varName[0]=='a':
            aout[eval(v.varName[5:len(v.varName)-1])]=int(v.x)
    s=''
    for r in range(simon.r+1):
        for i in range(block_size):
            if [x[(r*simon.block_size+i)*2],x[(r*simon.block_size+i)*2+1]]==[0,0]:
                s+='0'
            if [x[(r*simon.block_size+i)*2],x[(r*simon.block_size+i)*2+1]]==[0,1]:
                s+='1'
            if [x[(r*simon.block_size+i)*2],x[(r*simon.block_size+i)*2+1]]==[1,1]:
                s+='2'
        s+='\n'
    print(s)