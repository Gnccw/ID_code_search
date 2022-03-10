from cao_pylib import sbox_model
from gurobipy import *

'''
x=[x0,x1,x3,...,x63]
'''
p=[48, 1, 18, 35, 32, 49, 2, 19, 16, 33, 50, 3, 0, 17, 34, 51, 52, 5, 22, 39, 36, 53, 6, 23, 20, 37, 54, 7, 4, 21, 38, 55, 56, 9, 26, 43, 40, 57, 10, 27, 24, 41, 58, 11, 8, 25, 42, 59, 60, 13, 30, 47, 44, 61, 14, 31, 28, 45, 62, 15, 12, 29, 46, 63]
p_inv=[12, 1, 6, 11, 28, 17, 22, 27, 44, 33, 38, 43, 60, 49, 54, 59, 8, 13, 2, 7, 24, 29, 18, 23, 40, 45, 34, 39, 56, 61, 50, 55, 4, 9, 14, 3, 20, 25, 30, 19, 36, 41, 46, 35, 52, 57, 62, 51, 0, 5, 10, 15, 16, 21, 26, 31, 32, 37, 42, 47, 48, 53, 58, 63]

class Gift():
    def __init__(self,rounds,block_size,index,value):
        self.r=rounds
        self.block_size=block_size
        self.index=index
        if value==0:
            self.value=[0,0]
        if value==1:
            self.value=[0,1]
        
        self.m=Model()

        self.x=self.m.addVars([i for i in range(self.block_size*2*(self.r+1))],vtype=GRB.BINARY,name='x')
        self.sout=self.m.addVars([i for i in range(self.block_size*2*self.r)],vtype=GRB.BINARY,name='sout')

        constr=LinExpr()
        for i in range(self.block_size*2):
            constr+=self.x[i]
        self.m.addConstr(constr>=1)
    
    def gen_constr(self):
        for r in range(self.r):
            index=self.block_size*r

            x0=[]
            sout=[]
            x1=[]
            for i in range(64):
                x0.append([self.x[(index+i)*2],self.x[(index+i)*2+1]])
                sout.append([self.sout[(index+i)*2],self.sout[(index+i)*2+1]])
                x1.append([self.x[(index+i+self.block_size)*2],self.x[(index+i+self.block_size)*2+1]])
            
            for i in range(int(self.block_size/4)):
                tmpx=[]
                tmpy=[]
                for j in range(4):
                    tmpx+=x0[i*4+j]
                    tmpy+=sout[i*4+j]
                sbox_model.gen_model_from_ine(self.m,tmpx+tmpy,'gift/gift_sbox_ine_linear.txt')
            
            for i in range(64):
                self.m.addConstr(sout[i][0]==x1[p[i]][0])
                self.m.addConstr(sout[i][1]==x1[p[i]][1])
            
        index=self.r*self.block_size+self.index
        self.m.addConstr(self.x[index*2]==self.value[0],name='output0')
        self.m.addConstr(self.x[index*2+1]==self.value[1],name='output1')

class Gift_inv():
    def __init__(self,rounds,block_size,index,value):
        self.r=rounds
        self.block_size=block_size
        self.index=index
        if value==0:
            self.value=[0,0]
        if value==1:
            self.value=[0,1]
        
        self.m=Model()

        self.x=self.m.addVars([i for i in range(self.block_size*2*(self.r+1))],vtype=GRB.BINARY,name='x')
        self.sout=self.m.addVars([i for i in range(self.block_size*2*self.r)],vtype=GRB.BINARY,name='sout')

        constr=LinExpr()
        for i in range(self.block_size*2):
            constr+=self.x[i]
        self.m.addConstr(constr>=1)
    
    def gen_constr(self):
        for r in range(self.r):
            index=self.block_size*r

            x0=[]
            sout=[]
            x1=[]
            for i in range(64):
                x0.append([self.x[(index+i)*2],self.x[(index+i)*2+1]])
                sout.append([self.sout[(index+i)*2],self.sout[(index+i)*2+1]])
                x1.append([self.x[(index+i+self.block_size)*2],self.x[(index+i+self.block_size)*2+1]])
            
            for i in range(int(self.block_size/4)):
                tmpx=[]
                tmpy=[]
                for j in range(4):
                    tmpx+=x0[i*4+j]
                    tmpy+=sout[i*4+j]
                sbox_model.gen_model_from_ine(self.m,tmpx+tmpy,'gift/gift_sbox_ine_linear_inv.txt')
            
            for i in range(64):
                self.m.addConstr(sout[i][0]==x1[p_inv[i]][0])
                self.m.addConstr(sout[i][1]==x1[p_inv[i]][1])
            
        index=(self.r-1)*self.block_size+self.index
        self.m.addConstr(self.sout[index*2]==self.value[0],name='output0')
        self.m.addConstr(self.sout[index*2+1]==self.value[1],name='output1')


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


def _get_value(test):
    xnum=len(test.x)
    snum=len(test.sout)

    x=[0 for i in range(xnum)]
    s=[0 for i in range(snum)]

    for v in test.m.getVars():
        if v.varName[0]=='x':
            x[eval(v.varName[2:len(v.varName)-1])]=int(v.x)
        if v.varName[0]=='s':
            s[eval(v.varName[5:len(v.varName)-1])]=int(v.x)

    x=var2value(x)
    s=var2value(s)

    xstr=[]
    sstr=[]
    for r in range(test.r):
        tmpx=''
        tmps=''
        for i in range(test.block_size):
            if i%4==0:
                tmpx+=' '
                tmps+=' '
            tmpx+=str(x[r*test.block_size+i])
            tmps+=str(s[r*test.block_size+i])
        xstr.append(tmpx)
        sstr.append(tmps)
    tmpx=''
    for i in range(test.block_size):
        if i%4==0:
            tmpx+=' '
        tmpx+=str(x[test.r*test.block_size+i])
    xstr.append(tmpx)

    for r in range(test.r):
        print(xstr[r])
        print('-------sbox-------')
        print(sstr[r])
        print('-------player-----')
        print(xstr[r+1])
        print('\n')
        print()

        