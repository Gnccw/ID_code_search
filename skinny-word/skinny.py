from gurobipy import *
from cao_pylib import sbox_model
from cao_pylib import operation

def print_table(x):
    for xx in x:
        print(xx)

def list2value(x):
    if x==[0,0]:
        return 0
    if x==[0,1]:
        return 1
    if x==[1,1]:
        return 2
    else:
        return 'err'

def _get_value(skinny):
    xnum=(skinny.r*16+16)*2
    x=[0 for i in range(xnum)]
    sr=[0 for i in range(xnum-32)]

    for v in skinny.m.getVars():
        if v.varName[0]=='x':
            x[eval(v.varName[2:len(v.varName)-1])]=int(v.x)
        if v.varName[0]=='s':
            sr[eval(v.varName[5:len(v.varName)-1])]=int(v.x)
    
    xx=[]
    for r in range(skinny.r+1):
        tmpx=[[0 for i in range(4)] for j in range(4)]
        for row in range(4):
            for col in range(4):
                tmpx[row][col]=list2value([x[(r*16+row*4+col)*2],x[(r*16+row*4+col)*2+1]])
        xx.append(tmpx)

    ss=[]    
    for r in range(skinny.r):
        tmpx=[[0 for i in range(4)] for j in range(4)]
        for row in range(4):
            for col in range(4):
                tmpx[row][col]=list2value([sr[(r*16+row*4+col)*2],sr[(r*16+row*4+col)*2+1]])
        ss.append(tmpx)
    
    for r in range(skinny.r):
        print_table(xx[r])
        print('SR')
        print_table(ss[r])
        print('MC')
        print_table(xx[r+1])
        print('-----------------------------------------')


def _get_value_inv(skinny):
    xnum=(skinny.r*16+16)*2
    x=[0 for i in range(xnum)]
    sr=[0 for i in range(xnum-32)]

    for v in skinny.m.getVars():
        if v.varName[0]=='x':
            x[eval(v.varName[2:len(v.varName)-1])]=int(v.x)
        if v.varName[0]=='s':
            sr[eval(v.varName[5:len(v.varName)-1])]=int(v.x)
    
    xx=[]
    for r in range(skinny.r+1):
        tmpx=[[0 for i in range(4)] for j in range(4)]
        for row in range(4):
            for col in range(4):
                tmpx[row][col]=list2value([x[(r*16+row*4+col)*2],x[(r*16+row*4+col)*2+1]])
        xx.append(tmpx)

    ss=[]    
    for r in range(skinny.r):
        tmpx=[[0 for i in range(4)] for j in range(4)]
        for row in range(4):
            for col in range(4):
                tmpx[row][col]=list2value([sr[(r*16+row*4+col)*2],sr[(r*16+row*4+col)*2+1]])
        ss.append(tmpx)
    
    for r in range(skinny.r):
        print_table(xx[r])
        print('MC')
        print_table(ss[r])
        print('SR')
        print_table(xx[r+1])
        print('-----------------------------------------')


class Skinny():
    def __init__(self,rounds,block_size,index,value):
        self.r=rounds
        self.index=index
        self.block_size=16
        if value==0:
            self.value=[0,0]
        if value==1:
            self.value=[0,1]
        
        self.m=Model()
        self.x=self.m.addVars([i for i in range(16*(self.r+1)*2)],vtype=GRB.BINARY,name='x')
        self.sr_x=self.m.addVars([i for i in range(16*self.r*2)],vtype=GRB.BINARY,name='sr_x')

        constr=LinExpr()
        for i in range(32):
            constr+=self.x[i]
        self.m.addConstr(constr>=1)

    def gen_constr(self):
        for r in range(self.r):
            x=[[0 for i in range(4)] for j in range(4)]
            sr_x=[[0 for i in range(4)] for j in range(4)]
            y=[[0 for i in range(4)] for j in range(4)]
            for row in range(4):
                for col in range(4):
                    x[row][col]=[self.x[(r*16+row*4+col)*2],self.x[(r*16+row*4+col)*2+1]]
                    sr_x[row][col]=[self.sr_x[(r*16+row*4+col)*2],self.sr_x[(r*16+row*4+col)*2+1]]
                    y[row][col]=[self.x[((r+1)*16+row*4+col)*2],self.x[((r+1)*16+row*4+col)*2+1]]
            for row in range(4):
                x[row]=operation.roate_right(x[row],row)
            
            for row in range(4):
                for col in range(4):
                    self.m.addConstr(x[row][col][0]==sr_x[row][col][0])
                    self.m.addConstr(x[row][col][1]==sr_x[row][col][1])

            for col in range(4):
                sbox_model.gen_model_from_ine(self.m,sr_x[0][col]+sr_x[2][col]+sr_x[3][col]+y[0][col],'ine//word_3xor_ine.txt')
                self.m.addConstr(sr_x[0][col][0]==y[1][col][0])
                self.m.addConstr(sr_x[0][col][1]==y[1][col][1])
                sbox_model.gen_model_from_ine(self.m,sr_x[1][col]+sr_x[2][col]+y[2][col],'ine//word_xor_ine.txt')
                sbox_model.gen_model_from_ine(self.m,sr_x[0][col]+sr_x[2][col]+y[3][col],'ine//word_xor_ine.txt')
        index=(self.index+self.r*16)*2
        self.m.addConstr(self.x[index]==self.value[0])
        self.m.addConstr(self.x[index+1]==self.value[1])


class Skinny_inv():
    def __init__(self,rounds,block_size,index,value):
        self.r=rounds
        self.index=index
        self.block_size=16
        if value==0:
            self.value=[0,0]
        if value==1:
            self.value=[0,1]
        
        self.m=Model()
        self.x=self.m.addVars([i for i in range(16*(self.r+1)*2)],vtype=GRB.BINARY,name='x')
        self.sr_x=self.m.addVars([i for i in range(16*self.r*2)],vtype=GRB.BINARY,name='sr_x')

        constr=LinExpr()
        for i in range(32):
            constr+=self.x[i]
        self.m.addConstr(constr>=1)

    def gen_constr(self):
        for r in range(self.r):
            x=[[0 for i in range(4)] for j in range(4)]
            y=[[0 for i in range(4)] for j in range(4)]
            sr_x=[[0 for i in range(4)] for j in range(4)]

            for row in range(4):
                for col in range(4):
                    x[row][col]=[self.x[(r*16+row*4+col)*2],self.x[(r*16+row*4+col)*2+1]]
                    y[row][col]=[self.x[((r+1)*16+row*4+col)*2],self.x[((r+1)*16+row*4+col)*2+1]]
                    sr_x[row][col]=[self.sr_x[(r*16+row*4+col)*2],self.sr_x[(r*16+row*4+col)*2+1]]
            
            for col in range(4):
                sbox_model.gen_model_from_ine(self.m,x[1][col]+x[2][col]+x[3][col]+sr_x[1][col],'ine//word_3xor_ine.txt')
                self.m.addConstr(x[1][col][0]==sr_x[0][col][0])
                self.m.addConstr(x[1][col][1]==sr_x[0][col][1])
                sbox_model.gen_model_from_ine(self.m,x[1][col]+x[3][col]+sr_x[2][col],'ine//word_xor_ine.txt')
                sbox_model.gen_model_from_ine(self.m,x[0][col]+x[3][col]+sr_x[3][col],'ine//word_xor_ine.txt')
            
            for row in range(4):
                sr_x[row]=operation.roate_left(sr_x[row],row)

            for row in range(4):
                for col in range(4):
                    self.m.addConstr(y[row][col][0]==sr_x[row][col][0])
                    self.m.addConstr(y[row][col][1]==sr_x[row][col][1])
            
        index=(self.index+self.r*16)*2
        self.m.addConstr(self.x[index]==self.value[0])
        self.m.addConstr(self.x[index+1]==self.value[1])


for i in range(16):
    test=Skinny(6,16,i,0)
    test.gen_constr()
    test.m.Params.Outputflag=0
    test.m.write('test.lp')
    test.m.optimize()
    
    if test.m.status!=GRB.Status.INFEASIBLE:
        _get_value_inv(test)
        print(i)
        print('##################################################')
