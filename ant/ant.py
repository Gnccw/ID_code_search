from gurobipy import *
from cao_pylib import sbox_model

'''
x=[Xn,Xn-1,...,X0]
'''
def _perm128(var):
    p=[30, 3, 60, 57, 34, 7, 0, 61, 38, 11, 4, 1, 42, 15, 8, 5, 46, 19, 12, 9, 50, 23, 16, 13, 54, 27, 20, 17, 58, 31, 24, 21, 62, 35, 28, 25, 2, 39, 32, 29, 6, 43, 36, 33, 10, 47, 40, 37, 14, 51, 44, 41, 18, 55, 48, 45, 22, 59, 52, 49, 26, 63, 56, 53]
    res=[0 for i in range(64)]
    for i in range(64):
        res[i]=var[p[i]]
    return res


def _perm256(var):
    p=[62, 3, 124, 121, 66, 7, 0, 125, 70, 11, 4, 1, 74, 15, 8, 5, 78, 19, 12, 9, 82, 23, 16, 13, 86, 27, 20, 17, 90, 31, 24, 21, 94, 35, 28, 25, 98, 39, 32, 29, 102, 43, 36, 33, 106, 47, 40, 37, 110, 51, 44, 41, 114, 55, 48, 45, 118, 59, 52, 49, 122, 63, 56, 53, 126, 67, 60, 57, 2, 71, 64, 61, 6, 75, 68, 65, 10, 79, 72, 69, 14, 83, 76, 73, 18, 87, 80, 77, 22, 91, 84, 81, 26, 95, 88, 85, 30, 99, 92, 89, 34, 103, 96, 93, 38, 107, 100, 97, 42, 111, 104, 101, 46, 115, 108, 105, 50, 119, 112, 109, 54, 123, 116, 113, 58, 127, 120, 117]
    res=[0 for i in range(128)]
    for i in range(128):
        res[i]=var[p[i]]
    return res


def _G0(m,x,andout,xorout,G0_out):
    x_var=[]
    andout_var=[]
    xorout_var=[]
    G0_out_var=[]
    for i in range(int(len(x)/2)):
        x_var.append([x[i*2],x[i*2+1]])
    for i in range(int(len(andout)/2)):
        andout_var.append([andout[i*2],andout[i*2+1]])
    for i in range(int(len(xorout)/2)):
        xorout_var.append([xorout[i*2],xorout[i*2+1]])
    for i in range(int(len(G0_out)/2)):
        G0_out_var.append([G0_out[i*2],G0_out[i*2+1]])
    
    tmp_out=[]
    for i in range(len(xorout_var)):
        sbox_model.gen_model_from_ine(m,x_var[i*4+0]+x_var[i*4+1]+andout_var[i],'ine/xor_ine.txt')
        sbox_model.gen_model_from_ine(m,andout_var[i]+x_var[i*4+3]+xorout_var[i],'ine/2and_ine.txt')
        tmp_out+=[x_var[i*4+0],x_var[i*4+1],x_var[i*4+2],xorout_var[i]]
    
    tmp_in=[]
    if len(x_var)==64:
        tmp_in=_perm128(tmp_out)
    if len(x_var)==128:
        tmp_in=_perm256(tmp_out)
    
    for i in range(len(xorout_var)):
        sbox_model.gen_model_from_ine(m,tmp_in[i*4+0]+tmp_in[i*4+1]+andout_var[len(xorout_var)+i],'ine/2and_ine.txt')
        sbox_model.gen_model_from_ine(m,andout_var[len(xorout_var)+i]+tmp_in[i*4+3]+G0_out_var[i*4+3],'ine/xor_ine.txt')
        m.addConstr(tmp_in[i*4+0]==G0_out_var[i*4+0])
        m.addConstr(tmp_in[i*4+1]==G0_out_var[i*4+1])
        m.addConstr(tmp_in[i*4+2]==G0_out_var[i*4+2])


def _G1(m,x,andout,xorout,G0_out):
    x_var=[]
    andout_var=[]
    xorout_var=[]
    G0_out_var=[]
    for i in range(int(len(x)/2)):
        x_var.append([x[i*2],x[i*2+1]])
    for i in range(int(len(andout)/2)):
        andout_var.append([andout[i*2],andout[i*2+1]])
    for i in range(int(len(xorout)/2)):
        xorout_var.append([xorout[i*2],xorout[i*2+1]])
    for i in range(int(len(G0_out)/2)):
        G0_out_var.append([G0_out[i*2],G0_out[i*2+1]])
    
    tmp_out=[]
    for i in range(len(xorout_var)):
        sbox_model.gen_model_from_ine(m,x_var[i*4+0]+x_var[i*4+1]+andout_var[i],'ine/xor_ine.txt')
        sbox_model.gen_model_from_ine(m,andout_var[i]+x_var[i*4+2]+xorout_var[i],'ine/2and_ine.txt')
        tmp_out+=[x_var[i*4+0],x_var[i*4+1],xorout_var[i],x_var[i*4+3]]
    
    tmp_in=[]
    if len(x_var)==64:
        tmp_in=_perm128(tmp_out)
    if len(x_var)==128:
        tmp_in=_perm256(tmp_out)
    
    for i in range(len(xorout_var)):
        sbox_model.gen_model_from_ine(m,tmp_in[i*4+0]+tmp_in[i*4+1]+andout_var[len(xorout_var)+i],'ine/2and_ine.txt')
        sbox_model.gen_model_from_ine(m,andout_var[len(xorout_var)+i]+tmp_in[i*4+2]+G0_out_var[i*4+2],'ine/xor_ine.txt')
        m.addConstr(tmp_in[i*4+0]==G0_out_var[i*4+0])
        m.addConstr(tmp_in[i*4+1]==G0_out_var[i*4+1])
        m.addConstr(tmp_in[i*4+3]==G0_out_var[i*4+3])


class Ant():
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
        self.and_out_G0=self.m.addVars([i for i in range(self.r*int(self.block_size/2))],vtype=GRB.BINARY,name='and_out_G0')
        self.xor_out_G0=self.m.addVars([i for i in range(self.r*int(self.block_size/4))],vtype=GRB.BINARY,name='xor_out_G0')
        self.and_out_G1=self.m.addVars([i for i in range(self.r*int(self.block_size/2))],vtype=GRB.BINARY,name='and_out_G1')
        self.xor_out_G1=self.m.addVars([i for i in range(self.r*int(self.block_size/4))],vtype=GRB.BINARY,name='xor_out_G1')
        self.out_G0=self.m.addVars([i for i in range(self.block_size*self.r)],vtype=GRB.BINARY,name='out_G0')
        self.out_G1=self.m.addVars([i for i in range(self.block_size*self.r)],vtype=GRB.BINARY,name='out_G1')

        self.m.update()
        self.x_var=[]
        for i in range((self.r+1)*self.block_size):
            self.x_var.append([self.x[i*2],self.x[i*2+1]])
        
        
        self.and_out_G0_var=[]
        for i in range(self.r*int(self.block_size/4)):
            self.and_out_G0_var.append([self.and_out_G0[i*2],self.and_out_G0[i*2+1]])
        

        self.and_out_G1_var=[]
        for i in range(self.r*int(self.block_size/4)):
            self.and_out_G1_var.append([self.and_out_G1[i*2],self.and_out_G1[i*2+1]])

        self.xor_out_G0_var=[]
        for i in range(self.r*int(self.block_size/8)):
            self.xor_out_G0_var.append([self.xor_out_G0[i*2],self.xor_out_G0[i*2+1]])
        
        
        self.xor_out_G1_var=[]
        for i in range(self.r*int(self.block_size/8)):
            self.xor_out_G1_var.append([self.xor_out_G1[i*2],self.xor_out_G1[i*2+1]])
        
        self.out_G0_var=[]
        for i in range(self.r*int(self.block_size/2)):
            self.out_G0_var.append([self.out_G0[i*2],self.out_G0[i*2+1]])
        
        self.out_G1_var=[]
        for i in range(self.r*int(self.block_size/2)):
            self.out_G1_var.append([self.out_G1[i*2],self.out_G1[i*2+1]])

        constr=LinExpr()
        for i in range(self.block_size*2):
            constr+=self.x[i]
        self.m.addConstr(constr>=1)


    def gen_constr(self):
        for r in range(self.r):
            left_state=[]
            for i in range(int(self.block_size/2)):
                left_state.append(self.x[r*self.block_size*2+self.block_size-i])


Ant(2,128,2,1)