from gurobipy import *
from cao_pylib import operation
from cao_pylib import sbox_model
import copy 
import sys

def print_table(x):
    for i in x:
        print(i)

def _state(x):
    res=[]
    for row in range(5):
        tmpx=[]
        for col in range(64):
            tmpx.append([x[(row*64+col)*2],x[(row*64+col)*2+1]])
        res.append(tmpx)
    return res


def print_state(x):
    for i in x:
        s=''
        for j in i:
            s+=str(j)
        print(s)

def _get_value(ascon):
    xnum=len(ascon.x)
    snum=len(ascon.sout)
    rounds=int(snum/640)
    x=[0 for i in range(xnum)]
    sout=[0 for i in range(snum)]

    for v in ascon.m.getVars():
        if v.varName[0]=='x':
            x[eval(v.varName[2:len(v.varName)-1])]=int(v.x)
        if v.varName[0]=='s':
            sout[eval(v.varName[5:len(v.varName)-1])]=int(v.x)
    
    xx=[]
    for r in range(int(xnum/640)):
        tmpx=[]
        for row in range(5):
            tmpxx=[]
            for col in range(64):
                tmp=[x[r*640+(row*64+col)*2],x[r*640+(row*64+col)*2+1]]
                if tmp==[0,0]:
                    tmpxx.append(0)
                if tmp==[0,1]:
                    tmpxx.append(1)
                if tmp==[1,1]:
                    tmpxx.append('?')
            tmpx.append(tmpxx)
        xx.append(tmpx)
    
    ss=[]
    for r in range(int(snum/640)):
        tmps=[]
        for row in range(5):
            tmpss=[]
            for col in range(64):
                tmp=[sout[r*640+(row*64+col)*2],sout[r*640+(row*64+col)*2+1]]
                if tmp==[0,0]:
                    tmpss.append(0)
                if tmp==[0,1]:
                    tmpss.append(1)
                if tmp==[1,1]:
                    tmpss.append('?')
            tmps.append(tmpss)
        ss.append(tmps)
    
    for r in range(rounds):
        print_state(xx[r])
        print()
        print_state(ss[r])
        print()
        print()
    if xnum>snum:
        print_state(xx[rounds])

def _get_value_inv(ascon):
    xnum=len(ascon.x)
    rounds=int(xnum/640)
    x=[0 for i in range(xnum)]

    for v in ascon.m.getVars():
        if v.varName[0]=='x':
            x[eval(v.varName[2:len(v.varName)-1])]=int(v.x)
    
    xx=[]
    for r in range(int(xnum/640)):
        tmpx=[]
        for row in range(5):
            tmpxx=[]
            for col in range(64):
                tmp=[x[r*640+(row*64+col)*2],x[r*640+(row*64+col)*2+1]]
                if tmp==[0,0]:
                    tmpxx.append(0)
                if tmp==[0,1]:
                    tmpxx.append(1)
                if tmp==[1,1]:
                    tmpxx.append(2)
            tmpx.append(tmpxx)
        xx.append(tmpx)   
    for i in range(len(xx)):
        print_state(xx[i])
        print()

x=['x'+str(i) for i in range(31)]
t=['t'+str(i) for i in range(9)]

def inv_xor_constr(m,x,y,t):
    remain_x=(len(x)-1)%3
    ine_num=int((len(x)-1)/3)

    tx=[]
    ty=[]
    for pos in range(4):
        tmpx=[x[pos]]
        if pos==0:
            for i in range(ine_num-1):
                tmpx.append(t[i])
        else:
            for i in range(1,ine_num):
                tmpx.append(x[i*3+pos])
        tx.append(tmpx)
    
    for i in range(ine_num-1):
        ty.append(t[i])
    if remain_x==0:
        ty.append(y)
    if remain_x!=0:
        ty.append(t[ine_num-1])
    
    for i in range(len(tx[0])):
        #print(tx[0][i]+tx[1][i]+tx[2][i]+tx[3][i]+ty[i])
        sbox_model.gen_model_from_ine(m,tx[0][i]+tx[1][i]+tx[2][i]+tx[3][i]+ty[i],'ine/_4xor_ine.txt')
    
    xx=[]
    yy=[]
    if remain_x==2:
        xx.append(t[ine_num-1])
        xx.append(x[len(x)-2])
        xx.append(x[len(x)-1])
        yy.append(y)
        sbox_model.gen_model_from_ine(m,xx[0]+xx[1]+xx[2]+yy[0],'ine/3xor_ine.txt')
        #print(xx[0]+xx[1]+xx[2]+yy[0])
    if remain_x==1:
        xx.append(t[ine_num-1])
        xx.append(x[len(x)-1])
        yy.append(y)
        sbox_model.gen_model_from_ine(m,xx[0]+xx[1]+yy[0],'ine/xor_ine.txt')
        #print(xx[0]+xx[1]+yy[0])

class Ascon():
    def __init__(self,round,block_size,index,value):
        self.r=round
        self.block_size=block_size
        #self.flag=flag  #flag=1,线性层多于非线性层，flag=0,线性层等于非线性层
        self.flag=1
        self.index=index
        if value==1:
            self.value=[0,1]
        if value==0:
            self.value=[0,0]
        self.rotate_right=[[0,19,28],[0,61,39],[0,1,6],[0,10,17],[0,7,41]]
        
        self.m=Model()
        if self.flag==0:
            self.x=self.m.addVars([i for i in  range(640*(self.r+1))],vtype=GRB.BINARY,name='x')
            self.sout=self.m.addVars([i for i in  range(640*(self.r))],vtype=GRB.BINARY,name='sout')
        if self.flag==1:
            self.x=self.m.addVars([i for i in  range(640*(self.r))],vtype=GRB.BINARY,name='x')
            self.sout=self.m.addVars([i for i in  range(640*self.r)],vtype=GRB.BINARY,name='sout')
        
        constr=LinExpr()
        for i in range(640):
            constr+=self.x[i]
        self.m.addConstr(constr>=1)

    def gen_constr(self):
        if self.flag==0:
            for r in range(self.r):
                indexStart=640*r
                xstate=[]
                sout_state=[]
                xorout_state=[]
                for row in range(5):
                    tmpx=[]
                    tmpxor=[]
                    tmps=[]
                    for col in range(64):
                        tmpx.append([self.x[(row*64+col)*2+indexStart],self.x[(row*64+col)*2+indexStart+1]])
                        tmps.append([self.sout[(row*64+col)*2+indexStart],self.sout[(row*64+col)*2+indexStart+1]])
                        tmpxor.append([self.x[(row*64+col)*2+indexStart+640],self.x[(row*64+col)*2+indexStart+641]])
                    xstate.append(tmpx)
                    xorout_state.append(tmpxor)
                    sout_state.append(tmps)

                #非线性层
                for col in range(64):
                    tmpx=[]
                    tmps=[]
                    for row in range(5):
                        tmpx+=xstate[row][col]
                        tmps.append(sout_state[row][col])
                    for i in range(5):
                        tmp=tmpx+tmps[i]
                        sbox_model.gen_model_from_ine(self.m,tmp,'ascon/ascon_ine_'+str(i)+'.txt')
                
                #线性层
                for row in range(5):
                    tmp1=operation.roate_right(sout_state[row],self.rotate_right[row][1])
                    tmp2=operation.roate_right(sout_state[row],self.rotate_right[row][2])
                    for col in range(64):
                        sbox_model.gen_model_from_ine(self.m,sout_state[row][col]+tmp1[col]+tmp2[col]+xorout_state[row][col],
                        'ine/3xor_ine.txt')
                
            index=640*self.r+self.index*2
            self.m.addConstr(self.x[index]==self.value[0],name='output0')
            self.m.addConstr(self.x[index+1]==self.value[1],name='output1')
        if self.flag==1:
            for r in range(self.r-1):
                indexStart=640*r
                xstate=[]
                sout_state=[]
                xorout_state=[]
                for row in range(5):
                    tmpx=[]
                    tmpxor=[]
                    tmps=[]
                    for col in range(64):
                        tmpx.append([self.x[(row*64+col)*2+indexStart],self.x[(row*64+col)*2+indexStart+1]])
                        tmps.append([self.sout[(row*64+col)*2+indexStart],self.sout[(row*64+col)*2+indexStart+1]])
                        tmpxor.append([self.x[(row*64+col)*2+indexStart+640],self.x[(row*64+col)*2+indexStart+641]])
                    xstate.append(tmpx)
                    xorout_state.append(tmpxor)
                    sout_state.append(tmps)

                #非线性层
                for col in range(64):
                    tmpx=[]
                    tmps=[]
                    for row in range(5):
                        tmpx+=xstate[row][col]
                        tmps.append(sout_state[row][col])
                    for i in range(5):
                        tmp=tmpx+tmps[i]
                        sbox_model.gen_model_from_ine(self.m,tmp,'ascon/ascon_ine_'+str(i)+'.txt')
                
                #线性层
                for row in range(5):
                    tmp1=operation.roate_right(sout_state[row],self.rotate_right[row][1])
                    tmp2=operation.roate_right(sout_state[row],self.rotate_right[row][2])
                    for col in range(64):
                        sbox_model.gen_model_from_ine(self.m,sout_state[row][col]+tmp1[col]+tmp2[col]+xorout_state[row][col],
                        'ine/3xor_ine.txt')
            
            #额外的非线性层
            xstate=[]
            sout_state=[]
            for row in range(5):
                tmpx=[]
                tmps=[]
                for col in range(64):
                    tmpx.append([self.x[(self.r-1)*640+(row*64+col)*2],self.x[(self.r-1)*640+(row*64+col)*2+1]])
                    tmps.append([self.sout[(self.r-1)*640+(row*64+col)*2],self.sout[(self.r-1)*640+(row*64+col)*2+1]])
                xstate.append(tmpx)
                sout_state.append(tmps)

            for col in range(64):
                tmpx=[]
                tmps=[]
                for row in range(5):
                    tmpx+=xstate[row][col]
                    tmps.append(sout_state[row][col])
                for i in range(5):
                    tmp=tmpx+tmps[i]
                    sbox_model.gen_model_from_ine(self.m,tmp,'ascon/ascon_ine_'+str(i)+'.txt')
            index=640*(self.r-1)+self.index*2
            self.m.addConstr(self.sout[index]==self.value[0],name='output0')
            self.m.addConstr(self.sout[index+1]==self.value[1],name='output1')


class Ascon_inv():
    def __init__(self,rounds,block_size,index,value):
        self.r=rounds
        self.block_size=block_size
        #self.flag=flag  
        self.flag=1          #flag=0,
        self.index=index
        if value==0:
            self.value=[0,0]
        if value==1:
            self.value=[0,1]
        self.rotate_right=[ [0, 3, 6, 9, 11, 12, 14, 15, 17, 18, 19, 21, 22, 24, 25, 27, 30, 33, 36, 38, 39, 41, 42, 44, 45, 47, 50, 53, 57, 60, 63],
                            [0, 1, 2, 3, 4, 8, 11, 13, 14, 16, 19, 21, 23, 24, 25, 27, 28, 29, 30, 35, 39, 43, 44, 45, 47, 48, 51, 53, 54, 55, 57, 60, 61],
                            [0, 2, 4, 6, 7, 10, 11, 13, 14, 15, 17, 18, 20, 23, 26, 27, 28, 32, 34, 35, 36, 37, 40, 42, 46, 47, 52, 58, 59, 60, 61, 62, 63],
                            [1, 2, 4, 6, 7, 9, 12, 17, 18, 21, 22, 23, 24, 26, 27, 28, 29, 31, 32, 33, 35, 36, 37, 40, 42, 44, 47, 48, 49, 53, 58, 61, 63],
                            [0, 1, 2, 3, 4, 5, 9, 10, 11, 13, 16, 20, 21, 22, 24, 25, 28, 29, 30, 31, 35, 36, 40, 41, 44, 45, 46, 47, 48, 50, 53, 55, 60, 61, 63]]
        self.rotate_num=[9,10,10,10,11]
        
        self.m=Model()
        if self.flag==0:
            self.x=self.m.addVars([i for i in range(640*(2*self.r))],vtype=GRB.BINARY,name='x')
            self.t=self.m.addVars([i for i in range(50*64*(self.r-1)*2)],vtype=GRB.BINARY,name='t')
            self.d=self.m.addVars([i for i in range(640)],vtype=GRB.BINARY,name='d')
        if self.flag==1:
            self.x=self.m.addVars([i for i in range(640*(2*self.r+1))],vtype=GRB.BINARY,name='x')
            self.t=self.m.addVars([i for i in range(50*64*(self.r)*2)],vtype=GRB.BINARY,name='t')
            self.d=self.m.addVars([i for i in range(320)],vtype=GRB.BINARY,name='d')
            #self.in=self.m.addVars([i for in range(640)],vtype=GRB.BINARY,name='in')

        constr=LinExpr()
        for i in range(640):
            constr+=self.x[i]
        self.m.addConstr(constr>=1)

        constr=LinExpr()
        for i in range(self.block_size):
            constr+=self.d[i]
        self.m.setObjective(constr,GRB.MAXIMIZE)

        for i in range(self.block_size):
            self.m.addConstr(self.x[i*2]>=self.d[i])
            self.m.addConstr(self.x[i*2+1]>=self.d[i])
            self.m.addConstr(self.x[i*2]+self.x[i*2+1]-self.d[i]<=1)
    
    def gen_constr(self):
        if self.flag==0:
            R=self.r-1
        if self.flag==1:
            R=self.r
        for r in range(R):
            indexStart=640*r*2
            xstate=[]
            xor_out_state=[]
            sout_state=[]
            for row in range(5):
                tmpx=[]
                tmpxor=[] #异或输出
                tmps=[]
                for col in range(64):
                    tmpx.append([self.x[(row*64+col)*2+indexStart],self.x[(row*64+col)*2+indexStart+1]])
                    tmpxor.append([self.x[(row*64+col)*2+indexStart+640],self.x[(row*64+col)*2+indexStart+641]])
                    tmps.append(([self.x[(row*64+col)*2+indexStart+1280],self.x[(row*64+col)*2+indexStart+1281]]))
                xstate.append(tmpx)
                sout_state.append(tmpxor)
                xor_out_state.append(tmps)
            self.m.update()
        
            #非线性层
            for col in range(64):
                tmpx=[]
                tmps=[]
                for row in range(5):
                    tmpx+=xstate[row][col]
                    tmps.append(sout_state[row][col])
                
                for i in range(5):
                    tmp=tmpx+tmps[i]
                    sbox_model.gen_model_from_ine(self.m,tmp,'ascon/ascon_inv_ine_'+str(i)+'.txt')
            
            #线性层
            tindexstart=50*64*r
            for row in range(5):
                tindex=tindexstart
                for i in range(row):
                    tindex+=self.rotate_num[i]*64
                tmpx=[]#线性层一行输入的值
                for col in range(64):
                    tmpx.append(sout_state[row][col])
                for col in range(64):
                    tmpt=[]
                    tmpxxx=[] 
                    for j in range(len(self.rotate_right[row])):
                        tmpxxx.append(tmpx[(col-self.rotate_right[row][j]+64)%64])
                    for j in range(self.rotate_num[row]):
                        tmpt.append([self.t[(tindex+self.rotate_num[row]*col+j)*2],self.t[(tindex+self.rotate_num[row]*col+j)*2+1]])

                    inv_xor_constr(self.m,tmpxxx,xor_out_state[row][col],tmpt)
        
        if self.flag==0:
            indexstart=640*(self.r-1)*2
            for col in range(64):
                tmpx=[]
                tmps=[]
                for row in range(5):
                    tmpx.append(self.x[indexstart+(row*64+col)*2])
                    tmpx.append(self.x[indexstart+(row*64+col)*2+1])
                    tmps.append([self.x[indexstart+(row*64+col)*2+640],self.x[indexstart+(row*64+col)*2+641]])
                for i in range(5):
                    tmp=tmpx+tmps[i]
                    sbox_model.gen_model_from_ine(self.m,tmp,'ascon/ascon_inv_ine_'+str(i)+'.txt')
                
        index=len(self.x)-640+self.index*2
        self.m.addConstr(self.x[index]==self.value[0],name='output0')
        self.m.addConstr(self.x[index+1]==self.value[1],name='output1')
          

'''
text=Ascon_inv(3,0,1,1)
text.gen_constr()
text.m.write('ascon\\ascon_inv.lp')
'''


if __name__=='__main__':
    t1=Ascon(3,1 , 0, 0)
    t1.gen_constr()
    t1.m.optimize()
    
