from os import linesep
from gurobipy import *
from cao_pylib import operation
from cao_pylib import sbox_model
from cao_pylib import cipher_sbox
import time
import copy

def var2bin(x):
    res=[]
    num=int(len(x)/2)
    #print(x)
    for i in range(num):
        if [x[i*2],x[i*2+1]]==[0,0]:
            res.append(0)
        if [x[i*2],x[i*2+1]]==[0,1]:
            res.append(1)
        if [x[i*2],x[i*2+1]]==[1,1]:
            res.append(2)
    return res


def bin2var(x):
    res=[]
    for i in range(len(x)):
        if x[i]=='0':
            res+=[0,0]
        if x[i]=='1':
            res+=[0,1]
        if x[i]=='2':
            res+=[1,1]
    return res

class Simon():
    def __init__(self,block_size,constra,rounds,forward_rounds,backward_rounds,data='',c='',filter='',meme=''):
        self.constra=constra #中间轮的矛盾
        self.rounds=rounds
        self.forward_rounds=forward_rounds
        self.backward_rounds=backward_rounds
        self.block_size=block_size
        self.half_block=int(self.block_size/2)
        if meme:
            self.meme=meme
        
        self.m=Model()
        self.ix=self.m.addVars([i for i in range(self.block_size*(self.rounds+2)*2)],vtype=GRB.BINARY,name='ix')
        self.iaout=self.m.addVars([i for i in range(self.rounds*self.block_size)],vtype=GRB.BINARY,name='iaout')

        self.fx=self.m.addVars([i for i in range(self.block_size*(self.forward_rounds+1)*2)],vtype=GRB.BINARY,name='fx')
        self.faout=self.m.addVars([i for i in range(self.forward_rounds*self.block_size)],vtype=GRB.BINARY,name='faout')
        self.fxout=self.m.addVars([i for i in range(self.forward_rounds*self.block_size)],vtype=GRB.BINARY,name='fxout')
        self.fk=self.m.addVars([i for i in range((self.forward_rounds-1)*int(self.block_size/2))],vtype=GRB.BINARY,name='fk')
        self.fp=self.m.addVars([i for i in range((self.forward_rounds)*int(self.block_size/2))],vtype=GRB.BINARY,name='fp')

        self.bx=self.m.addVars([i for i in range(self.block_size*(self.backward_rounds+1)*2)],vtype=GRB.BINARY,name='bx')
        self.baout=self.m.addVars([i for i in range(self.backward_rounds*self.block_size)],vtype=GRB.BINARY,name='baout')
        self.bxout=self.m.addVars([i for i in range(self.backward_rounds*self.block_size)],vtype=GRB.BINARY,name='bxout')
        self.bk=self.m.addVars([i for i in range((self.backward_rounds-1)*int(self.block_size/2))],vtype=GRB.BINARY,name='bk')
        self.bp=self.m.addVars([i for i in range((self.backward_rounds)*int(self.block_size/2))],vtype=GRB.BINARY,name='bp')

        #存储
        self.mem=self.m.addVar(vtype=GRB.INTEGER,name='mem')
        #明密文对数量
        self.data=self.m.addVar(vtype=GRB.INTEGER,name='data')
        #复杂度
        self.c=self.m.addVar(vtype=GRB.CONTINUOUS,name='c')
        #筛选率
        self.filter=self.m.addVar(vtype=GRB.INTEGER,lb=-30,ub=90,name='filter')
        self.m.setObjective(self.filter,GRB.MAXIMIZE)
        
        self.m.addConstr(self.filter>=filter)
        self.m.addConstr(self.c<=c)
        self.m.addConstr(self.data<=data)
        
        #明文自由度
        self.num2_in=self.m.addVars([i for i in range(self.block_size)],vtype=GRB.BINARY,name='num2_in')
        for i in range(self.block_size):
            self.m.addConstr(self.fx[(self.block_size*self.forward_rounds+i)*2]>=self.num2_in[i])
            self.m.addConstr(self.fx[(self.block_size*self.forward_rounds+i)*2+1]>=self.num2_in[i])
            self.m.addConstr(self.fx[(self.block_size*self.forward_rounds+i)*2]+self.fx[(self.block_size*self.forward_rounds+i)*2+1]-self.num2_in[i]<=1)
        #密文自由度
        self.num2_out=self.m.addVars([i for i in range(self.block_size)],vtype=GRB.BINARY,name='num2_out')
        for i in range(self.block_size):
            self.m.addConstr(self.bx[(self.block_size*self.backward_rounds+i)*2]>=self.num2_out[i])
            self.m.addConstr(self.bx[(self.block_size*self.backward_rounds+i)*2+1]>=self.num2_out[i])
            self.m.addConstr(self.bx[(self.block_size*self.backward_rounds+i)*2]+self.bx[(self.block_size*self.backward_rounds+i)*2+1]-self.num2_out[i]<=1)
        
        
        #可用明文对数量--------------------------------------------------
        data_min=LinExpr()
        data_max=LinExpr()

        for i in range(len(self.num2_in)):
            data_max+=(self.num2_in[i]+self.num2_out[i])
        data_max=data_max-1
        for i in range(len(self.num2_in)):
            data_min+=2*self.num2_in[i]+self.num2_out[i]
        data_min=data_min-1-self.block_size   #data_min:一个完整structure的有效明密文对数量

        self.m.addConstr(self.data>=data_min)
        self.m.addConstr(self.data<=data_max)
        self.plain=self.m.addVar(vtype=GRB.INTEGER,name='plain')
        self.m.addConstr(self.plain==(data_max-self.data))
        self.m.addConstr(self.plain==0)
        self.m.addConstr(self.meme>=self.data)
        #-----------------------------------------------------------------
        
        
        #限制输入输出差分不为0---------------------------------------------
        #num1_in,num1_out:输入差分为1的比特数
        self.num1_in=self.m.addVars([i for i in range(self.block_size)],vtype=GRB.BINARY,name='num1_in')
        for i in range(self.block_size):
            sbox_model.gen_model_from_ine(self.m,[self.ix[i*2],self.ix[i*2+1],self.num1_in[i]],'ine//num1.txt')

        self.num1_out=self.m.addVars([i for i in range(self.block_size)],vtype=GRB.BINARY,name='num1_out')
        for i in range(self.block_size):
            sbox_model.gen_model_from_ine(self.m,[self.ix[((self.rounds+1)*self.block_size+i)*2],self.ix[((self.rounds+1)*self.block_size+i)*2+1],self.num1_out[i]],'ine//num1.txt')
        
        constr_in=LinExpr()
        for i in range(self.block_size):
            constr_in+=self.num1_in[i]
        self.m.addConstr(constr_in>=1)
        constr_out=LinExpr()
        for i in range(self.block_size):
            constr_out+=self.num1_out[i]
        self.m.addConstr(constr_out>=1)
       #--------------------------------------------------------------------
        
        #前后向部分和不可能差分的连接-----------
        for i in range(self.block_size*2):
            self.m.addConstr(self.fx[i]==self.ix[i])
            self.m.addConstr(self.ix[(self.rounds+1)*self.block_size*2+i]==self.bx[i])
        #-------------------------------------
        
        #猜测密钥的个数-------------------------------------------------------
        self.sum_k=LinExpr()
        for i in range(len(self.fk)):
            self.sum_k+=self.fk[i]
        for i in range(len(self.bk)):
            self.sum_k+=self.bk[i]
        self.m.addConstr(self.sum_k<=900)       #保证所猜测的密钥个数不超过密钥总数
        self.m.setObjective(self.sum_k,GRB.MINIMIZE)

        self.bk0=LinExpr()
        self.fk0=LinExpr()
        for i in range(self.half_block):
            self.bk0+=self.bk[i]
        for i in range(self.half_block):
            self.fk0+=self.fk[i]
        
        #----------------------------------------------------------------------

        
        #存储空间---------------------------------------------------------
        self.sum_bp0=LinExpr()
        self.sum_fp0=LinExpr()
        for i in range(self.half_block):
            self.sum_bp0+=self.bp[i]
            self.sum_fp0+=self.fp[i]
        #self.m.addConstr(self.mem==(self.data-sum_bp0-sum_fp0))
        #----------------------------------------------------------------
        
        #总的概率--------------------------------------------------------
        self.sum_p=LinExpr()
        for i in range(len(self.fp)):
            self.sum_p+=self.fp[i]
        for i in range(len(self.bp)):
            self.sum_p+=self.bp[i]         
        #----------------------------------------------------------------
        
        
        #复杂度----------------------------------------------------------
        self.q=self.m.addVar(vtype=GRB.BINARY,name='q')
        c_0=LinExpr()
        c_1=LinExpr()
        c_0=self.data+self.sum_k-self.sum_p+self.sum_bp0
        c_1=self.data+self.sum_k-self.sum_p+self.sum_fp0
        self.m.addConstr(c_0+(1-self.q)*c_1==self.c)
        #----------------------------------------------------------------
    
        #筛选率
        #self.m.addConstr(self.filter==((self.data-self.sum_p)*1.44-self.bk0))
        self.m.addConstr(self.filter==(self.data-self.sum_p)-self.bk0)
        
    def gen_constr(self):
        self.m.update()
        
        #前向传播部分
        for r in range(self.forward_rounds):
            x0=[]
            x1=[]
            x0_L=[]
            x0_R=[]
            x1_L=[]
            x1_R=[]
            faout=[]
            fxout=[]
            fp=[]
            for i in range(self.block_size):
                x0.append([self.fx[((r)*self.block_size+i)*2],self.fx[((r)*self.block_size+i)*2+1]])
                x1.append([self.fx[((r+1)*self.block_size+i)*2],self.fx[((r+1)*self.block_size+i)*2+1]])
            for i in range(self.half_block):
                faout.append([self.faout[((r)*self.half_block+i)*2],self.faout[((r)*self.half_block+i)*2+1]])
                fxout.append([self.fxout[((r)*self.half_block+i)*2],self.fxout[((r)*self.half_block+i)*2+1]])
                
                fp.append(self.fp[r*self.half_block+i])
                x0_L.append(x0[i])
                x0_R.append(x0[i+self.half_block])
                x1_L.append(x1[i])
                x1_R.append(x1[i+self.half_block])

            x0_L1=operation.roate_left(x0_R,1)
            x0_L2=operation.roate_left(x0_R,2)
            x0_L8=operation.roate_left(x0_R,8)

            for i in range(self.half_block):
                #与操作
                sbox_model.gen_model_from_ine(self.m,x0_L1[i]+x0_L8[i]+faout[i],'ine//2and_ine.txt')
                #前向或操作
                sbox_model.gen_model_from_ine(self.m,x0_L2[i]+faout[i]+x1_R[i]+fxout[i],'ine//3xor_ine.txt')
                #前向或操作的筛选率
                sbox_model.gen_model_from_ine(self.m,fxout[i]+x0_L[i]+[fp[i]],'ine//p_ine.txt')
                #逆向或操作
                sbox_model.gen_model_from_ine(self.m,x0_L2[i]+faout[i]+x0_L[i]+x1_R[i],'ine//3xor_ine.txt')
                #轮输入输出左右交换
                self.m.addConstr(x0_R[i][0]==x1_L[i][0])
                self.m.addConstr(x0_R[i][1]==x1_L[i][1])
            
            
            #前后轮密钥的影响
            if r<self.forward_rounds-2:
                fk0=[]
                fk1=[]
                for i in range(self.half_block):
                    fk0.append(self.fk[r*self.half_block+i])
                    fk1.append(self.fk[(r+1)*self.half_block+i])
                for i in range(self.half_block):
                    self.m.addConstr(fk1[(i+self.half_block+1)%self.half_block]>=fk0[i])
                    self.m.addConstr(fk1[(i+self.half_block+8)%self.half_block]>=fk0[i])
                    self.m.addConstr(fk1[(i+self.half_block+2)%self.half_block]>=fk0[i])
            
            #筛选与操作影响的密钥 x1_R和fxout1
            if r<self.forward_rounds-1:
                fk=[]
                for i in range(self.half_block):
                    fk.append(self.fk[r*self.half_block+i])
                for i in range(self.half_block):
                    #self.m.addGenConstrMax(fk[(i+8)%self.half_block],[x1_L[(i+1)%self.half_block][0],x1_L[(i+1)%self.half_block][1]])
                    #self.m.addGenConstrMax(fk[(i+1)%self.half_block],[x1_L[(i+8)%self.half_block][0],x1_L[(i+8)%self.half_block][1]])
                    self.m.addConstr(2*fk[(i+8)%self.half_block]-x1_L[(i+1)%self.half_block][0]-x1_L[(i+1)%self.half_block][1]>=0)
                    self.m.addConstr(2*fk[(i+1)%self.half_block]-x1_L[(i+8)%self.half_block][0]-x1_L[(i+8)%self.half_block][1]>=0)
        
        
        #不可能差分部分
        for r in range(self.constra[0]):
            x0=[]
            x1=[]
            x0_L=[]
            x0_R=[]
            x1_L=[]
            x1_R=[]
            aout=[]
            for i in range(self.block_size):
                x0.append([self.ix[(r*self.block_size+i)*2],self.ix[(r*self.block_size+i)*2+1]])
                x1.append([self.ix[((r+1)*self.block_size+i)*2],self.ix[((r+1)*self.block_size+i)*2+1]])
            for i in range(self.half_block):
                x0_L.append(x0[i])
                x0_R.append(x0[i+self.half_block])
                x1_L.append(x1[i])
                x1_R.append(x1[i+self.half_block])
                aout.append([self.iaout[(self.half_block*r+i)*2],self.iaout[(self.half_block*r+i)*2+1]])
            
            x0_L1=operation.roate_left(x0_L,1)
            x0_L8=operation.roate_left(x0_L,8)
            x0_L2=operation.roate_left(x0_L,2)

            for i in range(self.half_block):
                sbox_model.gen_model_from_ine(self.m,x0_L1[i]+x0_L8[i]+aout[i],'ine//2and_ine.txt')
                sbox_model.gen_model_from_ine(self.m,x0_R[i]+aout[i]+x0_L2[i]+x1_L[i],'ine//3xor_ine.txt')
                self.m.addConstr(x0_L[i][0]==x1_R[i][0])
                self.m.addConstr(x0_L[i][1]==x1_R[i][1])
        
        for r in range(self.constra[0],self.rounds):
            x0=[]
            x1=[]
            x0_L=[]
            x0_R=[]
            x1_L=[]
            x1_R=[]
            aout=[]
            for i in range(self.block_size):
                x0.append([self.ix[((r+1)*self.block_size+i)*2],self.ix[((r+1)*self.block_size+i)*2+1]])
                x1.append([self.ix[((r+2)*self.block_size+i)*2],self.ix[((r+2)*self.block_size+i)*2+1]])
            for i in range(self.half_block):
                x0_L.append(x0[i])
                x0_R.append(x0[i+self.half_block])
                x1_L.append(x1[i])
                x1_R.append(x1[i+self.half_block])
                aout.append([self.iaout[(self.half_block*r+i)*2],self.iaout[(self.half_block*r+i)*2+1]])
            x0_L1=operation.roate_left(x0_L,1)
            x0_L8=operation.roate_left(x0_L,8)
            x0_L2=operation.roate_left(x0_L,2)

            for i in range(self.half_block):
                sbox_model.gen_model_from_ine(self.m,x0_L1[i]+x0_L8[i]+aout[i],'ine//2and_ine.txt')
                sbox_model.gen_model_from_ine(self.m,x0_L2[i]+x1_L[i]+aout[i]+x0_R[i],'ine//3xor_ine.txt')
                self.m.addConstr(x0_L[i][0]==x1_R[i][0])
                self.m.addConstr(x0_L[i][1]==x1_R[i][1])
            
        
        index=self.constra[0]*self.block_size+self.constra[1]
        index1=(self.constra[0]+1)*self.block_size+self.constra[1]
        if self.constra[2]==0:
            value=[0,0]
            value1=[0,1]
        if self.constra[2]==1:
            value=[0,1]
            value1=[0,0]
        self.m.addConstr(self.ix[index*2]==value[0])
        self.m.addConstr(self.ix[index*2+1]==value[1])
        self.m.addConstr(self.ix[index1*2]==value1[0])
        self.m.addConstr(self.ix[index1*2+1]==value1[1])
    
            
        #后向传播部分
        for r in range(self.backward_rounds):
            x0=[]
            x1=[]
            x0_L=[]
            x0_R=[]
            x1_L=[]
            x1_R=[]
            aout=[]
            xout=[]
            bp=[]
            
            for i in range(self.block_size):
                x0.append([self.bx[(r*self.block_size+i)*2],self.bx[(r*self.block_size+i)*2+1]])
                x1.append([self.bx[((r+1)*self.block_size+i)*2],self.bx[((r+1)*self.block_size+i)*2+1]])
            
            for i in range(self.half_block):
                aout.append([self.baout[(r*self.half_block+i)*2],self.baout[(r*self.half_block+i)*2+1]])
                xout.append([self.bxout[(r*self.half_block+i)*2],self.bxout[(r*self.half_block+i)*2+1]])
                bp.append(self.bp[r*self.half_block+i])

                x0_L.append(x0[i])
                x0_R.append(x0[i+self.half_block])
                x1_L.append(x1[i])
                x1_R.append(x1[i+self.half_block])

            x0_L1=operation.roate_left(x0_L,1)
            x0_L2=operation.roate_left(x0_L,2)
            x0_L8=operation.roate_left(x0_L,8)

            for i in range(self.half_block):
                #与操作
                sbox_model.gen_model_from_ine(self.m,x0_L1[i]+x0_L8[i]+aout[i],'ine//2and_ine.txt')
                #解密或操作
                sbox_model.gen_model_from_ine(self.m,aout[i]+x1_L[i]+x0_L2[i]+xout[i],'ine//3xor_ine.txt')
                #筛选率
                sbox_model.gen_model_from_ine(self.m,xout[i]+x0_R[i]+[bp[i]],'ine//p_ine.txt')
                #加密方向或操作
                sbox_model.gen_model_from_ine(self.m,x0_L2[i]+aout[i]+x0_R[i]+x1_L[i],'ine//3xor_ine.txt')
                #轮输入输出左右交换
                self.m.addConstr(x1_R[i][0]==x0_L[i][0])
                self.m.addConstr(x1_R[i][1]==x0_L[i][1])
            
            #前后轮密钥的影响
            if r<self.backward_rounds-2:
                bk0=[]
                bk1=[]
                for i in range(self.half_block):
                    bk0.append(self.bk[r*self.half_block+i])
                    bk1.append(self.bk[(r+1)*self.half_block+i])
                for i in range(self.half_block):
                    self.m.addConstr(bk1[(i+self.half_block+1)%self.half_block]>=bk0[i])
                    self.m.addConstr(bk1[(i+self.half_block+8)%self.half_block]>=bk0[i])
                    self.m.addConstr(bk1[(i+self.half_block+2)%self.half_block]>=bk0[i])

            #筛选与操作影响的密钥 x1_R和fxout1
            if r<self.backward_rounds-1:
                bk=[]
                for i in range(self.half_block):
                    bk.append(self.bk[r*self.half_block+i])
                for i in range(self.half_block):
                    #self.m.addGenConstrMax(fk[(i+8)%self.half_block],[x1_L[(i+1)%self.half_block][0],x1_L[(i+1)%self.half_block][1]])
                    #self.m.addGenConstrMax(fk[(i+1)%self.half_block],[x1_L[(i+8)%self.half_block][0],x1_L[(i+8)%self.half_block][1]])
                    self.m.addConstr(2*bk[(i+8)%self.half_block]-x0_L[(i+1)%self.half_block][0]-x0_L[(i+1)%self.half_block][1]>=0)
                    self.m.addConstr(2*bk[(i+1)%self.half_block]-x0_L[(i+8)%self.half_block][0]-x0_L[(i+8)%self.half_block][1]>=0)

            
            '''
            constr=LinExpr()
            constr1=LinExpr()
            for i in range(32,64):
                constr+=self.fx[i]
            for i in range(32):
                constr1+=self.fx[i]
            self.m.addConstr(constr==3)
            #self.m.addConstr(constr==5)
            self.m.addConstr(constr1==0)
            '''
            
def get_var_f(test):
    """
    依次打印不可能差分迹之前 扩展的轮的输入输出差分、筛选率、密钥
    """
    block_size=test.block_size
    half_block=int(block_size/2)

    x={}
    k={}
    p={}

    xout={}
    aout={}

    for v in test.m.getVars():
        if v.varName[:2]=='fx' and v.varName[2]!='o':
            x[eval(v.varName[3:len(v.varName)-1])]=int(v.x)
        if v.varName[:2]=='fk':
            k[eval(v.varName[3:len(v.varName)-1])]=int(v.x)
        if v.varName[:2]=='fp':
            p[eval(v.varName[3:len(v.varName)-1])]=int(v.x)
        if v.varName[:5]=='fxout':
            xout[eval(v.varName[6:len(v.varName)-1])]=int(v.x)
        if v.varName[:5]=='faout':
            aout[eval(v.varName[6:len(v.varName)-1])]=int(v.x)
    x=var2bin(x)
    xout=var2bin(xout)
    aout=var2bin(aout)

    for r in range(test.forward_rounds+1):
        tmp='\n'
        for i in range(block_size):
            tmp+=str(x[r*block_size+i])
        print(tmp)
    
    for r in range(test.forward_rounds):
        tmp='\n'
        for i in range(half_block):
            tmp+=str(xout[r*half_block+i])
        print(tmp)
    
    print('p:')
    for r in range(test.forward_rounds):
        tmp='\n'
        num=0
        for i in range(half_block):
            tmp+=str(p[r*half_block+i])
            if p[r*half_block+i]==1:
                num+=1
        tmp+='    '+str(num)
        print(tmp)
    
    print('k:')
    for r in range(test.forward_rounds-1):
        tmp='\n'
        num=0
        for i in range(half_block):
            tmp+=str(k[r*half_block+i])
            if k[r*half_block+i]==1:
                num+=1
        tmp+='    '+str(num)
        print(tmp)

def get_var_i(test):
    block_size=test.block_size
    x={}
    for v in test.m.getVars():
        if v.varName[:2]=='ix' and v.varName[2]!='o':
            x[eval(v.varName[3:len(v.varName)-1])]=int(v.x)
    x=var2bin(x)

    for r in range(test.rounds+2):
        tmp='\n'
        for i in range(block_size):
            tmp+=str(x[r*block_size+i])
        print(tmp)

def get_var_b(test):
    """
    依次打印不可能差分迹后 扩展的轮的输入输出差分、筛选率、密钥
    """
    block_size=test.block_size
    half_block=int(block_size/2)

    x={}
    k={}
    p={}

    xout={}
    aout={}

    for v in test.m.getVars():
        if v.varName[:2]=='bx' and v.varName[2]!='o':
            x[eval(v.varName[3:len(v.varName)-1])]=int(v.x)
        if v.varName[:2]=='bk':
            k[eval(v.varName[3:len(v.varName)-1])]=int(v.x)
        if v.varName[:2]=='bp':
            p[eval(v.varName[3:len(v.varName)-1])]=int(v.x)
        if v.varName[:5]=='bxout':
            xout[eval(v.varName[6:len(v.varName)-1])]=int(v.x)
        if v.varName[:5]=='baout':
            aout[eval(v.varName[6:len(v.varName)-1])]=int(v.x)
    print(len(x))
    x=var2bin(x)
    xout=var2bin(xout)
    aout=var2bin(aout)
    print(len(x))

    print('x')
    for r in range(test.backward_rounds+1):
        tmp='\n'
        for i in range(block_size):
            tmp+=str(x[r*block_size+i])
        print(tmp)
    
    print('xout')
    for r in range(test.backward_rounds):
        tmp='\n'
        for i in range(half_block):
            tmp+=str(xout[r*half_block+i])
        print(tmp)
    
    print('p')
    for r in range(test.backward_rounds):
        tmp='\n'
        num=0
        for i in range(half_block):
            tmp+=str(p[r*half_block+i])
            if p[r*half_block+i]==1:
                num+=1
        tmp+='    '+str(num)
        print(tmp)
    
    print('k')
    for r in range(test.backward_rounds-1):
        tmp='\n'
        num=0
        for i in range(half_block):
            tmp+=str(k[r*half_block+i])
            if k[r*half_block+i]==1:
                num+=1
        tmp+='    '+str(num)
        print(tmp)


def xechange_half(x):
    right_half=''
    left_half=''
    num=int(len(x)/2)
    for i in range(num):
        right_half+=x[i]
        left_half+=x[i+num]
    return left_half+right_half

#给2*block_size个输入var添加约束
def add_input_constr(cla,inputdiff):
    constr=LinExpr()
    const=1
    for i in range(len(inputdiff)):
        if inputdiff[i]==1:
            constr-=cla.bx[i]
            const-=1
        else:
            constr+=cla.bx[i]
    cla.m.addConstr(constr>=const)


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
        if i==0 or i=='0':
            res+=[0,0]
        if i==1 or i=='1':
            res+=[0,1]
        if i==2 or i=='2':
            res+=[1,1]
    return res


#将输入差分限制在某个集合中
def add_constr_in(cla,x):
    for i in range(cla.block_size):
        if x[i]=='2':
            cla.m.addConstr(cla.bx[i*2]+cla.bx[i*2+1]<=eval(x[i]))
        if x[i]=='1':
            cla.m.addConstr(cla.bx[i*2]==0)
            cla.m.addConstr(cla.bx[i*2+1]==1)
        if x[i]=='0':
            cla.m.addConstr(cla.bx[i*2]==0)
            cla.m.addConstr(cla.bx[i*2+1]==0)

#返回后向传播部分的输入差分的值
def get_bx0(cla):
    res=[0 for i in range(cla.block_size*2)]
    for v in cla.m.getVars():
        if v.varName[:2]=='bx' and v.varName[2]!='o':
            index=eval(v.varName[3:len(v.varName)-1])
            if index<cla.block_size*2:
                res[index]=int(v.x)
    return res

def get_fx0(cla):
    res=[0 for i in range(cla.block_size*2)]
    for v in cla.m.getVars():
        if v.varName[:2]=='fx' and v.varName[2]!='o':
            index=eval(v.varName[3:len(v.varName)-1])
            if index<cla.block_size*2:
                res[index]=int(v.x)
    return res

#给2*block_size个输入var赋值
def add_input_value(cla,inputdiff):
    inputdiff_var=value2var(inputdiff)
    for i in range(len(inputdiff_var)):
        cla.m.addConstr(cla.fx[i]==inputdiff_var[i])

def add_output_value(cla,outputdiff):
    outputdiff_var=value2var(outputdiff)
    for i in range(len(outputdiff_var)):
        cla.m.addConstr(cla.bx[i]==outputdiff_var[i])



'''

xf0='00000000000000000000000000000001'
#xb0='10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100012'
#xb0='0000000000000000000000000000000000000000000000000000000000000000 2000000000000000000000000000000000000000000000000000000000000020'


xb0=xechange_half(xb0)
xf0_var=value2var(xf0)
    

res=[]
bx_all_input=[]
for i in range(1):
    test=Simon(32,[10,48,0],13,4,4,data=999,c=256,filter=-1)
    test.gen_constr()
    add_constr_in(test,xb0)
    add_input_value(test,xf0_var)
    if len(bx_all_input)!=0:
        for j in bx_all_input:
                add_input_constr(test,j)
        test.m.Params.OutputFlag=0
        test.m.optimize()

        res.append([test.data.x,test.c.x,test.filter.x,test.plain])
        print([test.data.x,test.c.x,test.filter.x,test.plain])
        #get_var_i(test)
        #get_var_f(test)
        
        #print('---------------------------------------')
        #get_var_b(test)
        
        input_bx0=var2value(get_bx0(test))
        tmp=''
        for i in range(len(input_bx0)):
            tmp+=str(input_bx0[i])
        print(tmp)
        all_input_bx0=all_poss_bin_list(input_bx0)
        for k in all_input_bx0:
            bx_all_input.append(value2var(k))
    print(res)
    
'''

def fix_fx_bx(fx,bx,block_size,fx_round,bx_round,data='',c='',filter=''):
    test=Simon(block_size,[11,0,0],19,fx_round,bx_round,data,c,filter)
    test.gen_constr()
    add_input_value(test,fx)
    add_output_value(test,bx)
    #test.m.Params.OutputFlag=0
    test.m.optimize()
    if test.m.Status!=GRB.status.INFEASIBLE:

        input_bx0=var2value(get_bx0(test))
        input_fx0=var2value(get_fx0(test))
            
        tmp1=''
        for i in range(len(input_fx0)):
            tmp1+=str(input_fx0[i])
        print(tmp1)

        tmp=''
        for i in range(len(input_bx0)):
            tmp+=str(input_bx0[i])
        print(tmp)
        '''
        get_var_f(test)
        print('---------------------------------------')
        get_var_i(test)
        print('---------------------------------------')
        get_var_b(test)
        '''
        #print(test.mem)
        print(test.data)
        print(test.c)
        print(test.filter)
        print(test.bk)

'''
fx0="00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001"
#fx0=xechange_half(fx0)
bx0="10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001012"
#bx0="10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001222"
#bx0="10000000000000000000000000000000000000000000000000000000000000002000000000000000000000000000000000000000000000000000000001000022"
fix_fx_bx(fx0,bx0,128,5,5,data=999,c=999,filter=0)

'''
block_size=32
rounds=11
fp=open('simon_'+str(block_size)+'_r'+str(rounds)+'_constradict.txt')
constra=fp.read()
constra=constra.split('\n')
constra=constra[:len(constra)-1]

xf0='000000000000000000000000000000000000000000000000000000000000000'
xb0_all=[]
for i in constra:
    i=i.split(',')
    a=[eval(i[0]),eval(i[1]),eval(i[2])]
    print(a)
    if a[0]<0:
        continue

    for i in range(1):

        test=Simon(block_size,a,rounds,4,4,data=999,c=72,filter=-1,meme=49)
        test.gen_constr()

        #add_input_value(test,xf0)
        if len(xb0_all)!=0:
            for xb0 in xb0_all:
                add_input_constr(test,xb0)

        test.m.Params.OutputFlag=0
        test.m.optimize()
        
        if test.m.Status!=GRB.status.INFEASIBLE:

            input_bx0=var2value(get_bx0(test))
            input_fx0=var2value(get_fx0(test))
            
            tmp1=''
            for i in range(len(input_fx0)):
                tmp1+=str(input_fx0[i])
            print(tmp1)

            tmp=''
            for i in range(len(input_bx0)):
                tmp+=str(input_bx0[i])
            print(tmp)

            all_input_bx0=all_poss_bin_list(input_bx0)
            for k in all_input_bx0:
                xb0_all.append(value2var(k))
            
            
            get_var_f(test)
            print('-----------------------------------------------------------------------------------------')
            get_var_i(test)
            print('-----------------------------------------------------------------------------------------')
            get_var_b(test)
            
            
            #print(test.mem)
            print(test.data)
            print(test.filter)
            print(test.c)
            #print(test.bk0)
            '''
            for v in test.m.getVars():
                if v.varName[:3]=='fp[':
                    print(v)
                if v.VarName[:3]=='bp[':
                    print(v)
            break
            '''
        else:
            break
