from gurobipy import *
from cao_pylib import sbox_model
from cao_pylib import operation

def undis2diff(x):
    return x[0]+x[1]

def remove_value(m,x,value):
    v=[]
    for i in range(len(value)):
        if value[i]!=' ':
            v.append(eval(value[i]))
    constr=LinExpr()
    for i in range(len(v)):
        if v[i]==1:
            constr+=1-x[i]
        if v[i]==0:
            constr+=x[i]
    m.addConstr(constr>=1)

def left_right_xchange(x):
    """
    将x分为左右两部分并交换
    """
    half=int(len(x)/2)
    L=x[:half]
    R=x[half:]
    return R+L

def print_x(m,var_name,num,length):
    """
    var:变量名，字符串的形式
    num:变量的数量
    length:一般是block_size或half_block
    """
    '''以length比特为一行输出某个变量的值,该值存入列表中'''

    tmp=[0 for i in range(num)]
    for v in m.getVars():
        if v.varName[:len(var_name)]==var_name and v.varName[len(var_name)]=='[':
            tmp[eval(v.varName[len(var_name)+1:len(v.varName)-1])]=int(v.x)
    
    res=[]
    for r in range(int((num/length)/2)):
        s=''
        for i in range(length):
            if i%4==0:
                s+=' '
            s+=str(tmp[(r*length+i)*2]+tmp[(r*length+i)*2+1])
        res.append(s)
    return res

def Lblock(m,rounds,block_size,index,value,fx,fs_in,fs_out,fxor_out,fk):
    if value==0:
        value=[0,0]
    if value==1:
        value=[0,1]
    half_block=int(block_size/2)
    
    #生成约束
    for r in range(rounds):
        x0=[]
        x1=[]
        s_in=[]
        s_out=[]
        xor_out=[]
        k=[]
        x0_R=[]
        x0_L=[]
        x1_R=[]
        x1_L=[]
        x0_R_LR8=[]

        for i in range(block_size):
            x0.append([fx[(block_size*r+i)*2],fx[(block_size*r+i)*2+1]])
            x1.append([fx[(block_size*r+i+block_size)*2],fx[(block_size*r+i+block_size)*2+1]])
            
        for i in range(half_block):
            s_in.append([fs_in[(half_block*r+i)*2],fs_in[(half_block*r+i)*2+1]])
            s_out.append([fs_out[(half_block*r+i)*2],fs_out[(half_block*r+i)*2+1]])
            xor_out.append([fxor_out[(half_block*r+i)*2],fxor_out[(half_block*r+i)*2+1]])
            k.append([fk[(half_block*r+i)*2],fk[(half_block*r+i)*2+1]])

        x0_L=x0[:half_block]
        x0_R=x0[half_block:]
        x1_L=x1[:half_block]
        x1_R=x1[half_block:]
        x0_R_LR8=operation.roate_left(x0_R, 8)

        #加轮密钥
        for i in range(half_block):
            sbox_model.gen_model_from_ine(m, x0_L[i]+k[i]+s_in[i],'ine//xor_ine.txt')
            
        #s层
        for i in range(8):
            tmp=[]
            for j in range(4):
                tmp+=s_in[i*4+j]
            for j in range(4):
                tmp+=s_out[i*4+j]
            m.update()
            sbox_model.gen_model_from_ine(m,tmp, 'lblock//ine//sbox_ine'+str(i)+'.txt')
        #s层置换
        s_out_list=[]
        for i in range(8):
            tmp=[]
            for j in range(4):
                tmp.append(s_out[i*4+j])
            s_out_list.append(tmp)
        ss_out=[]
        ss_out+=s_out_list[1]
        ss_out+=s_out_list[3]
        ss_out+=s_out_list[0]
        ss_out+=s_out_list[2]
        ss_out+=s_out_list[5]
        ss_out+=s_out_list[7]
        ss_out+=s_out_list[4]
        ss_out+=s_out_list[6]
            
        #异或层
        for i in range(half_block):
            sbox_model.gen_model_from_ine(m, ss_out[i]+x0_R_LR8[i]+x1_L[i], 'ine//xor_ine.txt')
            
        #交换
        for i in range(half_block):
            m.addConstr(x0_L[i][0]==x1_R[i][0])
            m.addConstr(x0_L[i][1]==x1_R[i][1])
        
    #给最后一轮添加值
    index=(rounds*block_size+index)*2
    m.addConstr(fx[index]==value[0])
    m.addConstr(fx[index+1]==value[1])
            

def Lblock_inv(m,rounds,block_size,index,value,bx,bs_in,bs_out,bxor_out,bk):
    half_block=int(block_size/2)
    if value==0:
        value=[0,0]
    if value==1:
        value=[0,1]
    
    #def gen_constr(self):
    for r in range(rounds):
        x0=[]
        x1=[]
        s_in=[]
        s_out=[]
        xor_out=[]
        k=[]
        x0_R=[]
        x0_L=[]
        x1_R=[] 
        x1_L=[]

        for i in range(block_size):
            x0.append([bx[(block_size*r+i)*2],bx[(block_size*r+i)*2+1]])
            x1.append([bx[(block_size*r+i+block_size)*2],bx[(block_size*r+i+block_size)*2+1]])
            
        for i in range(half_block):
            s_in.append([bs_in[(half_block*r+i)*2],bs_in[(half_block*r+i)*2+1]])
            s_out.append([bs_out[(half_block*r+i)*2],bs_out[(half_block*r+i)*2+1]])
            xor_out.append([bxor_out[(half_block*r+i)*2],bxor_out[(half_block*r+i)*2+1]])
            k.append([bk[(half_block*r+i)*2],bk[(half_block*r+i)*2+1]])

            x0_L=x0[:half_block]
            x0_R=x0[half_block:]
            x1_L=x1[:half_block]
            x1_R=x1[half_block:]

        #加轮密钥
        for i in range(half_block):
            sbox_model.gen_model_from_ine(m, x0_L[i]+k[i]+s_in[i],'ine//xor_ine.txt')
            
        #s层
        for i in range(8):
            tmp=[]
            for j in range(4):
                tmp+=s_in[i*4+j]
            for j in range(4):
                tmp+=s_out[i*4+j]
            m.update()
            sbox_model.gen_model_from_ine(m,tmp, 'lblock//ine//sbox_ine'+str(i)+'.txt')
        #s层置换
        s_out_list=[]
        for i in range(8):
            tmp=[]
            for j in range(4):
                tmp.append(s_out[i*4+j])
            s_out_list.append(tmp)
        ss_out=[]
        ss_out+=s_out_list[1]
        ss_out+=s_out_list[3]
        ss_out+=s_out_list[0]
        ss_out+=s_out_list[2]
        ss_out+=s_out_list[5]
        ss_out+=s_out_list[7]
        ss_out+=s_out_list[4]
        ss_out+=s_out_list[6]
            
        #异或层
        for i in range(half_block):
            sbox_model.gen_model_from_ine(m, ss_out[i]+x0_R[i]+xor_out[i], 'ine//xor_ine.txt')
            
        #交换
        for i in range(half_block):
            m.addConstr(x0_L[i][0]==x1_R[i][0])
            m.addConstr(x0_L[i][1]==x1_R[i][1])
            
        xor_out_RL8=operation.roate_right(xor_out, 8)
        for i in range(half_block):
            m.addConstr(x1_L[i][0]==xor_out_RL8[i][0])
            m.addConstr(x1_L[i][1]==xor_out_RL8[i][1])
            
    #给最后一轮添加值
    x_last=[]
    for i in range(block_size*2):
        x_last.append(bx[block_size*rounds*2+i])
    x_last_c=left_right_xchange(x_last)
    m.addConstr(x_last_c[index*2]==value[0])
    m.addConstr(x_last_c[index*2+1]==value[1])


def schedule(m,kk):
    
    for r in range(31):
        k0=[]
        k1=[]
        for i in range(80):
            k0.append(kk[i+r*80])
            k1.append(kk[r*80+80+i])

        k0_RL29=operation.roate_left(k0, 29)
        sbox_model.gen_model_from_ine(m, [k0_RL29[0],k0_RL29[1],k0_RL29[2],k0_RL29[3],k1[0],k1[1],k1[2],k1[3]], 'lblock//ine//sbox_ine9.txt')
        sbox_model.gen_model_from_ine(m, [k0_RL29[4],k0_RL29[5],k0_RL29[6],k0_RL29[7],k1[4],k1[5],k1[6],k1[7]], 'lblock//ine//sbox_ine8.txt')
        
        for i in range(8,80):
            m.addConstr(k0_RL29[i]==k1[i])

def print_k(m,var_name,num,length):
    tmp=[i for i in range(num)]
    for v in m.getVars():
        if v.varName[:len(var_name)]==var_name and v.varName[len(var_name)]=='[':
            tmp[eval(v.varName[len(var_name)+1:len(v.varName)-1])]=int(v.x)
    
    res=[]
    for r in range(int(num/length)):
        s=''
        for i in range(length):
            if i%4==0:
                s+=' '
            s+=str(tmp[r*length+i])
        res.append(s)
    return res


def main(fore,last):
    m=Model()

#---建立变量--------------------------------------------------------------------------------------
    #密码算法的变量
    fx=m.addVars([i for i in range((fore[0]+1)*64*2)],vtype=GRB.BINARY,name='fx')
    fs_in=m.addVars([i for i in range(fore[0]*64)],vtype=GRB.BINARY,name='fs_in')
    fs_out=m.addVars([i for i in range(fore[0]*64)],vtype=GRB.BINARY,name='fs_out')
    fxor_out=m.addVars([i for i in range(fore[0]*64)],vtype=GRB.BINARY,name='fxor_out')
    #密钥的变量
    fk=m.addVars([i for i in range(fore[0]*64)],vtype=GRB.BINARY,name='fk')

    #密码算法的变量
    bx=m.addVars([i for i in range((last[0]+1)*64*2)],vtype=GRB.BINARY,name='bx')
    bs_in=m.addVars([i for i in range(last[0]*64)],vtype=GRB.BINARY,name='bs_in')
    bs_out=m.addVars([i for i in range(last[0]*64)],vtype=GRB.BINARY,name='bs_out')
    bxor_out=m.addVars([i for i in range(last[0]*64)],vtype=GRB.BINARY,name='bxor_out')
    #密钥的变量
    bk=m.addVars([i for i in range(last[0]*64)],vtype=GRB.BINARY,name='bk')
    
    kk=m.addVars([i for i in range(80*32)],vtype=GRB.BINARY,name='kk')
#---建立变量---------------------------------------------------------------------------------------
    
#------输入约束，差分不为0--------------------------------------------------------------------------    
    constr=LinExpr()
    for i in range(128):
        constr+=fx[i]
    for i in range(80):
        constr+=kk[i]
    m.addConstr(constr>=1)
    m.addConstr(constr<=40)
    
#------输入约束，差分不为0--------------------------------------------------------------------------

#-----最后一轮的变量--------------------------------------------------------------------------------
    constr=LinExpr()
    for i in range(16):
        for j in range(8):
            constr+=kk[i*80+j]
    m.addConstr(constr<=2)

    Lblock(m, fore[0],fore[1], fore[2], fore[3],fx,fs_in,fs_out,fxor_out,fk)
    Lblock_inv(m, last[0],last[1], last[2], last[3],bx,bs_in,bs_out,bxor_out,bk)
    schedule(m,kk)
    
    for fr in range(fore[0]):
        for i in range(32):
            m.addConstr(fk[(fr*32+i)*2]==0)
            m.addConstr(fk[(fr*32+i)*2+1]==kk[fr*80+i])
    
    for bf in range(last[0]):
        for i in range(32):
            m.addConstr(bk[(bf*32+i)*2]==0)
            m.addConstr(bk[(bf*32+i)*2+1]==kk[(fore[0]+last[0]-bf-1)*80+i])
            #m.addConstr(bk[(bf*32+i)*2+1]==0)
    value='0000 0010'
    value1='0000 0001'
    kk_v=[]
    for i in range(len(kk)):
        kk_v.append(kk[i])
    remove_value(m, kk_v[400:408],value)
    remove_value(m, kk_v[400:408],value1)
    
    m.Params.OutputFlag=0
    m.optimize()
    print(fore[2])
    if m.Status!=GRB.status.INFEASIBLE:
        print(fore[2],'yes')

        fx_s=print_x(m, 'fx', len(fx), 64)
        fsin_s=print_x(m,'fs_in',len(fs_in),32)
        fsout_s=print_x(m,'fs_out',len(fs_out),32)
    
        bx_s=print_x(m, 'bx', len(bx), 64)
        bsin_s=print_x(m,'bs_in',len(bs_in),32)
        bsout_s=print_x(m,'bs_out',len(bs_out),32)
        kk_x=print_k(m, 'kk', 80*32, 80)


        for r in range(fore[0]):
            print('| '+fx_s[r][:40]+' | '+fx_s[r][40:]+' | '+kk_x[r][0:]+'|') 
            #print(kk_x[r][:40])
            #print(fsin_s[r])
            #print('')
        print('| '+fx_s[fore[0]][:40]+' | '+fx_s[fore[0]][40:]+' | |')
        print('---------------------------------')
        print('| '+bx_s[last[0]][40:]+' | '+bx_s[last[0]][:40]+' | |')
        for r in range(last[0]):
            print('| '+bx_s[last[0]-r-1][40:]+' | '+bx_s[last[0]-r-1][:40]+' | '+kk_x[r+fore[0]][0:]+'|')
            #print(kk_x[r+fore[0]][:40])
            
            #print(bx_s[last[0]-r-1])
            #print()
        for i in range(len(kk_x)):
            print(kk_x[i])

for i in range(46,47):
    main([9,64,i,1],[8,64,i,0])
    
