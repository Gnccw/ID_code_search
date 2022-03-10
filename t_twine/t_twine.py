from gurobipy import *
from cao_pylib import sbox_model

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

def equal(m,x,y):
    m.addConstr(x[0]==y[0])
    m.addConstr(x[1]==y[1])

def t_twine(m,fore,fx,fxor_out,fk):
    rounds=fore[0]
    index=fore[1]
    value=fore[2]
    if value==0:
        value=[0,0]
    if value==1:
        value=[0,1]

    for r in range(rounds):
        x0=[]
        x1=[]
        tk=[]
        xor_out=[]

        for i in range(16):
            x0.append([fx[(r*16+i)*2],fx[(r*16+i)*2+1]])
            x1.append([fx[(r*16+i+16)*2],fx[(r*16+i+16)*2+1]])
        
        for i in range(6):
            tk.append([fk[(r*6+i)*2],fk[(r*6+i)*2+1]])
        
        for i in range(8):
            xor_out.append([fxor_out[(r*8+i)*2],fxor_out[(r*8+i)*2+1]])
        
        sbox_model.gen_model_from_ine(m, x0[0]+x0[1]+tk[0]+x1[0], 'ine//word_3xor_ine.txt')
        equal(m, x0[0],x1[5])
        
        sbox_model.gen_model_from_ine(m, x0[2]+x0[3]+tk[1]+x1[4], 'ine//word_3xor_ine.txt')
        equal(m, x0[2], x1[1])

        sbox_model.gen_model_from_ine(m, x0[4]+x0[5]+x1[12], 'ine//word_xor_ine.txt')
        equal(m, x0[4], x1[7])

        sbox_model.gen_model_from_ine(m, x0[6]+x0[7]+tk[2]+x1[8], 'ine//word_3xor_ine.txt')
        equal(m, x0[6], x1[3])

        sbox_model.gen_model_from_ine(m, x0[8]+x0[9]+tk[3]+x1[6], 'ine//word_3xor_ine.txt')
        equal(m, x0[8], x1[13])

        sbox_model.gen_model_from_ine(m, x0[10]+x0[11]+tk[4]+x1[2], 'ine//word_3xor_ine.txt')
        equal(m, x0[10], x1[9])

        sbox_model.gen_model_from_ine(m, x0[12]+x0[13]+x1[10], 'ine//word_xor_ine.txt')
        equal(m, x0[12], x1[15])

        sbox_model.gen_model_from_ine(m, x0[14]+x0[15]+tk[5]+x1[14], 'ine//word_3xor_ine.txt')
        equal(m, x0[14], x1[11])

    index=rounds*16+index
    m.addConstr(fx[index*2]==value[0])
    m.addConstr(fx[index*2+1]==value[1])


def t_twine_inv(m,last,bx,bxor_out,bk):
    rounds=last[0]
    index=last[1]
    value=last[2]
    if value==0:
        value=[0,0]
    if value==1:
        value=[0,1]
    
    for r in range(rounds):
        x0=[]
        x1=[]
        tk=[]
        xor_out=[]

        for i in range(16):
            x0.append([bx[(r*16+i)*2],bx[(r*16+i)*2+1]])
            x1.append([bx[(r*16+i+16)*2],bx[(r*16+i+16)*2+1]])
        
        for i in range(6):
            tk.append([bk[(r*6+i)*2],bk[(r*6+i)*2+1]])
        
        for i in range(8):
            xor_out.append([bxor_out[(r*8+i)*2],bxor_out[(r*8+i)*2+1]])
        

        sbox_model.gen_model_from_ine(m, x0[5]+x0[0]+tk[0]+x1[1], 'ine//word_3xor_ine.txt')
        equal(m, x1[0],x0[5])
        
        sbox_model.gen_model_from_ine(m, x0[1]+x0[4]+tk[1]+x1[3], 'ine//word_3xor_ine.txt')
        equal(m, x1[2], x0[1])

        sbox_model.gen_model_from_ine(m, x0[7]+x0[12]+x1[5], 'ine//word_xor_ine.txt')
        equal(m, x1[4], x0[7])

        sbox_model.gen_model_from_ine(m, x0[3]+x0[8]+tk[2]+x1[7], 'ine//word_3xor_ine.txt')
        equal(m, x1[6], x0[3])

        sbox_model.gen_model_from_ine(m, x0[13]+x0[6]+tk[3]+x1[9], 'ine//word_3xor_ine.txt')
        equal(m, x1[8], x0[13])

        sbox_model.gen_model_from_ine(m, x0[9]+x0[2]+tk[4]+x1[11], 'ine//word_3xor_ine.txt')
        equal(m, x1[10], x0[9])

        sbox_model.gen_model_from_ine(m, x0[15]+x0[10]+x1[13], 'ine//word_xor_ine.txt')
        equal(m, x1[12], x0[15])

        sbox_model.gen_model_from_ine(m, x0[11]+x0[14]+tk[5]+x1[15], 'ine//word_3xor_ine.txt')
        equal(m, x1[14], x0[11])

    index=rounds*16+index
    m.addConstr(bx[index*2]==value[0])
    m.addConstr(bx[index*2+1]==value[1])


def schedule(m,k):
    for r in range(34):
        k0=[]
        k1=[]
        for i in range(16):
            k0.append([k[(r*16+i)*2],k[(r*16+i)*2+1]])
            k1.append([k[(r*16+i)*2+32],k[(r*16+i)*2+1+32]])
        tmp=[k0[1],k0[0],k0[3],k0[4],k0[2],k0[5]]
        k0=tmp+k0[6:]

        for i in range(16):
            equal(m,k1[(i+16-6)%16] ,k0[i] )


def main(fore,last):
    frounds=fore[0]
    lrounds=last[0]
    m=Model()
    
    fx=m.addVars([i for i in range(32*frounds+32)],vtype=GRB.BINARY,name='fx')
    fxor_out=m.addVars([i for i in range(16*frounds)],vtype=GRB.BINARY,name='fxor_out')
    fk=m.addVars([i for i in range(12*frounds)],name='fk')
    
    bx=m.addVars([i for i in range(32*lrounds+32)],vtype=GRB.BINARY,name='bx')
    bxor_out=m.addVars([i for i in range(16*lrounds)],vtype=GRB.BINARY,name='bxor_out')
    bk=m.addVars([i for i in range(12*lrounds)],name='bk')

    tk=m.addVars([i for i in range(32*35)],vtype=GRB.BINARY,name='tk')

    #限制输入差分不为0-----------------------------
    constr=LinExpr()
    for i in range(32):
        constr+=fx[i]
    for i in range(32):
        constr+=tk[i]
    m.addConstr(constr>=1)
    #------------------------------------------------

    #将密钥算法和密码算法的模型结合起来-----------------
    for fr in range(fore[0]):
        tmpk=[]
        for i in range(6):
            tmpk.append([tk[(fr*16+(5-i))*2],tk[(fr*16+(5-i))*2+1]])
        for i in range(6):
            m.addConstr(fk[(fr*6+i)*2]==tmpk[i][0])
            m.addConstr(fk[(fr*6+i)*2+1]==tmpk[i][1])
    
    for br in range(last[0]):
        tmpk=[]
        for i in range(6):
            tmpk.append([tk[((last[0]+fore[0]-1-br)*16-i+5)*2],tk[((last[0]+fore[0]-1-br)*16+(5-i))*2+1]])
        for i in range(6):
            m.addConstr(bk[(br*6+i)*2]==tmpk[i][0])
            m.addConstr(bk[(br*6+i)*2+1]==tmpk[i][1])

    #-------------------------------------------------

    t_twine(m, fore, fx, fxor_out, fk)
    t_twine_inv(m, last, bx, bxor_out, bk)
    schedule(m, tk)

    m.Params.OutputFlag=0
    m.optimize()
    if m.status!=GRB.Status.INFEASIBLE:
        print(fore[1],'yes')
        fx_x=print_x(m, 'fx', len(fx), 16)
        bx_x=print_x(m, 'bx', len(bx), 16)
        k_x=print_x(m,'tk',len(tk),16)
        fk_x=print_x(m, 'fk', len(fk), 6)
        bk_x=print_x(m, 'bk', len(bk), 6)

        for i in range(fore[0]):
            print(fx_x[i],'-----',k_x[i])
        print(fx_x[fore[0]])
        
        print('-------------------------------------------')
        
        print(bx_x[last[0]])

        for i in range(last[0]):
            print(bx_x[last[0]-1-i],'-----',k_x[i+fore[0]])


for i in range(16):
    main([8,i,0],[10,i,1])
    