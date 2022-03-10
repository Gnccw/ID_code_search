from cao_pylib import sbox_model
cnf="(4'+5)(0+1')(1+2')(2+3')(3+5')(0'+5)"
sbox_model.gen_ine_from_cnf(cnf,'simon/br_ine.txt')
