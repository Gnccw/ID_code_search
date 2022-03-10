from cao_pylib import cipher_sbox
from cao_pylib import sbox
from cao_pylib import sbox_model
import sys
sys.path.append('c:\\Users\\Gnccw\\Desktop\\imposs difference')
import main

rectangle_sbox=cipher_sbox.rectangle_sbox

posspoint=main.undisturb_sbox(rectangle_sbox)
inv_posspoint=main.undisturb_sbox_invert(rectangle_sbox)

binposspoint=main.tri_pset2bin_pset(posspoint)
bin_inv_posspoint=main.tri_pset2bin_pset(inv_posspoint)

'''
sbox.write_patterns_csv(binposspoint,'rectangle\\rectangle_sbox.csv')
sbox.write_patterns_csv(bin_inv_posspoint,'rectangle\\rectangle_sbox_inv.csv')

cnf="(6'+12)(0'+14)(5'+12)(7'+14)(1'+12)(5'+14)(0'+1)(2'+3)(4'+5)(6'+7)(8+9')(9+10')(10+11')(11+15')(14'+15)(12'+13)(2'+14)(3+8'+13)(3'+11)(3+13'+15)(0+2+5+7+14')(1+2+5+6+12')(3'+12+13'+14')(11'+12+15)(2'+13)(13+14+15')(3+12+13'+14)(3'+12'+14+15')"
inv_cnf="(2'+10)(2'+12)(0'+12)(7'+10)(6'+7)(7'+12)(0'+1)(2'+3)(4'+5)(8+9')(9+14')(14+15')(10'+11)(11'+15)(12'+13)(1'+10)(1+3'+12)(5'+15)(5+8'+13)(5+11+13')(4'+10)(10+13+15')(3+10'+12)(1+3+4+7+12')(0+2+3'+4+7+10'+12')(5'+10'+12+13')(5'+10+11'+12')(11+12+15')(5+10+12+15')(4'+13)"

sbox_model.gen_ine_from_cnf(cnf,'rectangle\\rectangle_sbox_ine.txt')
sbox_model.gen_ine_from_cnf(inv_cnf,'rectangle\\inv_rectangle_sbox_ine.txt')
'''
main.test_ine('rectangle\\rectangle_sbox_ine.txt',binposspoint)



