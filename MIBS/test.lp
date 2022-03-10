\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
 
Subject To
 R0: x[0] + x[1] + x[2] + x[3] + x[4] + x[5] + x[6] + x[7] + x[8] + x[9]
   + x[10] + x[11] + x[12] + x[13] + x[14] + x[15] + x[16] + x[17] + x[18]
   + x[19] + x[20] + x[21] + x[22] + x[23] + x[24] + x[25] + x[26] + x[27]
   + x[28] + x[29] + x[30] + x[31] + x[32] + x[33] + x[34] + x[35] + x[36]
   + x[37] + x[38] + x[39] + x[40] + x[41] + x[42] + x[43] + x[44] + x[45]
   + x[46] + x[47] + x[48] + x[49] + x[50] + x[51] + x[52] + x[53] + x[54]
   + x[55] + x[56] + x[57] + x[58] + x[59] + x[60] + x[61] + x[62] + x[63]
   + x[64] + x[65] + x[66] + x[67] + x[68] + x[69] + x[70] + x[71] + x[72]
   + x[73] + x[74] + x[75] + x[76] + x[77] + x[78] + x[79] + x[80] + x[81]
   + x[82] + x[83] + x[84] + x[85] + x[86] + x[87] + x[88] + x[89] + x[90]
   + x[91] + x[92] + x[93] + x[94] + x[95] + x[96] + x[97] + x[98] + x[99]
   + x[100] + x[101] + x[102] + x[103] + x[104] + x[105] + x[106] + x[107]
   + x[108] + x[109] + x[110] + x[111] + x[112] + x[113] + x[114] + x[115]
   + x[116] + x[117] + x[118] + x[119] + x[120] + x[121] + x[122] + x[123]
   + x[124] + x[125] + x[126] + x[127] >= 1
 R1: - x[0] + x[1] >= 0
 R2: - x[2] + x[3] >= 0
 R3: - x[4] + x[5] >= 0
 R4: - x[6] + x[7] >= 0
 R5: - x[7] + sout[7] >= 0
 R6: sout[0] - sout[1] >= 0
 R7: sout[1] - sout[2] >= 0
 R8: sout[2] - sout[3] >= 0
 R9: sout[3] - sout[4] >= 0
 R10: sout[4] - sout[5] >= 0
 R11: sout[5] - sout[6] >= 0
 R12: sout[6] - sout[7] >= 0
 R13: - x[5] + sout[7] >= 0
 R14: - x[3] + sout[7] >= 0
 R15: - x[1] + sout[7] >= 0
 R16: x[1] + x[3] + x[5] + x[7] - sout[0] >= 0
 R17: - x[8] + x[9] >= 0
 R18: - x[10] + x[11] >= 0
 R19: - x[12] + x[13] >= 0
 R20: - x[14] + x[15] >= 0
 R21: - x[15] + sout[15] >= 0
 R22: sout[8] - sout[9] >= 0
 R23: sout[9] - sout[10] >= 0
 R24: sout[10] - sout[11] >= 0
 R25: sout[11] - sout[12] >= 0
 R26: sout[12] - sout[13] >= 0
 R27: sout[13] - sout[14] >= 0
 R28: sout[14] - sout[15] >= 0
 R29: - x[13] + sout[15] >= 0
 R30: - x[11] + sout[15] >= 0
 R31: - x[9] + sout[15] >= 0
 R32: x[9] + x[11] + x[13] + x[15] - sout[8] >= 0
 R33: - x[16] + x[17] >= 0
 R34: - x[18] + x[19] >= 0
 R35: - x[20] + x[21] >= 0
 R36: - x[22] + x[23] >= 0
 R37: - x[23] + sout[23] >= 0
 R38: sout[16] - sout[17] >= 0
 R39: sout[17] - sout[18] >= 0
 R40: sout[18] - sout[19] >= 0
 R41: sout[19] - sout[20] >= 0
 R42: sout[20] - sout[21] >= 0
 R43: sout[21] - sout[22] >= 0
 R44: sout[22] - sout[23] >= 0
 R45: - x[21] + sout[23] >= 0
 R46: - x[19] + sout[23] >= 0
 R47: - x[17] + sout[23] >= 0
 R48: x[17] + x[19] + x[21] + x[23] - sout[16] >= 0
 R49: - x[24] + x[25] >= 0
 R50: - x[26] + x[27] >= 0
 R51: - x[28] + x[29] >= 0
 R52: - x[30] + x[31] >= 0
 R53: - x[31] + sout[31] >= 0
 R54: sout[24] - sout[25] >= 0
 R55: sout[25] - sout[26] >= 0
 R56: sout[26] - sout[27] >= 0
 R57: sout[27] - sout[28] >= 0
 R58: sout[28] - sout[29] >= 0
 R59: sout[29] - sout[30] >= 0
 R60: sout[30] - sout[31] >= 0
 R61: - x[29] + sout[31] >= 0
 R62: - x[27] + sout[31] >= 0
 R63: - x[25] + sout[31] >= 0
 R64: x[25] + x[27] + x[29] + x[31] - sout[24] >= 0
 R65: - x[32] + x[33] >= 0
 R66: - x[34] + x[35] >= 0
 R67: - x[36] + x[37] >= 0
 R68: - x[38] + x[39] >= 0
 R69: - x[39] + sout[39] >= 0
 R70: sout[32] - sout[33] >= 0
 R71: sout[33] - sout[34] >= 0
 R72: sout[34] - sout[35] >= 0
 R73: sout[35] - sout[36] >= 0
 R74: sout[36] - sout[37] >= 0
 R75: sout[37] - sout[38] >= 0
 R76: sout[38] - sout[39] >= 0
 R77: - x[37] + sout[39] >= 0
 R78: - x[35] + sout[39] >= 0
 R79: - x[33] + sout[39] >= 0
 R80: x[33] + x[35] + x[37] + x[39] - sout[32] >= 0
 R81: - x[40] + x[41] >= 0
 R82: - x[42] + x[43] >= 0
 R83: - x[44] + x[45] >= 0
 R84: - x[46] + x[47] >= 0
 R85: - x[47] + sout[47] >= 0
 R86: sout[40] - sout[41] >= 0
 R87: sout[41] - sout[42] >= 0
 R88: sout[42] - sout[43] >= 0
 R89: sout[43] - sout[44] >= 0
 R90: sout[44] - sout[45] >= 0
 R91: sout[45] - sout[46] >= 0
 R92: sout[46] - sout[47] >= 0
 R93: - x[45] + sout[47] >= 0
 R94: - x[43] + sout[47] >= 0
 R95: - x[41] + sout[47] >= 0
 R96: x[41] + x[43] + x[45] + x[47] - sout[40] >= 0
 R97: - x[48] + x[49] >= 0
 R98: - x[50] + x[51] >= 0
 R99: - x[52] + x[53] >= 0
 R100: - x[54] + x[55] >= 0
 R101: - x[55] + sout[55] >= 0
 R102: sout[48] - sout[49] >= 0
 R103: sout[49] - sout[50] >= 0
 R104: sout[50] - sout[51] >= 0
 R105: sout[51] - sout[52] >= 0
 R106: sout[52] - sout[53] >= 0
 R107: sout[53] - sout[54] >= 0
 R108: sout[54] - sout[55] >= 0
 R109: - x[53] + sout[55] >= 0
 R110: - x[51] + sout[55] >= 0
 R111: - x[49] + sout[55] >= 0
 R112: x[49] + x[51] + x[53] + x[55] - sout[48] >= 0
 R113: - x[56] + x[57] >= 0
 R114: - x[58] + x[59] >= 0
 R115: - x[60] + x[61] >= 0
 R116: - x[62] + x[63] >= 0
 R117: - x[63] + sout[63] >= 0
 R118: sout[56] - sout[57] >= 0
 R119: sout[57] - sout[58] >= 0
 R120: sout[58] - sout[59] >= 0
 R121: sout[59] - sout[60] >= 0
 R122: sout[60] - sout[61] >= 0
 R123: sout[61] - sout[62] >= 0
 R124: sout[62] - sout[63] >= 0
 R125: - x[61] + sout[63] >= 0
 R126: - x[59] + sout[63] >= 0
 R127: - x[57] + sout[63] >= 0
 R128: x[57] + x[59] + x[61] + x[63] - sout[56] >= 0
 R129: sout[9] + sout[17] + sout[25] + sout[33] + sout[41] - sout[49]
   + xorout[1] >= 0
 R130: sout[9] + sout[17] + sout[25] + sout[33] - sout[41] + sout[49]
   + xorout[1] >= 0
 R131: sout[9] + sout[17] + sout[25] - sout[33] + sout[41] + sout[49]
   + xorout[1] >= 0
 R132: sout[9] + sout[17] - sout[25] + sout[33] + sout[41] + sout[49]
   + xorout[1] >= 0
 R133: sout[9] - sout[17] + sout[25] + sout[33] + sout[41] + sout[49]
   + xorout[1] >= 0
 R134: - sout[9] + sout[17] + sout[25] + sout[33] + sout[41] + sout[49]
   + xorout[1] >= 0
 R135: sout[9] + sout[17] + sout[25] + sout[33] + sout[41] + sout[49]
   - xorout[1] >= 0
 R136: sout[9] + sout[17] + sout[25] - sout[33] - sout[41] - sout[49]
   + xorout[1] >= -2
 R137: sout[9] + sout[17] - sout[25] + sout[33] - sout[41] - sout[49]
   + xorout[1] >= -2
 R138: sout[9] - sout[17] + sout[25] + sout[33] - sout[41] - sout[49]
   + xorout[1] >= -2
 R139: - sout[9] + sout[17] + sout[25] + sout[33] - sout[41] - sout[49]
   + xorout[1] >= -2
 R140: sout[9] + sout[17] - sout[25] - sout[33] + sout[41] - sout[49]
   + xorout[1] >= -2
 R141: sout[9] - sout[17] + sout[25] - sout[33] + sout[41] - sout[49]
   + xorout[1] >= -2
 R142: - sout[9] + sout[17] + sout[25] - sout[33] + sout[41] - sout[49]
   + xorout[1] >= -2
 R143: sout[9] - sout[17] - sout[25] + sout[33] + sout[41] - sout[49]
   + xorout[1] >= -2
 R144: - sout[9] + sout[17] - sout[25] + sout[33] + sout[41] - sout[49]
   + xorout[1] >= -2
 R145: - sout[9] - sout[17] + sout[25] + sout[33] + sout[41] - sout[49]
   + xorout[1] >= -2
 R146: sout[9] + sout[17] - sout[25] - sout[33] - sout[41] + sout[49]
   + xorout[1] >= -2
 R147: sout[9] - sout[17] + sout[25] - sout[33] - sout[41] + sout[49]
   + xorout[1] >= -2
 R148: - sout[9] + sout[17] + sout[25] - sout[33] - sout[41] + sout[49]
   + xorout[1] >= -2
 R149: sout[9] - sout[17] - sout[25] + sout[33] - sout[41] + sout[49]
   + xorout[1] >= -2
 R150: - sout[9] + sout[17] - sout[25] + sout[33] - sout[41] + sout[49]
   + xorout[1] >= -2
 R151: - sout[9] - sout[17] + sout[25] + sout[33] - sout[41] + sout[49]
   + xorout[1] >= -2
 R152: sout[9] - sout[17] - sout[25] - sout[33] + sout[41] + sout[49]
   + xorout[1] >= -2
 R153: - sout[9] + sout[17] - sout[25] - sout[33] + sout[41] + sout[49]
   + xorout[1] >= -2
 R154: - sout[9] - sout[17] + sout[25] - sout[33] + sout[41] + sout[49]
   + xorout[1] >= -2
 R155: - sout[9] - sout[17] - sout[25] + sout[33] + sout[41] + sout[49]
   + xorout[1] >= -2
 R156: sout[9] - sout[17] - sout[25] - sout[33] - sout[41] - sout[49]
   + xorout[1] >= -4
 R157: - sout[9] + sout[17] - sout[25] - sout[33] - sout[41] - sout[49]
   + xorout[1] >= -4
 R158: - sout[9] - sout[17] + sout[25] - sout[33] - sout[41] - sout[49]
   + xorout[1] >= -4
 R159: - sout[9] - sout[17] - sout[25] + sout[33] - sout[41] - sout[49]
   + xorout[1] >= -4
 R160: - sout[9] - sout[17] - sout[25] - sout[33] + sout[41] - sout[49]
   + xorout[1] >= -4
 R161: - sout[9] - sout[17] - sout[25] - sout[33] - sout[41] + sout[49]
   + xorout[1] >= -4
 R162: - sout[48] + xorout[0] >= 0
 R163: - sout[40] + xorout[0] >= 0
 R164: - sout[32] + xorout[0] >= 0
 R165: - sout[24] + xorout[0] >= 0
 R166: - sout[16] + xorout[0] >= 0
 R167: - sout[8] + xorout[0] >= 0
 R168: - sout[8] + sout[9] >= 0
 R169: - sout[16] + sout[17] >= 0
 R170: - sout[24] + sout[25] >= 0
 R171: - sout[32] + sout[33] >= 0
 R172: - sout[40] + sout[41] >= 0
 R173: - sout[48] + sout[49] >= 0
 R174: sout[8] + sout[16] + sout[24] + sout[32] + sout[40] + sout[48]
   - xorout[0] >= 0
 R175: - xorout[0] + xorout[1] >= 0
 R176: - sout[9] - sout[17] - sout[25] - sout[33] - sout[41] - sout[49]
   + xorout[0] - xorout[1] >= -6
 R177: sout[9] + sout[17] - sout[25] - sout[33] - sout[41] - sout[49]
   + xorout[0] - xorout[1] >= -4
 R178: sout[9] - sout[17] + sout[25] - sout[33] - sout[41] - sout[49]
   + xorout[0] - xorout[1] >= -4
 R179: - sout[9] + sout[17] + sout[25] - sout[33] - sout[41] - sout[49]
   + xorout[0] - xorout[1] >= -4
 R180: sout[9] - sout[17] - sout[25] + sout[33] - sout[41] - sout[49]
   + xorout[0] - xorout[1] >= -4
 R181: - sout[9] + sout[17] - sout[25] + sout[33] - sout[41] - sout[49]
   + xorout[0] - xorout[1] >= -4
 R182: - sout[9] - sout[17] + sout[25] + sout[33] - sout[41] - sout[49]
   + xorout[0] - xorout[1] >= -4
 R183: sout[9] - sout[17] - sout[25] - sout[33] + sout[41] - sout[49]
   + xorout[0] - xorout[1] >= -4
 R184: - sout[9] + sout[17] - sout[25] - sout[33] + sout[41] - sout[49]
   + xorout[0] - xorout[1] >= -4
 R185: - sout[9] - sout[17] + sout[25] - sout[33] + sout[41] - sout[49]
   + xorout[0] - xorout[1] >= -4
 R186: - sout[9] - sout[17] - sout[25] + sout[33] + sout[41] - sout[49]
   + xorout[0] - xorout[1] >= -4
 R187: sout[9] - sout[17] - sout[25] - sout[33] - sout[41] + sout[49]
   + xorout[0] - xorout[1] >= -4
 R188: - sout[9] + sout[17] - sout[25] - sout[33] - sout[41] + sout[49]
   + xorout[0] - xorout[1] >= -4
 R189: - sout[9] - sout[17] + sout[25] - sout[33] - sout[41] + sout[49]
   + xorout[0] - xorout[1] >= -4
 R190: - sout[9] - sout[17] - sout[25] + sout[33] - sout[41] + sout[49]
   + xorout[0] - xorout[1] >= -4
 R191: - sout[9] - sout[17] - sout[25] - sout[33] + sout[41] + sout[49]
   + xorout[0] - xorout[1] >= -4
 R192: sout[9] + sout[17] + sout[25] + sout[33] - sout[41] - sout[49]
   + xorout[0] - xorout[1] >= -2
 R193: sout[9] + sout[17] + sout[25] - sout[33] + sout[41] - sout[49]
   + xorout[0] - xorout[1] >= -2
 R194: sout[9] + sout[17] - sout[25] + sout[33] + sout[41] - sout[49]
   + xorout[0] - xorout[1] >= -2
 R195: sout[9] - sout[17] + sout[25] + sout[33] + sout[41] - sout[49]
   + xorout[0] - xorout[1] >= -2
 R196: - sout[9] + sout[17] + sout[25] + sout[33] + sout[41] - sout[49]
   + xorout[0] - xorout[1] >= -2
 R197: sout[9] + sout[17] + sout[25] - sout[33] - sout[41] + sout[49]
   + xorout[0] - xorout[1] >= -2
 R198: sout[9] + sout[17] - sout[25] + sout[33] - sout[41] + sout[49]
   + xorout[0] - xorout[1] >= -2
 R199: sout[9] - sout[17] + sout[25] + sout[33] - sout[41] + sout[49]
   + xorout[0] - xorout[1] >= -2
 R200: - sout[9] + sout[17] + sout[25] + sout[33] - sout[41] + sout[49]
   + xorout[0] - xorout[1] >= -2
 R201: sout[9] + sout[17] - sout[25] - sout[33] + sout[41] + sout[49]
   + xorout[0] - xorout[1] >= -2
 R202: sout[9] - sout[17] + sout[25] - sout[33] + sout[41] + sout[49]
   + xorout[0] - xorout[1] >= -2
 R203: - sout[9] + sout[17] + sout[25] - sout[33] + sout[41] + sout[49]
   + xorout[0] - xorout[1] >= -2
 R204: sout[9] - sout[17] - sout[25] + sout[33] + sout[41] + sout[49]
   + xorout[0] - xorout[1] >= -2
 R205: - sout[9] + sout[17] - sout[25] + sout[33] + sout[41] + sout[49]
   + xorout[0] - xorout[1] >= -2
 R206: - sout[9] - sout[17] + sout[25] + sout[33] + sout[41] + sout[49]
   + xorout[0] - xorout[1] >= -2
 R207: sout[11] + sout[19] + sout[27] + sout[35] + sout[43] - sout[51]
   + xorout[3] >= 0
 R208: sout[11] + sout[19] + sout[27] + sout[35] - sout[43] + sout[51]
   + xorout[3] >= 0
 R209: sout[11] + sout[19] + sout[27] - sout[35] + sout[43] + sout[51]
   + xorout[3] >= 0
 R210: sout[11] + sout[19] - sout[27] + sout[35] + sout[43] + sout[51]
   + xorout[3] >= 0
 R211: sout[11] - sout[19] + sout[27] + sout[35] + sout[43] + sout[51]
   + xorout[3] >= 0
 R212: - sout[11] + sout[19] + sout[27] + sout[35] + sout[43] + sout[51]
   + xorout[3] >= 0
 R213: sout[11] + sout[19] + sout[27] + sout[35] + sout[43] + sout[51]
   - xorout[3] >= 0
 R214: sout[11] + sout[19] + sout[27] - sout[35] - sout[43] - sout[51]
   + xorout[3] >= -2
 R215: sout[11] + sout[19] - sout[27] + sout[35] - sout[43] - sout[51]
   + xorout[3] >= -2
 R216: sout[11] - sout[19] + sout[27] + sout[35] - sout[43] - sout[51]
   + xorout[3] >= -2
 R217: - sout[11] + sout[19] + sout[27] + sout[35] - sout[43] - sout[51]
   + xorout[3] >= -2
 R218: sout[11] + sout[19] - sout[27] - sout[35] + sout[43] - sout[51]
   + xorout[3] >= -2
 R219: sout[11] - sout[19] + sout[27] - sout[35] + sout[43] - sout[51]
   + xorout[3] >= -2
 R220: - sout[11] + sout[19] + sout[27] - sout[35] + sout[43] - sout[51]
   + xorout[3] >= -2
 R221: sout[11] - sout[19] - sout[27] + sout[35] + sout[43] - sout[51]
   + xorout[3] >= -2
 R222: - sout[11] + sout[19] - sout[27] + sout[35] + sout[43] - sout[51]
   + xorout[3] >= -2
 R223: - sout[11] - sout[19] + sout[27] + sout[35] + sout[43] - sout[51]
   + xorout[3] >= -2
 R224: sout[11] + sout[19] - sout[27] - sout[35] - sout[43] + sout[51]
   + xorout[3] >= -2
 R225: sout[11] - sout[19] + sout[27] - sout[35] - sout[43] + sout[51]
   + xorout[3] >= -2
 R226: - sout[11] + sout[19] + sout[27] - sout[35] - sout[43] + sout[51]
   + xorout[3] >= -2
 R227: sout[11] - sout[19] - sout[27] + sout[35] - sout[43] + sout[51]
   + xorout[3] >= -2
 R228: - sout[11] + sout[19] - sout[27] + sout[35] - sout[43] + sout[51]
   + xorout[3] >= -2
 R229: - sout[11] - sout[19] + sout[27] + sout[35] - sout[43] + sout[51]
   + xorout[3] >= -2
 R230: sout[11] - sout[19] - sout[27] - sout[35] + sout[43] + sout[51]
   + xorout[3] >= -2
 R231: - sout[11] + sout[19] - sout[27] - sout[35] + sout[43] + sout[51]
   + xorout[3] >= -2
 R232: - sout[11] - sout[19] + sout[27] - sout[35] + sout[43] + sout[51]
   + xorout[3] >= -2
 R233: - sout[11] - sout[19] - sout[27] + sout[35] + sout[43] + sout[51]
   + xorout[3] >= -2
 R234: sout[11] - sout[19] - sout[27] - sout[35] - sout[43] - sout[51]
   + xorout[3] >= -4
 R235: - sout[11] + sout[19] - sout[27] - sout[35] - sout[43] - sout[51]
   + xorout[3] >= -4
 R236: - sout[11] - sout[19] + sout[27] - sout[35] - sout[43] - sout[51]
   + xorout[3] >= -4
 R237: - sout[11] - sout[19] - sout[27] + sout[35] - sout[43] - sout[51]
   + xorout[3] >= -4
 R238: - sout[11] - sout[19] - sout[27] - sout[35] + sout[43] - sout[51]
   + xorout[3] >= -4
 R239: - sout[11] - sout[19] - sout[27] - sout[35] - sout[43] + sout[51]
   + xorout[3] >= -4
 R240: - sout[50] + xorout[2] >= 0
 R241: - sout[42] + xorout[2] >= 0
 R242: - sout[34] + xorout[2] >= 0
 R243: - sout[26] + xorout[2] >= 0
 R244: - sout[18] + xorout[2] >= 0
 R245: - sout[10] + xorout[2] >= 0
 R246: - sout[10] + sout[11] >= 0
 R247: - sout[18] + sout[19] >= 0
 R248: - sout[26] + sout[27] >= 0
 R249: - sout[34] + sout[35] >= 0
 R250: - sout[42] + sout[43] >= 0
 R251: - sout[50] + sout[51] >= 0
 R252: sout[10] + sout[18] + sout[26] + sout[34] + sout[42] + sout[50]
   - xorout[2] >= 0
 R253: - xorout[2] + xorout[3] >= 0
 R254: - sout[11] - sout[19] - sout[27] - sout[35] - sout[43] - sout[51]
   + xorout[2] - xorout[3] >= -6
 R255: sout[11] + sout[19] - sout[27] - sout[35] - sout[43] - sout[51]
   + xorout[2] - xorout[3] >= -4
 R256: sout[11] - sout[19] + sout[27] - sout[35] - sout[43] - sout[51]
   + xorout[2] - xorout[3] >= -4
 R257: - sout[11] + sout[19] + sout[27] - sout[35] - sout[43] - sout[51]
   + xorout[2] - xorout[3] >= -4
 R258: sout[11] - sout[19] - sout[27] + sout[35] - sout[43] - sout[51]
   + xorout[2] - xorout[3] >= -4
 R259: - sout[11] + sout[19] - sout[27] + sout[35] - sout[43] - sout[51]
   + xorout[2] - xorout[3] >= -4
 R260: - sout[11] - sout[19] + sout[27] + sout[35] - sout[43] - sout[51]
   + xorout[2] - xorout[3] >= -4
 R261: sout[11] - sout[19] - sout[27] - sout[35] + sout[43] - sout[51]
   + xorout[2] - xorout[3] >= -4
 R262: - sout[11] + sout[19] - sout[27] - sout[35] + sout[43] - sout[51]
   + xorout[2] - xorout[3] >= -4
 R263: - sout[11] - sout[19] + sout[27] - sout[35] + sout[43] - sout[51]
   + xorout[2] - xorout[3] >= -4
 R264: - sout[11] - sout[19] - sout[27] + sout[35] + sout[43] - sout[51]
   + xorout[2] - xorout[3] >= -4
 R265: sout[11] - sout[19] - sout[27] - sout[35] - sout[43] + sout[51]
   + xorout[2] - xorout[3] >= -4
 R266: - sout[11] + sout[19] - sout[27] - sout[35] - sout[43] + sout[51]
   + xorout[2] - xorout[3] >= -4
 R267: - sout[11] - sout[19] + sout[27] - sout[35] - sout[43] + sout[51]
   + xorout[2] - xorout[3] >= -4
 R268: - sout[11] - sout[19] - sout[27] + sout[35] - sout[43] + sout[51]
   + xorout[2] - xorout[3] >= -4
 R269: - sout[11] - sout[19] - sout[27] - sout[35] + sout[43] + sout[51]
   + xorout[2] - xorout[3] >= -4
 R270: sout[11] + sout[19] + sout[27] + sout[35] - sout[43] - sout[51]
   + xorout[2] - xorout[3] >= -2
 R271: sout[11] + sout[19] + sout[27] - sout[35] + sout[43] - sout[51]
   + xorout[2] - xorout[3] >= -2
 R272: sout[11] + sout[19] - sout[27] + sout[35] + sout[43] - sout[51]
   + xorout[2] - xorout[3] >= -2
 R273: sout[11] - sout[19] + sout[27] + sout[35] + sout[43] - sout[51]
   + xorout[2] - xorout[3] >= -2
 R274: - sout[11] + sout[19] + sout[27] + sout[35] + sout[43] - sout[51]
   + xorout[2] - xorout[3] >= -2
 R275: sout[11] + sout[19] + sout[27] - sout[35] - sout[43] + sout[51]
   + xorout[2] - xorout[3] >= -2
 R276: sout[11] + sout[19] - sout[27] + sout[35] - sout[43] + sout[51]
   + xorout[2] - xorout[3] >= -2
 R277: sout[11] - sout[19] + sout[27] + sout[35] - sout[43] + sout[51]
   + xorout[2] - xorout[3] >= -2
 R278: - sout[11] + sout[19] + sout[27] + sout[35] - sout[43] + sout[51]
   + xorout[2] - xorout[3] >= -2
 R279: sout[11] + sout[19] - sout[27] - sout[35] + sout[43] + sout[51]
   + xorout[2] - xorout[3] >= -2
 R280: sout[11] - sout[19] + sout[27] - sout[35] + sout[43] + sout[51]
   + xorout[2] - xorout[3] >= -2
 R281: - sout[11] + sout[19] + sout[27] - sout[35] + sout[43] + sout[51]
   + xorout[2] - xorout[3] >= -2
 R282: sout[11] - sout[19] - sout[27] + sout[35] + sout[43] + sout[51]
   + xorout[2] - xorout[3] >= -2
 R283: - sout[11] + sout[19] - sout[27] + sout[35] + sout[43] + sout[51]
   + xorout[2] - xorout[3] >= -2
 R284: - sout[11] - sout[19] + sout[27] + sout[35] + sout[43] + sout[51]
   + xorout[2] - xorout[3] >= -2
 R285: sout[13] + sout[21] + sout[29] + sout[37] + sout[45] - sout[53]
   + xorout[5] >= 0
 R286: sout[13] + sout[21] + sout[29] + sout[37] - sout[45] + sout[53]
   + xorout[5] >= 0
 R287: sout[13] + sout[21] + sout[29] - sout[37] + sout[45] + sout[53]
   + xorout[5] >= 0
 R288: sout[13] + sout[21] - sout[29] + sout[37] + sout[45] + sout[53]
   + xorout[5] >= 0
 R289: sout[13] - sout[21] + sout[29] + sout[37] + sout[45] + sout[53]
   + xorout[5] >= 0
 R290: - sout[13] + sout[21] + sout[29] + sout[37] + sout[45] + sout[53]
   + xorout[5] >= 0
 R291: sout[13] + sout[21] + sout[29] + sout[37] + sout[45] + sout[53]
   - xorout[5] >= 0
 R292: sout[13] + sout[21] + sout[29] - sout[37] - sout[45] - sout[53]
   + xorout[5] >= -2
 R293: sout[13] + sout[21] - sout[29] + sout[37] - sout[45] - sout[53]
   + xorout[5] >= -2
 R294: sout[13] - sout[21] + sout[29] + sout[37] - sout[45] - sout[53]
   + xorout[5] >= -2
 R295: - sout[13] + sout[21] + sout[29] + sout[37] - sout[45] - sout[53]
   + xorout[5] >= -2
 R296: sout[13] + sout[21] - sout[29] - sout[37] + sout[45] - sout[53]
   + xorout[5] >= -2
 R297: sout[13] - sout[21] + sout[29] - sout[37] + sout[45] - sout[53]
   + xorout[5] >= -2
 R298: - sout[13] + sout[21] + sout[29] - sout[37] + sout[45] - sout[53]
   + xorout[5] >= -2
 R299: sout[13] - sout[21] - sout[29] + sout[37] + sout[45] - sout[53]
   + xorout[5] >= -2
 R300: - sout[13] + sout[21] - sout[29] + sout[37] + sout[45] - sout[53]
   + xorout[5] >= -2
 R301: - sout[13] - sout[21] + sout[29] + sout[37] + sout[45] - sout[53]
   + xorout[5] >= -2
 R302: sout[13] + sout[21] - sout[29] - sout[37] - sout[45] + sout[53]
   + xorout[5] >= -2
 R303: sout[13] - sout[21] + sout[29] - sout[37] - sout[45] + sout[53]
   + xorout[5] >= -2
 R304: - sout[13] + sout[21] + sout[29] - sout[37] - sout[45] + sout[53]
   + xorout[5] >= -2
 R305: sout[13] - sout[21] - sout[29] + sout[37] - sout[45] + sout[53]
   + xorout[5] >= -2
 R306: - sout[13] + sout[21] - sout[29] + sout[37] - sout[45] + sout[53]
   + xorout[5] >= -2
 R307: - sout[13] - sout[21] + sout[29] + sout[37] - sout[45] + sout[53]
   + xorout[5] >= -2
 R308: sout[13] - sout[21] - sout[29] - sout[37] + sout[45] + sout[53]
   + xorout[5] >= -2
 R309: - sout[13] + sout[21] - sout[29] - sout[37] + sout[45] + sout[53]
   + xorout[5] >= -2
 R310: - sout[13] - sout[21] + sout[29] - sout[37] + sout[45] + sout[53]
   + xorout[5] >= -2
 R311: - sout[13] - sout[21] - sout[29] + sout[37] + sout[45] + sout[53]
   + xorout[5] >= -2
 R312: sout[13] - sout[21] - sout[29] - sout[37] - sout[45] - sout[53]
   + xorout[5] >= -4
 R313: - sout[13] + sout[21] - sout[29] - sout[37] - sout[45] - sout[53]
   + xorout[5] >= -4
 R314: - sout[13] - sout[21] + sout[29] - sout[37] - sout[45] - sout[53]
   + xorout[5] >= -4
 R315: - sout[13] - sout[21] - sout[29] + sout[37] - sout[45] - sout[53]
   + xorout[5] >= -4
 R316: - sout[13] - sout[21] - sout[29] - sout[37] + sout[45] - sout[53]
   + xorout[5] >= -4
 R317: - sout[13] - sout[21] - sout[29] - sout[37] - sout[45] + sout[53]
   + xorout[5] >= -4
 R318: - sout[52] + xorout[4] >= 0
 R319: - sout[44] + xorout[4] >= 0
 R320: - sout[36] + xorout[4] >= 0
 R321: - sout[28] + xorout[4] >= 0
 R322: - sout[20] + xorout[4] >= 0
 R323: - sout[12] + xorout[4] >= 0
 R324: - sout[12] + sout[13] >= 0
 R325: - sout[20] + sout[21] >= 0
 R326: - sout[28] + sout[29] >= 0
 R327: - sout[36] + sout[37] >= 0
 R328: - sout[44] + sout[45] >= 0
 R329: - sout[52] + sout[53] >= 0
 R330: sout[12] + sout[20] + sout[28] + sout[36] + sout[44] + sout[52]
   - xorout[4] >= 0
 R331: - xorout[4] + xorout[5] >= 0
 R332: - sout[13] - sout[21] - sout[29] - sout[37] - sout[45] - sout[53]
   + xorout[4] - xorout[5] >= -6
 R333: sout[13] + sout[21] - sout[29] - sout[37] - sout[45] - sout[53]
   + xorout[4] - xorout[5] >= -4
 R334: sout[13] - sout[21] + sout[29] - sout[37] - sout[45] - sout[53]
   + xorout[4] - xorout[5] >= -4
 R335: - sout[13] + sout[21] + sout[29] - sout[37] - sout[45] - sout[53]
   + xorout[4] - xorout[5] >= -4
 R336: sout[13] - sout[21] - sout[29] + sout[37] - sout[45] - sout[53]
   + xorout[4] - xorout[5] >= -4
 R337: - sout[13] + sout[21] - sout[29] + sout[37] - sout[45] - sout[53]
   + xorout[4] - xorout[5] >= -4
 R338: - sout[13] - sout[21] + sout[29] + sout[37] - sout[45] - sout[53]
   + xorout[4] - xorout[5] >= -4
 R339: sout[13] - sout[21] - sout[29] - sout[37] + sout[45] - sout[53]
   + xorout[4] - xorout[5] >= -4
 R340: - sout[13] + sout[21] - sout[29] - sout[37] + sout[45] - sout[53]
   + xorout[4] - xorout[5] >= -4
 R341: - sout[13] - sout[21] + sout[29] - sout[37] + sout[45] - sout[53]
   + xorout[4] - xorout[5] >= -4
 R342: - sout[13] - sout[21] - sout[29] + sout[37] + sout[45] - sout[53]
   + xorout[4] - xorout[5] >= -4
 R343: sout[13] - sout[21] - sout[29] - sout[37] - sout[45] + sout[53]
   + xorout[4] - xorout[5] >= -4
 R344: - sout[13] + sout[21] - sout[29] - sout[37] - sout[45] + sout[53]
   + xorout[4] - xorout[5] >= -4
 R345: - sout[13] - sout[21] + sout[29] - sout[37] - sout[45] + sout[53]
   + xorout[4] - xorout[5] >= -4
 R346: - sout[13] - sout[21] - sout[29] + sout[37] - sout[45] + sout[53]
   + xorout[4] - xorout[5] >= -4
 R347: - sout[13] - sout[21] - sout[29] - sout[37] + sout[45] + sout[53]
   + xorout[4] - xorout[5] >= -4
 R348: sout[13] + sout[21] + sout[29] + sout[37] - sout[45] - sout[53]
   + xorout[4] - xorout[5] >= -2
 R349: sout[13] + sout[21] + sout[29] - sout[37] + sout[45] - sout[53]
   + xorout[4] - xorout[5] >= -2
 R350: sout[13] + sout[21] - sout[29] + sout[37] + sout[45] - sout[53]
   + xorout[4] - xorout[5] >= -2
 R351: sout[13] - sout[21] + sout[29] + sout[37] + sout[45] - sout[53]
   + xorout[4] - xorout[5] >= -2
 R352: - sout[13] + sout[21] + sout[29] + sout[37] + sout[45] - sout[53]
   + xorout[4] - xorout[5] >= -2
 R353: sout[13] + sout[21] + sout[29] - sout[37] - sout[45] + sout[53]
   + xorout[4] - xorout[5] >= -2
 R354: sout[13] + sout[21] - sout[29] + sout[37] - sout[45] + sout[53]
   + xorout[4] - xorout[5] >= -2
 R355: sout[13] - sout[21] + sout[29] + sout[37] - sout[45] + sout[53]
   + xorout[4] - xorout[5] >= -2
 R356: - sout[13] + sout[21] + sout[29] + sout[37] - sout[45] + sout[53]
   + xorout[4] - xorout[5] >= -2
 R357: sout[13] + sout[21] - sout[29] - sout[37] + sout[45] + sout[53]
   + xorout[4] - xorout[5] >= -2
 R358: sout[13] - sout[21] + sout[29] - sout[37] + sout[45] + sout[53]
   + xorout[4] - xorout[5] >= -2
 R359: - sout[13] + sout[21] + sout[29] - sout[37] + sout[45] + sout[53]
   + xorout[4] - xorout[5] >= -2
 R360: sout[13] - sout[21] - sout[29] + sout[37] + sout[45] + sout[53]
   + xorout[4] - xorout[5] >= -2
 R361: - sout[13] + sout[21] - sout[29] + sout[37] + sout[45] + sout[53]
   + xorout[4] - xorout[5] >= -2
 R362: - sout[13] - sout[21] + sout[29] + sout[37] + sout[45] + sout[53]
   + xorout[4] - xorout[5] >= -2
 R363: sout[15] + sout[23] + sout[31] + sout[39] + sout[47] - sout[55]
   + xorout[7] >= 0
 R364: sout[15] + sout[23] + sout[31] + sout[39] - sout[47] + sout[55]
   + xorout[7] >= 0
 R365: sout[15] + sout[23] + sout[31] - sout[39] + sout[47] + sout[55]
   + xorout[7] >= 0
 R366: sout[15] + sout[23] - sout[31] + sout[39] + sout[47] + sout[55]
   + xorout[7] >= 0
 R367: sout[15] - sout[23] + sout[31] + sout[39] + sout[47] + sout[55]
   + xorout[7] >= 0
 R368: - sout[15] + sout[23] + sout[31] + sout[39] + sout[47] + sout[55]
   + xorout[7] >= 0
 R369: sout[15] + sout[23] + sout[31] + sout[39] + sout[47] + sout[55]
   - xorout[7] >= 0
 R370: sout[15] + sout[23] + sout[31] - sout[39] - sout[47] - sout[55]
   + xorout[7] >= -2
 R371: sout[15] + sout[23] - sout[31] + sout[39] - sout[47] - sout[55]
   + xorout[7] >= -2
 R372: sout[15] - sout[23] + sout[31] + sout[39] - sout[47] - sout[55]
   + xorout[7] >= -2
 R373: - sout[15] + sout[23] + sout[31] + sout[39] - sout[47] - sout[55]
   + xorout[7] >= -2
 R374: sout[15] + sout[23] - sout[31] - sout[39] + sout[47] - sout[55]
   + xorout[7] >= -2
 R375: sout[15] - sout[23] + sout[31] - sout[39] + sout[47] - sout[55]
   + xorout[7] >= -2
 R376: - sout[15] + sout[23] + sout[31] - sout[39] + sout[47] - sout[55]
   + xorout[7] >= -2
 R377: sout[15] - sout[23] - sout[31] + sout[39] + sout[47] - sout[55]
   + xorout[7] >= -2
 R378: - sout[15] + sout[23] - sout[31] + sout[39] + sout[47] - sout[55]
   + xorout[7] >= -2
 R379: - sout[15] - sout[23] + sout[31] + sout[39] + sout[47] - sout[55]
   + xorout[7] >= -2
 R380: sout[15] + sout[23] - sout[31] - sout[39] - sout[47] + sout[55]
   + xorout[7] >= -2
 R381: sout[15] - sout[23] + sout[31] - sout[39] - sout[47] + sout[55]
   + xorout[7] >= -2
 R382: - sout[15] + sout[23] + sout[31] - sout[39] - sout[47] + sout[55]
   + xorout[7] >= -2
 R383: sout[15] - sout[23] - sout[31] + sout[39] - sout[47] + sout[55]
   + xorout[7] >= -2
 R384: - sout[15] + sout[23] - sout[31] + sout[39] - sout[47] + sout[55]
   + xorout[7] >= -2
 R385: - sout[15] - sout[23] + sout[31] + sout[39] - sout[47] + sout[55]
   + xorout[7] >= -2
 R386: sout[15] - sout[23] - sout[31] - sout[39] + sout[47] + sout[55]
   + xorout[7] >= -2
 R387: - sout[15] + sout[23] - sout[31] - sout[39] + sout[47] + sout[55]
   + xorout[7] >= -2
 R388: - sout[15] - sout[23] + sout[31] - sout[39] + sout[47] + sout[55]
   + xorout[7] >= -2
 R389: - sout[15] - sout[23] - sout[31] + sout[39] + sout[47] + sout[55]
   + xorout[7] >= -2
 R390: sout[15] - sout[23] - sout[31] - sout[39] - sout[47] - sout[55]
   + xorout[7] >= -4
 R391: - sout[15] + sout[23] - sout[31] - sout[39] - sout[47] - sout[55]
   + xorout[7] >= -4
 R392: - sout[15] - sout[23] + sout[31] - sout[39] - sout[47] - sout[55]
   + xorout[7] >= -4
 R393: - sout[15] - sout[23] - sout[31] + sout[39] - sout[47] - sout[55]
   + xorout[7] >= -4
 R394: - sout[15] - sout[23] - sout[31] - sout[39] + sout[47] - sout[55]
   + xorout[7] >= -4
 R395: - sout[15] - sout[23] - sout[31] - sout[39] - sout[47] + sout[55]
   + xorout[7] >= -4
 R396: - sout[54] + xorout[6] >= 0
 R397: - sout[46] + xorout[6] >= 0
 R398: - sout[38] + xorout[6] >= 0
 R399: - sout[30] + xorout[6] >= 0
 R400: - sout[22] + xorout[6] >= 0
 R401: - sout[14] + xorout[6] >= 0
 R402: - sout[14] + sout[15] >= 0
 R403: - sout[22] + sout[23] >= 0
 R404: - sout[30] + sout[31] >= 0
 R405: - sout[38] + sout[39] >= 0
 R406: - sout[46] + sout[47] >= 0
 R407: - sout[54] + sout[55] >= 0
 R408: sout[14] + sout[22] + sout[30] + sout[38] + sout[46] + sout[54]
   - xorout[6] >= 0
 R409: - xorout[6] + xorout[7] >= 0
 R410: - sout[15] - sout[23] - sout[31] - sout[39] - sout[47] - sout[55]
   + xorout[6] - xorout[7] >= -6
 R411: sout[15] + sout[23] - sout[31] - sout[39] - sout[47] - sout[55]
   + xorout[6] - xorout[7] >= -4
 R412: sout[15] - sout[23] + sout[31] - sout[39] - sout[47] - sout[55]
   + xorout[6] - xorout[7] >= -4
 R413: - sout[15] + sout[23] + sout[31] - sout[39] - sout[47] - sout[55]
   + xorout[6] - xorout[7] >= -4
 R414: sout[15] - sout[23] - sout[31] + sout[39] - sout[47] - sout[55]
   + xorout[6] - xorout[7] >= -4
 R415: - sout[15] + sout[23] - sout[31] + sout[39] - sout[47] - sout[55]
   + xorout[6] - xorout[7] >= -4
 R416: - sout[15] - sout[23] + sout[31] + sout[39] - sout[47] - sout[55]
   + xorout[6] - xorout[7] >= -4
 R417: sout[15] - sout[23] - sout[31] - sout[39] + sout[47] - sout[55]
   + xorout[6] - xorout[7] >= -4
 R418: - sout[15] + sout[23] - sout[31] - sout[39] + sout[47] - sout[55]
   + xorout[6] - xorout[7] >= -4
 R419: - sout[15] - sout[23] + sout[31] - sout[39] + sout[47] - sout[55]
   + xorout[6] - xorout[7] >= -4
 R420: - sout[15] - sout[23] - sout[31] + sout[39] + sout[47] - sout[55]
   + xorout[6] - xorout[7] >= -4
 R421: sout[15] - sout[23] - sout[31] - sout[39] - sout[47] + sout[55]
   + xorout[6] - xorout[7] >= -4
 R422: - sout[15] + sout[23] - sout[31] - sout[39] - sout[47] + sout[55]
   + xorout[6] - xorout[7] >= -4
 R423: - sout[15] - sout[23] + sout[31] - sout[39] - sout[47] + sout[55]
   + xorout[6] - xorout[7] >= -4
 R424: - sout[15] - sout[23] - sout[31] + sout[39] - sout[47] + sout[55]
   + xorout[6] - xorout[7] >= -4
 R425: - sout[15] - sout[23] - sout[31] - sout[39] + sout[47] + sout[55]
   + xorout[6] - xorout[7] >= -4
 R426: sout[15] + sout[23] + sout[31] + sout[39] - sout[47] - sout[55]
   + xorout[6] - xorout[7] >= -2
 R427: sout[15] + sout[23] + sout[31] - sout[39] + sout[47] - sout[55]
   + xorout[6] - xorout[7] >= -2
 R428: sout[15] + sout[23] - sout[31] + sout[39] + sout[47] - sout[55]
   + xorout[6] - xorout[7] >= -2
 R429: sout[15] - sout[23] + sout[31] + sout[39] + sout[47] - sout[55]
   + xorout[6] - xorout[7] >= -2
 R430: - sout[15] + sout[23] + sout[31] + sout[39] + sout[47] - sout[55]
   + xorout[6] - xorout[7] >= -2
 R431: sout[15] + sout[23] + sout[31] - sout[39] - sout[47] + sout[55]
   + xorout[6] - xorout[7] >= -2
 R432: sout[15] + sout[23] - sout[31] + sout[39] - sout[47] + sout[55]
   + xorout[6] - xorout[7] >= -2
 R433: sout[15] - sout[23] + sout[31] + sout[39] - sout[47] + sout[55]
   + xorout[6] - xorout[7] >= -2
 R434: - sout[15] + sout[23] + sout[31] + sout[39] - sout[47] + sout[55]
   + xorout[6] - xorout[7] >= -2
 R435: sout[15] + sout[23] - sout[31] - sout[39] + sout[47] + sout[55]
   + xorout[6] - xorout[7] >= -2
 R436: sout[15] - sout[23] + sout[31] - sout[39] + sout[47] + sout[55]
   + xorout[6] - xorout[7] >= -2
 R437: - sout[15] + sout[23] + sout[31] - sout[39] + sout[47] + sout[55]
   + xorout[6] - xorout[7] >= -2
 R438: sout[15] - sout[23] - sout[31] + sout[39] + sout[47] + sout[55]
   + xorout[6] - xorout[7] >= -2
 R439: - sout[15] + sout[23] - sout[31] + sout[39] + sout[47] + sout[55]
   + xorout[6] - xorout[7] >= -2
 R440: - sout[15] - sout[23] + sout[31] + sout[39] + sout[47] + sout[55]
   + xorout[6] - xorout[7] >= -2
 R441: sout[1] + sout[17] + sout[25] + sout[41] + sout[49] - sout[57]
   + xorout[9] >= 0
 R442: sout[1] + sout[17] + sout[25] + sout[41] - sout[49] + sout[57]
   + xorout[9] >= 0
 R443: sout[1] + sout[17] + sout[25] - sout[41] + sout[49] + sout[57]
   + xorout[9] >= 0
 R444: sout[1] + sout[17] - sout[25] + sout[41] + sout[49] + sout[57]
   + xorout[9] >= 0
 R445: sout[1] - sout[17] + sout[25] + sout[41] + sout[49] + sout[57]
   + xorout[9] >= 0
 R446: - sout[1] + sout[17] + sout[25] + sout[41] + sout[49] + sout[57]
   + xorout[9] >= 0
 R447: sout[1] + sout[17] + sout[25] + sout[41] + sout[49] + sout[57]
   - xorout[9] >= 0
 R448: sout[1] + sout[17] + sout[25] - sout[41] - sout[49] - sout[57]
   + xorout[9] >= -2
 R449: sout[1] + sout[17] - sout[25] + sout[41] - sout[49] - sout[57]
   + xorout[9] >= -2
 R450: sout[1] - sout[17] + sout[25] + sout[41] - sout[49] - sout[57]
   + xorout[9] >= -2
 R451: - sout[1] + sout[17] + sout[25] + sout[41] - sout[49] - sout[57]
   + xorout[9] >= -2
 R452: sout[1] + sout[17] - sout[25] - sout[41] + sout[49] - sout[57]
   + xorout[9] >= -2
 R453: sout[1] - sout[17] + sout[25] - sout[41] + sout[49] - sout[57]
   + xorout[9] >= -2
 R454: - sout[1] + sout[17] + sout[25] - sout[41] + sout[49] - sout[57]
   + xorout[9] >= -2
 R455: sout[1] - sout[17] - sout[25] + sout[41] + sout[49] - sout[57]
   + xorout[9] >= -2
 R456: - sout[1] + sout[17] - sout[25] + sout[41] + sout[49] - sout[57]
   + xorout[9] >= -2
 R457: - sout[1] - sout[17] + sout[25] + sout[41] + sout[49] - sout[57]
   + xorout[9] >= -2
 R458: sout[1] + sout[17] - sout[25] - sout[41] - sout[49] + sout[57]
   + xorout[9] >= -2
 R459: sout[1] - sout[17] + sout[25] - sout[41] - sout[49] + sout[57]
   + xorout[9] >= -2
 R460: - sout[1] + sout[17] + sout[25] - sout[41] - sout[49] + sout[57]
   + xorout[9] >= -2
 R461: sout[1] - sout[17] - sout[25] + sout[41] - sout[49] + sout[57]
   + xorout[9] >= -2
 R462: - sout[1] + sout[17] - sout[25] + sout[41] - sout[49] + sout[57]
   + xorout[9] >= -2
 R463: - sout[1] - sout[17] + sout[25] + sout[41] - sout[49] + sout[57]
   + xorout[9] >= -2
 R464: sout[1] - sout[17] - sout[25] - sout[41] + sout[49] + sout[57]
   + xorout[9] >= -2
 R465: - sout[1] + sout[17] - sout[25] - sout[41] + sout[49] + sout[57]
   + xorout[9] >= -2
 R466: - sout[1] - sout[17] + sout[25] - sout[41] + sout[49] + sout[57]
   + xorout[9] >= -2
 R467: - sout[1] - sout[17] - sout[25] + sout[41] + sout[49] + sout[57]
   + xorout[9] >= -2
 R468: sout[1] - sout[17] - sout[25] - sout[41] - sout[49] - sout[57]
   + xorout[9] >= -4
 R469: - sout[1] + sout[17] - sout[25] - sout[41] - sout[49] - sout[57]
   + xorout[9] >= -4
 R470: - sout[1] - sout[17] + sout[25] - sout[41] - sout[49] - sout[57]
   + xorout[9] >= -4
 R471: - sout[1] - sout[17] - sout[25] + sout[41] - sout[49] - sout[57]
   + xorout[9] >= -4
 R472: - sout[1] - sout[17] - sout[25] - sout[41] + sout[49] - sout[57]
   + xorout[9] >= -4
 R473: - sout[1] - sout[17] - sout[25] - sout[41] - sout[49] + sout[57]
   + xorout[9] >= -4
 R474: - sout[56] + xorout[8] >= 0
 R475: - sout[48] + xorout[8] >= 0
 R476: - sout[40] + xorout[8] >= 0
 R477: - sout[24] + xorout[8] >= 0
 R478: - sout[16] + xorout[8] >= 0
 R479: - sout[0] + xorout[8] >= 0
 R480: - sout[0] + sout[1] >= 0
 R481: - sout[16] + sout[17] >= 0
 R482: - sout[24] + sout[25] >= 0
 R483: - sout[40] + sout[41] >= 0
 R484: - sout[48] + sout[49] >= 0
 R485: - sout[56] + sout[57] >= 0
 R486: sout[0] + sout[16] + sout[24] + sout[40] + sout[48] + sout[56]
   - xorout[8] >= 0
 R487: - xorout[8] + xorout[9] >= 0
 R488: - sout[1] - sout[17] - sout[25] - sout[41] - sout[49] - sout[57]
   + xorout[8] - xorout[9] >= -6
 R489: sout[1] + sout[17] - sout[25] - sout[41] - sout[49] - sout[57]
   + xorout[8] - xorout[9] >= -4
 R490: sout[1] - sout[17] + sout[25] - sout[41] - sout[49] - sout[57]
   + xorout[8] - xorout[9] >= -4
 R491: - sout[1] + sout[17] + sout[25] - sout[41] - sout[49] - sout[57]
   + xorout[8] - xorout[9] >= -4
 R492: sout[1] - sout[17] - sout[25] + sout[41] - sout[49] - sout[57]
   + xorout[8] - xorout[9] >= -4
 R493: - sout[1] + sout[17] - sout[25] + sout[41] - sout[49] - sout[57]
   + xorout[8] - xorout[9] >= -4
 R494: - sout[1] - sout[17] + sout[25] + sout[41] - sout[49] - sout[57]
   + xorout[8] - xorout[9] >= -4
 R495: sout[1] - sout[17] - sout[25] - sout[41] + sout[49] - sout[57]
   + xorout[8] - xorout[9] >= -4
 R496: - sout[1] + sout[17] - sout[25] - sout[41] + sout[49] - sout[57]
   + xorout[8] - xorout[9] >= -4
 R497: - sout[1] - sout[17] + sout[25] - sout[41] + sout[49] - sout[57]
   + xorout[8] - xorout[9] >= -4
 R498: - sout[1] - sout[17] - sout[25] + sout[41] + sout[49] - sout[57]
   + xorout[8] - xorout[9] >= -4
 R499: sout[1] - sout[17] - sout[25] - sout[41] - sout[49] + sout[57]
   + xorout[8] - xorout[9] >= -4
 R500: - sout[1] + sout[17] - sout[25] - sout[41] - sout[49] + sout[57]
   + xorout[8] - xorout[9] >= -4
 R501: - sout[1] - sout[17] + sout[25] - sout[41] - sout[49] + sout[57]
   + xorout[8] - xorout[9] >= -4
 R502: - sout[1] - sout[17] - sout[25] + sout[41] - sout[49] + sout[57]
   + xorout[8] - xorout[9] >= -4
 R503: - sout[1] - sout[17] - sout[25] - sout[41] + sout[49] + sout[57]
   + xorout[8] - xorout[9] >= -4
 R504: sout[1] + sout[17] + sout[25] + sout[41] - sout[49] - sout[57]
   + xorout[8] - xorout[9] >= -2
 R505: sout[1] + sout[17] + sout[25] - sout[41] + sout[49] - sout[57]
   + xorout[8] - xorout[9] >= -2
 R506: sout[1] + sout[17] - sout[25] + sout[41] + sout[49] - sout[57]
   + xorout[8] - xorout[9] >= -2
 R507: sout[1] - sout[17] + sout[25] + sout[41] + sout[49] - sout[57]
   + xorout[8] - xorout[9] >= -2
 R508: - sout[1] + sout[17] + sout[25] + sout[41] + sout[49] - sout[57]
   + xorout[8] - xorout[9] >= -2
 R509: sout[1] + sout[17] + sout[25] - sout[41] - sout[49] + sout[57]
   + xorout[8] - xorout[9] >= -2
 R510: sout[1] + sout[17] - sout[25] + sout[41] - sout[49] + sout[57]
   + xorout[8] - xorout[9] >= -2
 R511: sout[1] - sout[17] + sout[25] + sout[41] - sout[49] + sout[57]
   + xorout[8] - xorout[9] >= -2
 R512: - sout[1] + sout[17] + sout[25] + sout[41] - sout[49] + sout[57]
   + xorout[8] - xorout[9] >= -2
 R513: sout[1] + sout[17] - sout[25] - sout[41] + sout[49] + sout[57]
   + xorout[8] - xorout[9] >= -2
 R514: sout[1] - sout[17] + sout[25] - sout[41] + sout[49] + sout[57]
   + xorout[8] - xorout[9] >= -2
 R515: - sout[1] + sout[17] + sout[25] - sout[41] + sout[49] + sout[57]
   + xorout[8] - xorout[9] >= -2
 R516: sout[1] - sout[17] - sout[25] + sout[41] + sout[49] + sout[57]
   + xorout[8] - xorout[9] >= -2
 R517: - sout[1] + sout[17] - sout[25] + sout[41] + sout[49] + sout[57]
   + xorout[8] - xorout[9] >= -2
 R518: - sout[1] - sout[17] + sout[25] + sout[41] + sout[49] + sout[57]
   + xorout[8] - xorout[9] >= -2
 R519: sout[3] + sout[19] + sout[27] + sout[43] + sout[51] - sout[59]
   + xorout[11] >= 0
 R520: sout[3] + sout[19] + sout[27] + sout[43] - sout[51] + sout[59]
   + xorout[11] >= 0
 R521: sout[3] + sout[19] + sout[27] - sout[43] + sout[51] + sout[59]
   + xorout[11] >= 0
 R522: sout[3] + sout[19] - sout[27] + sout[43] + sout[51] + sout[59]
   + xorout[11] >= 0
 R523: sout[3] - sout[19] + sout[27] + sout[43] + sout[51] + sout[59]
   + xorout[11] >= 0
 R524: - sout[3] + sout[19] + sout[27] + sout[43] + sout[51] + sout[59]
   + xorout[11] >= 0
 R525: sout[3] + sout[19] + sout[27] + sout[43] + sout[51] + sout[59]
   - xorout[11] >= 0
 R526: sout[3] + sout[19] + sout[27] - sout[43] - sout[51] - sout[59]
   + xorout[11] >= -2
 R527: sout[3] + sout[19] - sout[27] + sout[43] - sout[51] - sout[59]
   + xorout[11] >= -2
 R528: sout[3] - sout[19] + sout[27] + sout[43] - sout[51] - sout[59]
   + xorout[11] >= -2
 R529: - sout[3] + sout[19] + sout[27] + sout[43] - sout[51] - sout[59]
   + xorout[11] >= -2
 R530: sout[3] + sout[19] - sout[27] - sout[43] + sout[51] - sout[59]
   + xorout[11] >= -2
 R531: sout[3] - sout[19] + sout[27] - sout[43] + sout[51] - sout[59]
   + xorout[11] >= -2
 R532: - sout[3] + sout[19] + sout[27] - sout[43] + sout[51] - sout[59]
   + xorout[11] >= -2
 R533: sout[3] - sout[19] - sout[27] + sout[43] + sout[51] - sout[59]
   + xorout[11] >= -2
 R534: - sout[3] + sout[19] - sout[27] + sout[43] + sout[51] - sout[59]
   + xorout[11] >= -2
 R535: - sout[3] - sout[19] + sout[27] + sout[43] + sout[51] - sout[59]
   + xorout[11] >= -2
 R536: sout[3] + sout[19] - sout[27] - sout[43] - sout[51] + sout[59]
   + xorout[11] >= -2
 R537: sout[3] - sout[19] + sout[27] - sout[43] - sout[51] + sout[59]
   + xorout[11] >= -2
 R538: - sout[3] + sout[19] + sout[27] - sout[43] - sout[51] + sout[59]
   + xorout[11] >= -2
 R539: sout[3] - sout[19] - sout[27] + sout[43] - sout[51] + sout[59]
   + xorout[11] >= -2
 R540: - sout[3] + sout[19] - sout[27] + sout[43] - sout[51] + sout[59]
   + xorout[11] >= -2
 R541: - sout[3] - sout[19] + sout[27] + sout[43] - sout[51] + sout[59]
   + xorout[11] >= -2
 R542: sout[3] - sout[19] - sout[27] - sout[43] + sout[51] + sout[59]
   + xorout[11] >= -2
 R543: - sout[3] + sout[19] - sout[27] - sout[43] + sout[51] + sout[59]
   + xorout[11] >= -2
 R544: - sout[3] - sout[19] + sout[27] - sout[43] + sout[51] + sout[59]
   + xorout[11] >= -2
 R545: - sout[3] - sout[19] - sout[27] + sout[43] + sout[51] + sout[59]
   + xorout[11] >= -2
 R546: sout[3] - sout[19] - sout[27] - sout[43] - sout[51] - sout[59]
   + xorout[11] >= -4
 R547: - sout[3] + sout[19] - sout[27] - sout[43] - sout[51] - sout[59]
   + xorout[11] >= -4
 R548: - sout[3] - sout[19] + sout[27] - sout[43] - sout[51] - sout[59]
   + xorout[11] >= -4
 R549: - sout[3] - sout[19] - sout[27] + sout[43] - sout[51] - sout[59]
   + xorout[11] >= -4
 R550: - sout[3] - sout[19] - sout[27] - sout[43] + sout[51] - sout[59]
   + xorout[11] >= -4
 R551: - sout[3] - sout[19] - sout[27] - sout[43] - sout[51] + sout[59]
   + xorout[11] >= -4
 R552: - sout[58] + xorout[10] >= 0
 R553: - sout[50] + xorout[10] >= 0
 R554: - sout[42] + xorout[10] >= 0
 R555: - sout[26] + xorout[10] >= 0
 R556: - sout[18] + xorout[10] >= 0
 R557: - sout[2] + xorout[10] >= 0
 R558: - sout[2] + sout[3] >= 0
 R559: - sout[18] + sout[19] >= 0
 R560: - sout[26] + sout[27] >= 0
 R561: - sout[42] + sout[43] >= 0
 R562: - sout[50] + sout[51] >= 0
 R563: - sout[58] + sout[59] >= 0
 R564: sout[2] + sout[18] + sout[26] + sout[42] + sout[50] + sout[58]
   - xorout[10] >= 0
 R565: - xorout[10] + xorout[11] >= 0
 R566: - sout[3] - sout[19] - sout[27] - sout[43] - sout[51] - sout[59]
   + xorout[10] - xorout[11] >= -6
 R567: sout[3] + sout[19] - sout[27] - sout[43] - sout[51] - sout[59]
   + xorout[10] - xorout[11] >= -4
 R568: sout[3] - sout[19] + sout[27] - sout[43] - sout[51] - sout[59]
   + xorout[10] - xorout[11] >= -4
 R569: - sout[3] + sout[19] + sout[27] - sout[43] - sout[51] - sout[59]
   + xorout[10] - xorout[11] >= -4
 R570: sout[3] - sout[19] - sout[27] + sout[43] - sout[51] - sout[59]
   + xorout[10] - xorout[11] >= -4
 R571: - sout[3] + sout[19] - sout[27] + sout[43] - sout[51] - sout[59]
   + xorout[10] - xorout[11] >= -4
 R572: - sout[3] - sout[19] + sout[27] + sout[43] - sout[51] - sout[59]
   + xorout[10] - xorout[11] >= -4
 R573: sout[3] - sout[19] - sout[27] - sout[43] + sout[51] - sout[59]
   + xorout[10] - xorout[11] >= -4
 R574: - sout[3] + sout[19] - sout[27] - sout[43] + sout[51] - sout[59]
   + xorout[10] - xorout[11] >= -4
 R575: - sout[3] - sout[19] + sout[27] - sout[43] + sout[51] - sout[59]
   + xorout[10] - xorout[11] >= -4
 R576: - sout[3] - sout[19] - sout[27] + sout[43] + sout[51] - sout[59]
   + xorout[10] - xorout[11] >= -4
 R577: sout[3] - sout[19] - sout[27] - sout[43] - sout[51] + sout[59]
   + xorout[10] - xorout[11] >= -4
 R578: - sout[3] + sout[19] - sout[27] - sout[43] - sout[51] + sout[59]
   + xorout[10] - xorout[11] >= -4
 R579: - sout[3] - sout[19] + sout[27] - sout[43] - sout[51] + sout[59]
   + xorout[10] - xorout[11] >= -4
 R580: - sout[3] - sout[19] - sout[27] + sout[43] - sout[51] + sout[59]
   + xorout[10] - xorout[11] >= -4
 R581: - sout[3] - sout[19] - sout[27] - sout[43] + sout[51] + sout[59]
   + xorout[10] - xorout[11] >= -4
 R582: sout[3] + sout[19] + sout[27] + sout[43] - sout[51] - sout[59]
   + xorout[10] - xorout[11] >= -2
 R583: sout[3] + sout[19] + sout[27] - sout[43] + sout[51] - sout[59]
   + xorout[10] - xorout[11] >= -2
 R584: sout[3] + sout[19] - sout[27] + sout[43] + sout[51] - sout[59]
   + xorout[10] - xorout[11] >= -2
 R585: sout[3] - sout[19] + sout[27] + sout[43] + sout[51] - sout[59]
   + xorout[10] - xorout[11] >= -2
 R586: - sout[3] + sout[19] + sout[27] + sout[43] + sout[51] - sout[59]
   + xorout[10] - xorout[11] >= -2
 R587: sout[3] + sout[19] + sout[27] - sout[43] - sout[51] + sout[59]
   + xorout[10] - xorout[11] >= -2
 R588: sout[3] + sout[19] - sout[27] + sout[43] - sout[51] + sout[59]
   + xorout[10] - xorout[11] >= -2
 R589: sout[3] - sout[19] + sout[27] + sout[43] - sout[51] + sout[59]
   + xorout[10] - xorout[11] >= -2
 R590: - sout[3] + sout[19] + sout[27] + sout[43] - sout[51] + sout[59]
   + xorout[10] - xorout[11] >= -2
 R591: sout[3] + sout[19] - sout[27] - sout[43] + sout[51] + sout[59]
   + xorout[10] - xorout[11] >= -2
 R592: sout[3] - sout[19] + sout[27] - sout[43] + sout[51] + sout[59]
   + xorout[10] - xorout[11] >= -2
 R593: - sout[3] + sout[19] + sout[27] - sout[43] + sout[51] + sout[59]
   + xorout[10] - xorout[11] >= -2
 R594: sout[3] - sout[19] - sout[27] + sout[43] + sout[51] + sout[59]
   + xorout[10] - xorout[11] >= -2
 R595: - sout[3] + sout[19] - sout[27] + sout[43] + sout[51] + sout[59]
   + xorout[10] - xorout[11] >= -2
 R596: - sout[3] - sout[19] + sout[27] + sout[43] + sout[51] + sout[59]
   + xorout[10] - xorout[11] >= -2
 R597: sout[5] + sout[21] + sout[29] + sout[45] + sout[53] - sout[61]
   + xorout[13] >= 0
 R598: sout[5] + sout[21] + sout[29] + sout[45] - sout[53] + sout[61]
   + xorout[13] >= 0
 R599: sout[5] + sout[21] + sout[29] - sout[45] + sout[53] + sout[61]
   + xorout[13] >= 0
 R600: sout[5] + sout[21] - sout[29] + sout[45] + sout[53] + sout[61]
   + xorout[13] >= 0
 R601: sout[5] - sout[21] + sout[29] + sout[45] + sout[53] + sout[61]
   + xorout[13] >= 0
 R602: - sout[5] + sout[21] + sout[29] + sout[45] + sout[53] + sout[61]
   + xorout[13] >= 0
 R603: sout[5] + sout[21] + sout[29] + sout[45] + sout[53] + sout[61]
   - xorout[13] >= 0
 R604: sout[5] + sout[21] + sout[29] - sout[45] - sout[53] - sout[61]
   + xorout[13] >= -2
 R605: sout[5] + sout[21] - sout[29] + sout[45] - sout[53] - sout[61]
   + xorout[13] >= -2
 R606: sout[5] - sout[21] + sout[29] + sout[45] - sout[53] - sout[61]
   + xorout[13] >= -2
 R607: - sout[5] + sout[21] + sout[29] + sout[45] - sout[53] - sout[61]
   + xorout[13] >= -2
 R608: sout[5] + sout[21] - sout[29] - sout[45] + sout[53] - sout[61]
   + xorout[13] >= -2
 R609: sout[5] - sout[21] + sout[29] - sout[45] + sout[53] - sout[61]
   + xorout[13] >= -2
 R610: - sout[5] + sout[21] + sout[29] - sout[45] + sout[53] - sout[61]
   + xorout[13] >= -2
 R611: sout[5] - sout[21] - sout[29] + sout[45] + sout[53] - sout[61]
   + xorout[13] >= -2
 R612: - sout[5] + sout[21] - sout[29] + sout[45] + sout[53] - sout[61]
   + xorout[13] >= -2
 R613: - sout[5] - sout[21] + sout[29] + sout[45] + sout[53] - sout[61]
   + xorout[13] >= -2
 R614: sout[5] + sout[21] - sout[29] - sout[45] - sout[53] + sout[61]
   + xorout[13] >= -2
 R615: sout[5] - sout[21] + sout[29] - sout[45] - sout[53] + sout[61]
   + xorout[13] >= -2
 R616: - sout[5] + sout[21] + sout[29] - sout[45] - sout[53] + sout[61]
   + xorout[13] >= -2
 R617: sout[5] - sout[21] - sout[29] + sout[45] - sout[53] + sout[61]
   + xorout[13] >= -2
 R618: - sout[5] + sout[21] - sout[29] + sout[45] - sout[53] + sout[61]
   + xorout[13] >= -2
 R619: - sout[5] - sout[21] + sout[29] + sout[45] - sout[53] + sout[61]
   + xorout[13] >= -2
 R620: sout[5] - sout[21] - sout[29] - sout[45] + sout[53] + sout[61]
   + xorout[13] >= -2
 R621: - sout[5] + sout[21] - sout[29] - sout[45] + sout[53] + sout[61]
   + xorout[13] >= -2
 R622: - sout[5] - sout[21] + sout[29] - sout[45] + sout[53] + sout[61]
   + xorout[13] >= -2
 R623: - sout[5] - sout[21] - sout[29] + sout[45] + sout[53] + sout[61]
   + xorout[13] >= -2
 R624: sout[5] - sout[21] - sout[29] - sout[45] - sout[53] - sout[61]
   + xorout[13] >= -4
 R625: - sout[5] + sout[21] - sout[29] - sout[45] - sout[53] - sout[61]
   + xorout[13] >= -4
 R626: - sout[5] - sout[21] + sout[29] - sout[45] - sout[53] - sout[61]
   + xorout[13] >= -4
 R627: - sout[5] - sout[21] - sout[29] + sout[45] - sout[53] - sout[61]
   + xorout[13] >= -4
 R628: - sout[5] - sout[21] - sout[29] - sout[45] + sout[53] - sout[61]
   + xorout[13] >= -4
 R629: - sout[5] - sout[21] - sout[29] - sout[45] - sout[53] + sout[61]
   + xorout[13] >= -4
 R630: - sout[60] + xorout[12] >= 0
 R631: - sout[52] + xorout[12] >= 0
 R632: - sout[44] + xorout[12] >= 0
 R633: - sout[28] + xorout[12] >= 0
 R634: - sout[20] + xorout[12] >= 0
 R635: - sout[4] + xorout[12] >= 0
 R636: - sout[4] + sout[5] >= 0
 R637: - sout[20] + sout[21] >= 0
 R638: - sout[28] + sout[29] >= 0
 R639: - sout[44] + sout[45] >= 0
 R640: - sout[52] + sout[53] >= 0
 R641: - sout[60] + sout[61] >= 0
 R642: sout[4] + sout[20] + sout[28] + sout[44] + sout[52] + sout[60]
   - xorout[12] >= 0
 R643: - xorout[12] + xorout[13] >= 0
 R644: - sout[5] - sout[21] - sout[29] - sout[45] - sout[53] - sout[61]
   + xorout[12] - xorout[13] >= -6
 R645: sout[5] + sout[21] - sout[29] - sout[45] - sout[53] - sout[61]
   + xorout[12] - xorout[13] >= -4
 R646: sout[5] - sout[21] + sout[29] - sout[45] - sout[53] - sout[61]
   + xorout[12] - xorout[13] >= -4
 R647: - sout[5] + sout[21] + sout[29] - sout[45] - sout[53] - sout[61]
   + xorout[12] - xorout[13] >= -4
 R648: sout[5] - sout[21] - sout[29] + sout[45] - sout[53] - sout[61]
   + xorout[12] - xorout[13] >= -4
 R649: - sout[5] + sout[21] - sout[29] + sout[45] - sout[53] - sout[61]
   + xorout[12] - xorout[13] >= -4
 R650: - sout[5] - sout[21] + sout[29] + sout[45] - sout[53] - sout[61]
   + xorout[12] - xorout[13] >= -4
 R651: sout[5] - sout[21] - sout[29] - sout[45] + sout[53] - sout[61]
   + xorout[12] - xorout[13] >= -4
 R652: - sout[5] + sout[21] - sout[29] - sout[45] + sout[53] - sout[61]
   + xorout[12] - xorout[13] >= -4
 R653: - sout[5] - sout[21] + sout[29] - sout[45] + sout[53] - sout[61]
   + xorout[12] - xorout[13] >= -4
 R654: - sout[5] - sout[21] - sout[29] + sout[45] + sout[53] - sout[61]
   + xorout[12] - xorout[13] >= -4
 R655: sout[5] - sout[21] - sout[29] - sout[45] - sout[53] + sout[61]
   + xorout[12] - xorout[13] >= -4
 R656: - sout[5] + sout[21] - sout[29] - sout[45] - sout[53] + sout[61]
   + xorout[12] - xorout[13] >= -4
 R657: - sout[5] - sout[21] + sout[29] - sout[45] - sout[53] + sout[61]
   + xorout[12] - xorout[13] >= -4
 R658: - sout[5] - sout[21] - sout[29] + sout[45] - sout[53] + sout[61]
   + xorout[12] - xorout[13] >= -4
 R659: - sout[5] - sout[21] - sout[29] - sout[45] + sout[53] + sout[61]
   + xorout[12] - xorout[13] >= -4
 R660: sout[5] + sout[21] + sout[29] + sout[45] - sout[53] - sout[61]
   + xorout[12] - xorout[13] >= -2
 R661: sout[5] + sout[21] + sout[29] - sout[45] + sout[53] - sout[61]
   + xorout[12] - xorout[13] >= -2
 R662: sout[5] + sout[21] - sout[29] + sout[45] + sout[53] - sout[61]
   + xorout[12] - xorout[13] >= -2
 R663: sout[5] - sout[21] + sout[29] + sout[45] + sout[53] - sout[61]
   + xorout[12] - xorout[13] >= -2
 R664: - sout[5] + sout[21] + sout[29] + sout[45] + sout[53] - sout[61]
   + xorout[12] - xorout[13] >= -2
 R665: sout[5] + sout[21] + sout[29] - sout[45] - sout[53] + sout[61]
   + xorout[12] - xorout[13] >= -2
 R666: sout[5] + sout[21] - sout[29] + sout[45] - sout[53] + sout[61]
   + xorout[12] - xorout[13] >= -2
 R667: sout[5] - sout[21] + sout[29] + sout[45] - sout[53] + sout[61]
   + xorout[12] - xorout[13] >= -2
 R668: - sout[5] + sout[21] + sout[29] + sout[45] - sout[53] + sout[61]
   + xorout[12] - xorout[13] >= -2
 R669: sout[5] + sout[21] - sout[29] - sout[45] + sout[53] + sout[61]
   + xorout[12] - xorout[13] >= -2
 R670: sout[5] - sout[21] + sout[29] - sout[45] + sout[53] + sout[61]
   + xorout[12] - xorout[13] >= -2
 R671: - sout[5] + sout[21] + sout[29] - sout[45] + sout[53] + sout[61]
   + xorout[12] - xorout[13] >= -2
 R672: sout[5] - sout[21] - sout[29] + sout[45] + sout[53] + sout[61]
   + xorout[12] - xorout[13] >= -2
 R673: - sout[5] + sout[21] - sout[29] + sout[45] + sout[53] + sout[61]
   + xorout[12] - xorout[13] >= -2
 R674: - sout[5] - sout[21] + sout[29] + sout[45] + sout[53] + sout[61]
   + xorout[12] - xorout[13] >= -2
 R675: sout[7] + sout[23] + sout[31] + sout[47] + sout[55] - sout[63]
   + xorout[15] >= 0
 R676: sout[7] + sout[23] + sout[31] + sout[47] - sout[55] + sout[63]
   + xorout[15] >= 0
 R677: sout[7] + sout[23] + sout[31] - sout[47] + sout[55] + sout[63]
   + xorout[15] >= 0
 R678: sout[7] + sout[23] - sout[31] + sout[47] + sout[55] + sout[63]
   + xorout[15] >= 0
 R679: sout[7] - sout[23] + sout[31] + sout[47] + sout[55] + sout[63]
   + xorout[15] >= 0
 R680: - sout[7] + sout[23] + sout[31] + sout[47] + sout[55] + sout[63]
   + xorout[15] >= 0
 R681: sout[7] + sout[23] + sout[31] + sout[47] + sout[55] + sout[63]
   - xorout[15] >= 0
 R682: sout[7] + sout[23] + sout[31] - sout[47] - sout[55] - sout[63]
   + xorout[15] >= -2
 R683: sout[7] + sout[23] - sout[31] + sout[47] - sout[55] - sout[63]
   + xorout[15] >= -2
 R684: sout[7] - sout[23] + sout[31] + sout[47] - sout[55] - sout[63]
   + xorout[15] >= -2
 R685: - sout[7] + sout[23] + sout[31] + sout[47] - sout[55] - sout[63]
   + xorout[15] >= -2
 R686: sout[7] + sout[23] - sout[31] - sout[47] + sout[55] - sout[63]
   + xorout[15] >= -2
 R687: sout[7] - sout[23] + sout[31] - sout[47] + sout[55] - sout[63]
   + xorout[15] >= -2
 R688: - sout[7] + sout[23] + sout[31] - sout[47] + sout[55] - sout[63]
   + xorout[15] >= -2
 R689: sout[7] - sout[23] - sout[31] + sout[47] + sout[55] - sout[63]
   + xorout[15] >= -2
 R690: - sout[7] + sout[23] - sout[31] + sout[47] + sout[55] - sout[63]
   + xorout[15] >= -2
 R691: - sout[7] - sout[23] + sout[31] + sout[47] + sout[55] - sout[63]
   + xorout[15] >= -2
 R692: sout[7] + sout[23] - sout[31] - sout[47] - sout[55] + sout[63]
   + xorout[15] >= -2
 R693: sout[7] - sout[23] + sout[31] - sout[47] - sout[55] + sout[63]
   + xorout[15] >= -2
 R694: - sout[7] + sout[23] + sout[31] - sout[47] - sout[55] + sout[63]
   + xorout[15] >= -2
 R695: sout[7] - sout[23] - sout[31] + sout[47] - sout[55] + sout[63]
   + xorout[15] >= -2
 R696: - sout[7] + sout[23] - sout[31] + sout[47] - sout[55] + sout[63]
   + xorout[15] >= -2
 R697: - sout[7] - sout[23] + sout[31] + sout[47] - sout[55] + sout[63]
   + xorout[15] >= -2
 R698: sout[7] - sout[23] - sout[31] - sout[47] + sout[55] + sout[63]
   + xorout[15] >= -2
 R699: - sout[7] + sout[23] - sout[31] - sout[47] + sout[55] + sout[63]
   + xorout[15] >= -2
 R700: - sout[7] - sout[23] + sout[31] - sout[47] + sout[55] + sout[63]
   + xorout[15] >= -2
 R701: - sout[7] - sout[23] - sout[31] + sout[47] + sout[55] + sout[63]
   + xorout[15] >= -2
 R702: sout[7] - sout[23] - sout[31] - sout[47] - sout[55] - sout[63]
   + xorout[15] >= -4
 R703: - sout[7] + sout[23] - sout[31] - sout[47] - sout[55] - sout[63]
   + xorout[15] >= -4
 R704: - sout[7] - sout[23] + sout[31] - sout[47] - sout[55] - sout[63]
   + xorout[15] >= -4
 R705: - sout[7] - sout[23] - sout[31] + sout[47] - sout[55] - sout[63]
   + xorout[15] >= -4
 R706: - sout[7] - sout[23] - sout[31] - sout[47] + sout[55] - sout[63]
   + xorout[15] >= -4
 R707: - sout[7] - sout[23] - sout[31] - sout[47] - sout[55] + sout[63]
   + xorout[15] >= -4
 R708: - sout[62] + xorout[14] >= 0
 R709: - sout[54] + xorout[14] >= 0
 R710: - sout[46] + xorout[14] >= 0
 R711: - sout[30] + xorout[14] >= 0
 R712: - sout[22] + xorout[14] >= 0
 R713: - sout[6] + xorout[14] >= 0
 R714: - sout[6] + sout[7] >= 0
 R715: - sout[22] + sout[23] >= 0
 R716: - sout[30] + sout[31] >= 0
 R717: - sout[46] + sout[47] >= 0
 R718: - sout[54] + sout[55] >= 0
 R719: - sout[62] + sout[63] >= 0
 R720: sout[6] + sout[22] + sout[30] + sout[46] + sout[54] + sout[62]
   - xorout[14] >= 0
 R721: - xorout[14] + xorout[15] >= 0
 R722: - sout[7] - sout[23] - sout[31] - sout[47] - sout[55] - sout[63]
   + xorout[14] - xorout[15] >= -6
 R723: sout[7] + sout[23] - sout[31] - sout[47] - sout[55] - sout[63]
   + xorout[14] - xorout[15] >= -4
 R724: sout[7] - sout[23] + sout[31] - sout[47] - sout[55] - sout[63]
   + xorout[14] - xorout[15] >= -4
 R725: - sout[7] + sout[23] + sout[31] - sout[47] - sout[55] - sout[63]
   + xorout[14] - xorout[15] >= -4
 R726: sout[7] - sout[23] - sout[31] + sout[47] - sout[55] - sout[63]
   + xorout[14] - xorout[15] >= -4
 R727: - sout[7] + sout[23] - sout[31] + sout[47] - sout[55] - sout[63]
   + xorout[14] - xorout[15] >= -4
 R728: - sout[7] - sout[23] + sout[31] + sout[47] - sout[55] - sout[63]
   + xorout[14] - xorout[15] >= -4
 R729: sout[7] - sout[23] - sout[31] - sout[47] + sout[55] - sout[63]
   + xorout[14] - xorout[15] >= -4
 R730: - sout[7] + sout[23] - sout[31] - sout[47] + sout[55] - sout[63]
   + xorout[14] - xorout[15] >= -4
 R731: - sout[7] - sout[23] + sout[31] - sout[47] + sout[55] - sout[63]
   + xorout[14] - xorout[15] >= -4
 R732: - sout[7] - sout[23] - sout[31] + sout[47] + sout[55] - sout[63]
   + xorout[14] - xorout[15] >= -4
 R733: sout[7] - sout[23] - sout[31] - sout[47] - sout[55] + sout[63]
   + xorout[14] - xorout[15] >= -4
 R734: - sout[7] + sout[23] - sout[31] - sout[47] - sout[55] + sout[63]
   + xorout[14] - xorout[15] >= -4
 R735: - sout[7] - sout[23] + sout[31] - sout[47] - sout[55] + sout[63]
   + xorout[14] - xorout[15] >= -4
 R736: - sout[7] - sout[23] - sout[31] + sout[47] - sout[55] + sout[63]
   + xorout[14] - xorout[15] >= -4
 R737: - sout[7] - sout[23] - sout[31] - sout[47] + sout[55] + sout[63]
   + xorout[14] - xorout[15] >= -4
 R738: sout[7] + sout[23] + sout[31] + sout[47] - sout[55] - sout[63]
   + xorout[14] - xorout[15] >= -2
 R739: sout[7] + sout[23] + sout[31] - sout[47] + sout[55] - sout[63]
   + xorout[14] - xorout[15] >= -2
 R740: sout[7] + sout[23] - sout[31] + sout[47] + sout[55] - sout[63]
   + xorout[14] - xorout[15] >= -2
 R741: sout[7] - sout[23] + sout[31] + sout[47] + sout[55] - sout[63]
   + xorout[14] - xorout[15] >= -2
 R742: - sout[7] + sout[23] + sout[31] + sout[47] + sout[55] - sout[63]
   + xorout[14] - xorout[15] >= -2
 R743: sout[7] + sout[23] + sout[31] - sout[47] - sout[55] + sout[63]
   + xorout[14] - xorout[15] >= -2
 R744: sout[7] + sout[23] - sout[31] + sout[47] - sout[55] + sout[63]
   + xorout[14] - xorout[15] >= -2
 R745: sout[7] - sout[23] + sout[31] + sout[47] - sout[55] + sout[63]
   + xorout[14] - xorout[15] >= -2
 R746: - sout[7] + sout[23] + sout[31] + sout[47] - sout[55] + sout[63]
   + xorout[14] - xorout[15] >= -2
 R747: sout[7] + sout[23] - sout[31] - sout[47] + sout[55] + sout[63]
   + xorout[14] - xorout[15] >= -2
 R748: sout[7] - sout[23] + sout[31] - sout[47] + sout[55] + sout[63]
   + xorout[14] - xorout[15] >= -2
 R749: - sout[7] + sout[23] + sout[31] - sout[47] + sout[55] + sout[63]
   + xorout[14] - xorout[15] >= -2
 R750: sout[7] - sout[23] - sout[31] + sout[47] + sout[55] + sout[63]
   + xorout[14] - xorout[15] >= -2
 R751: - sout[7] + sout[23] - sout[31] + sout[47] + sout[55] + sout[63]
   + xorout[14] - xorout[15] >= -2
 R752: - sout[7] - sout[23] + sout[31] + sout[47] + sout[55] + sout[63]
   + xorout[14] - xorout[15] >= -2
 R753: sout[1] + sout[9] + sout[25] + sout[33] + sout[49] - sout[57]
   + xorout[17] >= 0
 R754: sout[1] + sout[9] + sout[25] + sout[33] - sout[49] + sout[57]
   + xorout[17] >= 0
 R755: sout[1] + sout[9] + sout[25] - sout[33] + sout[49] + sout[57]
   + xorout[17] >= 0
 R756: sout[1] + sout[9] - sout[25] + sout[33] + sout[49] + sout[57]
   + xorout[17] >= 0
 R757: sout[1] - sout[9] + sout[25] + sout[33] + sout[49] + sout[57]
   + xorout[17] >= 0
 R758: - sout[1] + sout[9] + sout[25] + sout[33] + sout[49] + sout[57]
   + xorout[17] >= 0
 R759: sout[1] + sout[9] + sout[25] + sout[33] + sout[49] + sout[57]
   - xorout[17] >= 0
 R760: sout[1] + sout[9] + sout[25] - sout[33] - sout[49] - sout[57]
   + xorout[17] >= -2
 R761: sout[1] + sout[9] - sout[25] + sout[33] - sout[49] - sout[57]
   + xorout[17] >= -2
 R762: sout[1] - sout[9] + sout[25] + sout[33] - sout[49] - sout[57]
   + xorout[17] >= -2
 R763: - sout[1] + sout[9] + sout[25] + sout[33] - sout[49] - sout[57]
   + xorout[17] >= -2
 R764: sout[1] + sout[9] - sout[25] - sout[33] + sout[49] - sout[57]
   + xorout[17] >= -2
 R765: sout[1] - sout[9] + sout[25] - sout[33] + sout[49] - sout[57]
   + xorout[17] >= -2
 R766: - sout[1] + sout[9] + sout[25] - sout[33] + sout[49] - sout[57]
   + xorout[17] >= -2
 R767: sout[1] - sout[9] - sout[25] + sout[33] + sout[49] - sout[57]
   + xorout[17] >= -2
 R768: - sout[1] + sout[9] - sout[25] + sout[33] + sout[49] - sout[57]
   + xorout[17] >= -2
 R769: - sout[1] - sout[9] + sout[25] + sout[33] + sout[49] - sout[57]
   + xorout[17] >= -2
 R770: sout[1] + sout[9] - sout[25] - sout[33] - sout[49] + sout[57]
   + xorout[17] >= -2
 R771: sout[1] - sout[9] + sout[25] - sout[33] - sout[49] + sout[57]
   + xorout[17] >= -2
 R772: - sout[1] + sout[9] + sout[25] - sout[33] - sout[49] + sout[57]
   + xorout[17] >= -2
 R773: sout[1] - sout[9] - sout[25] + sout[33] - sout[49] + sout[57]
   + xorout[17] >= -2
 R774: - sout[1] + sout[9] - sout[25] + sout[33] - sout[49] + sout[57]
   + xorout[17] >= -2
 R775: - sout[1] - sout[9] + sout[25] + sout[33] - sout[49] + sout[57]
   + xorout[17] >= -2
 R776: sout[1] - sout[9] - sout[25] - sout[33] + sout[49] + sout[57]
   + xorout[17] >= -2
 R777: - sout[1] + sout[9] - sout[25] - sout[33] + sout[49] + sout[57]
   + xorout[17] >= -2
 R778: - sout[1] - sout[9] + sout[25] - sout[33] + sout[49] + sout[57]
   + xorout[17] >= -2
 R779: - sout[1] - sout[9] - sout[25] + sout[33] + sout[49] + sout[57]
   + xorout[17] >= -2
 R780: sout[1] - sout[9] - sout[25] - sout[33] - sout[49] - sout[57]
   + xorout[17] >= -4
 R781: - sout[1] + sout[9] - sout[25] - sout[33] - sout[49] - sout[57]
   + xorout[17] >= -4
 R782: - sout[1] - sout[9] + sout[25] - sout[33] - sout[49] - sout[57]
   + xorout[17] >= -4
 R783: - sout[1] - sout[9] - sout[25] + sout[33] - sout[49] - sout[57]
   + xorout[17] >= -4
 R784: - sout[1] - sout[9] - sout[25] - sout[33] + sout[49] - sout[57]
   + xorout[17] >= -4
 R785: - sout[1] - sout[9] - sout[25] - sout[33] - sout[49] + sout[57]
   + xorout[17] >= -4
 R786: - sout[56] + xorout[16] >= 0
 R787: - sout[48] + xorout[16] >= 0
 R788: - sout[32] + xorout[16] >= 0
 R789: - sout[24] + xorout[16] >= 0
 R790: - sout[8] + xorout[16] >= 0
 R791: - sout[0] + xorout[16] >= 0
 R792: - sout[0] + sout[1] >= 0
 R793: - sout[8] + sout[9] >= 0
 R794: - sout[24] + sout[25] >= 0
 R795: - sout[32] + sout[33] >= 0
 R796: - sout[48] + sout[49] >= 0
 R797: - sout[56] + sout[57] >= 0
 R798: sout[0] + sout[8] + sout[24] + sout[32] + sout[48] + sout[56]
   - xorout[16] >= 0
 R799: - xorout[16] + xorout[17] >= 0
 R800: - sout[1] - sout[9] - sout[25] - sout[33] - sout[49] - sout[57]
   + xorout[16] - xorout[17] >= -6
 R801: sout[1] + sout[9] - sout[25] - sout[33] - sout[49] - sout[57]
   + xorout[16] - xorout[17] >= -4
 R802: sout[1] - sout[9] + sout[25] - sout[33] - sout[49] - sout[57]
   + xorout[16] - xorout[17] >= -4
 R803: - sout[1] + sout[9] + sout[25] - sout[33] - sout[49] - sout[57]
   + xorout[16] - xorout[17] >= -4
 R804: sout[1] - sout[9] - sout[25] + sout[33] - sout[49] - sout[57]
   + xorout[16] - xorout[17] >= -4
 R805: - sout[1] + sout[9] - sout[25] + sout[33] - sout[49] - sout[57]
   + xorout[16] - xorout[17] >= -4
 R806: - sout[1] - sout[9] + sout[25] + sout[33] - sout[49] - sout[57]
   + xorout[16] - xorout[17] >= -4
 R807: sout[1] - sout[9] - sout[25] - sout[33] + sout[49] - sout[57]
   + xorout[16] - xorout[17] >= -4
 R808: - sout[1] + sout[9] - sout[25] - sout[33] + sout[49] - sout[57]
   + xorout[16] - xorout[17] >= -4
 R809: - sout[1] - sout[9] + sout[25] - sout[33] + sout[49] - sout[57]
   + xorout[16] - xorout[17] >= -4
 R810: - sout[1] - sout[9] - sout[25] + sout[33] + sout[49] - sout[57]
   + xorout[16] - xorout[17] >= -4
 R811: sout[1] - sout[9] - sout[25] - sout[33] - sout[49] + sout[57]
   + xorout[16] - xorout[17] >= -4
 R812: - sout[1] + sout[9] - sout[25] - sout[33] - sout[49] + sout[57]
   + xorout[16] - xorout[17] >= -4
 R813: - sout[1] - sout[9] + sout[25] - sout[33] - sout[49] + sout[57]
   + xorout[16] - xorout[17] >= -4
 R814: - sout[1] - sout[9] - sout[25] + sout[33] - sout[49] + sout[57]
   + xorout[16] - xorout[17] >= -4
 R815: - sout[1] - sout[9] - sout[25] - sout[33] + sout[49] + sout[57]
   + xorout[16] - xorout[17] >= -4
 R816: sout[1] + sout[9] + sout[25] + sout[33] - sout[49] - sout[57]
   + xorout[16] - xorout[17] >= -2
 R817: sout[1] + sout[9] + sout[25] - sout[33] + sout[49] - sout[57]
   + xorout[16] - xorout[17] >= -2
 R818: sout[1] + sout[9] - sout[25] + sout[33] + sout[49] - sout[57]
   + xorout[16] - xorout[17] >= -2
 R819: sout[1] - sout[9] + sout[25] + sout[33] + sout[49] - sout[57]
   + xorout[16] - xorout[17] >= -2
 R820: - sout[1] + sout[9] + sout[25] + sout[33] + sout[49] - sout[57]
   + xorout[16] - xorout[17] >= -2
 R821: sout[1] + sout[9] + sout[25] - sout[33] - sout[49] + sout[57]
   + xorout[16] - xorout[17] >= -2
 R822: sout[1] + sout[9] - sout[25] + sout[33] - sout[49] + sout[57]
   + xorout[16] - xorout[17] >= -2
 R823: sout[1] - sout[9] + sout[25] + sout[33] - sout[49] + sout[57]
   + xorout[16] - xorout[17] >= -2
 R824: - sout[1] + sout[9] + sout[25] + sout[33] - sout[49] + sout[57]
   + xorout[16] - xorout[17] >= -2
 R825: sout[1] + sout[9] - sout[25] - sout[33] + sout[49] + sout[57]
   + xorout[16] - xorout[17] >= -2
 R826: sout[1] - sout[9] + sout[25] - sout[33] + sout[49] + sout[57]
   + xorout[16] - xorout[17] >= -2
 R827: - sout[1] + sout[9] + sout[25] - sout[33] + sout[49] + sout[57]
   + xorout[16] - xorout[17] >= -2
 R828: sout[1] - sout[9] - sout[25] + sout[33] + sout[49] + sout[57]
   + xorout[16] - xorout[17] >= -2
 R829: - sout[1] + sout[9] - sout[25] + sout[33] + sout[49] + sout[57]
   + xorout[16] - xorout[17] >= -2
 R830: - sout[1] - sout[9] + sout[25] + sout[33] + sout[49] + sout[57]
   + xorout[16] - xorout[17] >= -2
 R831: sout[3] + sout[11] + sout[27] + sout[35] + sout[51] - sout[59]
   + xorout[19] >= 0
 R832: sout[3] + sout[11] + sout[27] + sout[35] - sout[51] + sout[59]
   + xorout[19] >= 0
 R833: sout[3] + sout[11] + sout[27] - sout[35] + sout[51] + sout[59]
   + xorout[19] >= 0
 R834: sout[3] + sout[11] - sout[27] + sout[35] + sout[51] + sout[59]
   + xorout[19] >= 0
 R835: sout[3] - sout[11] + sout[27] + sout[35] + sout[51] + sout[59]
   + xorout[19] >= 0
 R836: - sout[3] + sout[11] + sout[27] + sout[35] + sout[51] + sout[59]
   + xorout[19] >= 0
 R837: sout[3] + sout[11] + sout[27] + sout[35] + sout[51] + sout[59]
   - xorout[19] >= 0
 R838: sout[3] + sout[11] + sout[27] - sout[35] - sout[51] - sout[59]
   + xorout[19] >= -2
 R839: sout[3] + sout[11] - sout[27] + sout[35] - sout[51] - sout[59]
   + xorout[19] >= -2
 R840: sout[3] - sout[11] + sout[27] + sout[35] - sout[51] - sout[59]
   + xorout[19] >= -2
 R841: - sout[3] + sout[11] + sout[27] + sout[35] - sout[51] - sout[59]
   + xorout[19] >= -2
 R842: sout[3] + sout[11] - sout[27] - sout[35] + sout[51] - sout[59]
   + xorout[19] >= -2
 R843: sout[3] - sout[11] + sout[27] - sout[35] + sout[51] - sout[59]
   + xorout[19] >= -2
 R844: - sout[3] + sout[11] + sout[27] - sout[35] + sout[51] - sout[59]
   + xorout[19] >= -2
 R845: sout[3] - sout[11] - sout[27] + sout[35] + sout[51] - sout[59]
   + xorout[19] >= -2
 R846: - sout[3] + sout[11] - sout[27] + sout[35] + sout[51] - sout[59]
   + xorout[19] >= -2
 R847: - sout[3] - sout[11] + sout[27] + sout[35] + sout[51] - sout[59]
   + xorout[19] >= -2
 R848: sout[3] + sout[11] - sout[27] - sout[35] - sout[51] + sout[59]
   + xorout[19] >= -2
 R849: sout[3] - sout[11] + sout[27] - sout[35] - sout[51] + sout[59]
   + xorout[19] >= -2
 R850: - sout[3] + sout[11] + sout[27] - sout[35] - sout[51] + sout[59]
   + xorout[19] >= -2
 R851: sout[3] - sout[11] - sout[27] + sout[35] - sout[51] + sout[59]
   + xorout[19] >= -2
 R852: - sout[3] + sout[11] - sout[27] + sout[35] - sout[51] + sout[59]
   + xorout[19] >= -2
 R853: - sout[3] - sout[11] + sout[27] + sout[35] - sout[51] + sout[59]
   + xorout[19] >= -2
 R854: sout[3] - sout[11] - sout[27] - sout[35] + sout[51] + sout[59]
   + xorout[19] >= -2
 R855: - sout[3] + sout[11] - sout[27] - sout[35] + sout[51] + sout[59]
   + xorout[19] >= -2
 R856: - sout[3] - sout[11] + sout[27] - sout[35] + sout[51] + sout[59]
   + xorout[19] >= -2
 R857: - sout[3] - sout[11] - sout[27] + sout[35] + sout[51] + sout[59]
   + xorout[19] >= -2
 R858: sout[3] - sout[11] - sout[27] - sout[35] - sout[51] - sout[59]
   + xorout[19] >= -4
 R859: - sout[3] + sout[11] - sout[27] - sout[35] - sout[51] - sout[59]
   + xorout[19] >= -4
 R860: - sout[3] - sout[11] + sout[27] - sout[35] - sout[51] - sout[59]
   + xorout[19] >= -4
 R861: - sout[3] - sout[11] - sout[27] + sout[35] - sout[51] - sout[59]
   + xorout[19] >= -4
 R862: - sout[3] - sout[11] - sout[27] - sout[35] + sout[51] - sout[59]
   + xorout[19] >= -4
 R863: - sout[3] - sout[11] - sout[27] - sout[35] - sout[51] + sout[59]
   + xorout[19] >= -4
 R864: - sout[58] + xorout[18] >= 0
 R865: - sout[50] + xorout[18] >= 0
 R866: - sout[34] + xorout[18] >= 0
 R867: - sout[26] + xorout[18] >= 0
 R868: - sout[10] + xorout[18] >= 0
 R869: - sout[2] + xorout[18] >= 0
 R870: - sout[2] + sout[3] >= 0
 R871: - sout[10] + sout[11] >= 0
 R872: - sout[26] + sout[27] >= 0
 R873: - sout[34] + sout[35] >= 0
 R874: - sout[50] + sout[51] >= 0
 R875: - sout[58] + sout[59] >= 0
 R876: sout[2] + sout[10] + sout[26] + sout[34] + sout[50] + sout[58]
   - xorout[18] >= 0
 R877: - xorout[18] + xorout[19] >= 0
 R878: - sout[3] - sout[11] - sout[27] - sout[35] - sout[51] - sout[59]
   + xorout[18] - xorout[19] >= -6
 R879: sout[3] + sout[11] - sout[27] - sout[35] - sout[51] - sout[59]
   + xorout[18] - xorout[19] >= -4
 R880: sout[3] - sout[11] + sout[27] - sout[35] - sout[51] - sout[59]
   + xorout[18] - xorout[19] >= -4
 R881: - sout[3] + sout[11] + sout[27] - sout[35] - sout[51] - sout[59]
   + xorout[18] - xorout[19] >= -4
 R882: sout[3] - sout[11] - sout[27] + sout[35] - sout[51] - sout[59]
   + xorout[18] - xorout[19] >= -4
 R883: - sout[3] + sout[11] - sout[27] + sout[35] - sout[51] - sout[59]
   + xorout[18] - xorout[19] >= -4
 R884: - sout[3] - sout[11] + sout[27] + sout[35] - sout[51] - sout[59]
   + xorout[18] - xorout[19] >= -4
 R885: sout[3] - sout[11] - sout[27] - sout[35] + sout[51] - sout[59]
   + xorout[18] - xorout[19] >= -4
 R886: - sout[3] + sout[11] - sout[27] - sout[35] + sout[51] - sout[59]
   + xorout[18] - xorout[19] >= -4
 R887: - sout[3] - sout[11] + sout[27] - sout[35] + sout[51] - sout[59]
   + xorout[18] - xorout[19] >= -4
 R888: - sout[3] - sout[11] - sout[27] + sout[35] + sout[51] - sout[59]
   + xorout[18] - xorout[19] >= -4
 R889: sout[3] - sout[11] - sout[27] - sout[35] - sout[51] + sout[59]
   + xorout[18] - xorout[19] >= -4
 R890: - sout[3] + sout[11] - sout[27] - sout[35] - sout[51] + sout[59]
   + xorout[18] - xorout[19] >= -4
 R891: - sout[3] - sout[11] + sout[27] - sout[35] - sout[51] + sout[59]
   + xorout[18] - xorout[19] >= -4
 R892: - sout[3] - sout[11] - sout[27] + sout[35] - sout[51] + sout[59]
   + xorout[18] - xorout[19] >= -4
 R893: - sout[3] - sout[11] - sout[27] - sout[35] + sout[51] + sout[59]
   + xorout[18] - xorout[19] >= -4
 R894: sout[3] + sout[11] + sout[27] + sout[35] - sout[51] - sout[59]
   + xorout[18] - xorout[19] >= -2
 R895: sout[3] + sout[11] + sout[27] - sout[35] + sout[51] - sout[59]
   + xorout[18] - xorout[19] >= -2
 R896: sout[3] + sout[11] - sout[27] + sout[35] + sout[51] - sout[59]
   + xorout[18] - xorout[19] >= -2
 R897: sout[3] - sout[11] + sout[27] + sout[35] + sout[51] - sout[59]
   + xorout[18] - xorout[19] >= -2
 R898: - sout[3] + sout[11] + sout[27] + sout[35] + sout[51] - sout[59]
   + xorout[18] - xorout[19] >= -2
 R899: sout[3] + sout[11] + sout[27] - sout[35] - sout[51] + sout[59]
   + xorout[18] - xorout[19] >= -2
 R900: sout[3] + sout[11] - sout[27] + sout[35] - sout[51] + sout[59]
   + xorout[18] - xorout[19] >= -2
 R901: sout[3] - sout[11] + sout[27] + sout[35] - sout[51] + sout[59]
   + xorout[18] - xorout[19] >= -2
 R902: - sout[3] + sout[11] + sout[27] + sout[35] - sout[51] + sout[59]
   + xorout[18] - xorout[19] >= -2
 R903: sout[3] + sout[11] - sout[27] - sout[35] + sout[51] + sout[59]
   + xorout[18] - xorout[19] >= -2
 R904: sout[3] - sout[11] + sout[27] - sout[35] + sout[51] + sout[59]
   + xorout[18] - xorout[19] >= -2
 R905: - sout[3] + sout[11] + sout[27] - sout[35] + sout[51] + sout[59]
   + xorout[18] - xorout[19] >= -2
 R906: sout[3] - sout[11] - sout[27] + sout[35] + sout[51] + sout[59]
   + xorout[18] - xorout[19] >= -2
 R907: - sout[3] + sout[11] - sout[27] + sout[35] + sout[51] + sout[59]
   + xorout[18] - xorout[19] >= -2
 R908: - sout[3] - sout[11] + sout[27] + sout[35] + sout[51] + sout[59]
   + xorout[18] - xorout[19] >= -2
 R909: sout[5] + sout[13] + sout[29] + sout[37] + sout[53] - sout[61]
   + xorout[21] >= 0
 R910: sout[5] + sout[13] + sout[29] + sout[37] - sout[53] + sout[61]
   + xorout[21] >= 0
 R911: sout[5] + sout[13] + sout[29] - sout[37] + sout[53] + sout[61]
   + xorout[21] >= 0
 R912: sout[5] + sout[13] - sout[29] + sout[37] + sout[53] + sout[61]
   + xorout[21] >= 0
 R913: sout[5] - sout[13] + sout[29] + sout[37] + sout[53] + sout[61]
   + xorout[21] >= 0
 R914: - sout[5] + sout[13] + sout[29] + sout[37] + sout[53] + sout[61]
   + xorout[21] >= 0
 R915: sout[5] + sout[13] + sout[29] + sout[37] + sout[53] + sout[61]
   - xorout[21] >= 0
 R916: sout[5] + sout[13] + sout[29] - sout[37] - sout[53] - sout[61]
   + xorout[21] >= -2
 R917: sout[5] + sout[13] - sout[29] + sout[37] - sout[53] - sout[61]
   + xorout[21] >= -2
 R918: sout[5] - sout[13] + sout[29] + sout[37] - sout[53] - sout[61]
   + xorout[21] >= -2
 R919: - sout[5] + sout[13] + sout[29] + sout[37] - sout[53] - sout[61]
   + xorout[21] >= -2
 R920: sout[5] + sout[13] - sout[29] - sout[37] + sout[53] - sout[61]
   + xorout[21] >= -2
 R921: sout[5] - sout[13] + sout[29] - sout[37] + sout[53] - sout[61]
   + xorout[21] >= -2
 R922: - sout[5] + sout[13] + sout[29] - sout[37] + sout[53] - sout[61]
   + xorout[21] >= -2
 R923: sout[5] - sout[13] - sout[29] + sout[37] + sout[53] - sout[61]
   + xorout[21] >= -2
 R924: - sout[5] + sout[13] - sout[29] + sout[37] + sout[53] - sout[61]
   + xorout[21] >= -2
 R925: - sout[5] - sout[13] + sout[29] + sout[37] + sout[53] - sout[61]
   + xorout[21] >= -2
 R926: sout[5] + sout[13] - sout[29] - sout[37] - sout[53] + sout[61]
   + xorout[21] >= -2
 R927: sout[5] - sout[13] + sout[29] - sout[37] - sout[53] + sout[61]
   + xorout[21] >= -2
 R928: - sout[5] + sout[13] + sout[29] - sout[37] - sout[53] + sout[61]
   + xorout[21] >= -2
 R929: sout[5] - sout[13] - sout[29] + sout[37] - sout[53] + sout[61]
   + xorout[21] >= -2
 R930: - sout[5] + sout[13] - sout[29] + sout[37] - sout[53] + sout[61]
   + xorout[21] >= -2
 R931: - sout[5] - sout[13] + sout[29] + sout[37] - sout[53] + sout[61]
   + xorout[21] >= -2
 R932: sout[5] - sout[13] - sout[29] - sout[37] + sout[53] + sout[61]
   + xorout[21] >= -2
 R933: - sout[5] + sout[13] - sout[29] - sout[37] + sout[53] + sout[61]
   + xorout[21] >= -2
 R934: - sout[5] - sout[13] + sout[29] - sout[37] + sout[53] + sout[61]
   + xorout[21] >= -2
 R935: - sout[5] - sout[13] - sout[29] + sout[37] + sout[53] + sout[61]
   + xorout[21] >= -2
 R936: sout[5] - sout[13] - sout[29] - sout[37] - sout[53] - sout[61]
   + xorout[21] >= -4
 R937: - sout[5] + sout[13] - sout[29] - sout[37] - sout[53] - sout[61]
   + xorout[21] >= -4
 R938: - sout[5] - sout[13] + sout[29] - sout[37] - sout[53] - sout[61]
   + xorout[21] >= -4
 R939: - sout[5] - sout[13] - sout[29] + sout[37] - sout[53] - sout[61]
   + xorout[21] >= -4
 R940: - sout[5] - sout[13] - sout[29] - sout[37] + sout[53] - sout[61]
   + xorout[21] >= -4
 R941: - sout[5] - sout[13] - sout[29] - sout[37] - sout[53] + sout[61]
   + xorout[21] >= -4
 R942: - sout[60] + xorout[20] >= 0
 R943: - sout[52] + xorout[20] >= 0
 R944: - sout[36] + xorout[20] >= 0
 R945: - sout[28] + xorout[20] >= 0
 R946: - sout[12] + xorout[20] >= 0
 R947: - sout[4] + xorout[20] >= 0
 R948: - sout[4] + sout[5] >= 0
 R949: - sout[12] + sout[13] >= 0
 R950: - sout[28] + sout[29] >= 0
 R951: - sout[36] + sout[37] >= 0
 R952: - sout[52] + sout[53] >= 0
 R953: - sout[60] + sout[61] >= 0
 R954: sout[4] + sout[12] + sout[28] + sout[36] + sout[52] + sout[60]
   - xorout[20] >= 0
 R955: - xorout[20] + xorout[21] >= 0
 R956: - sout[5] - sout[13] - sout[29] - sout[37] - sout[53] - sout[61]
   + xorout[20] - xorout[21] >= -6
 R957: sout[5] + sout[13] - sout[29] - sout[37] - sout[53] - sout[61]
   + xorout[20] - xorout[21] >= -4
 R958: sout[5] - sout[13] + sout[29] - sout[37] - sout[53] - sout[61]
   + xorout[20] - xorout[21] >= -4
 R959: - sout[5] + sout[13] + sout[29] - sout[37] - sout[53] - sout[61]
   + xorout[20] - xorout[21] >= -4
 R960: sout[5] - sout[13] - sout[29] + sout[37] - sout[53] - sout[61]
   + xorout[20] - xorout[21] >= -4
 R961: - sout[5] + sout[13] - sout[29] + sout[37] - sout[53] - sout[61]
   + xorout[20] - xorout[21] >= -4
 R962: - sout[5] - sout[13] + sout[29] + sout[37] - sout[53] - sout[61]
   + xorout[20] - xorout[21] >= -4
 R963: sout[5] - sout[13] - sout[29] - sout[37] + sout[53] - sout[61]
   + xorout[20] - xorout[21] >= -4
 R964: - sout[5] + sout[13] - sout[29] - sout[37] + sout[53] - sout[61]
   + xorout[20] - xorout[21] >= -4
 R965: - sout[5] - sout[13] + sout[29] - sout[37] + sout[53] - sout[61]
   + xorout[20] - xorout[21] >= -4
 R966: - sout[5] - sout[13] - sout[29] + sout[37] + sout[53] - sout[61]
   + xorout[20] - xorout[21] >= -4
 R967: sout[5] - sout[13] - sout[29] - sout[37] - sout[53] + sout[61]
   + xorout[20] - xorout[21] >= -4
 R968: - sout[5] + sout[13] - sout[29] - sout[37] - sout[53] + sout[61]
   + xorout[20] - xorout[21] >= -4
 R969: - sout[5] - sout[13] + sout[29] - sout[37] - sout[53] + sout[61]
   + xorout[20] - xorout[21] >= -4
 R970: - sout[5] - sout[13] - sout[29] + sout[37] - sout[53] + sout[61]
   + xorout[20] - xorout[21] >= -4
 R971: - sout[5] - sout[13] - sout[29] - sout[37] + sout[53] + sout[61]
   + xorout[20] - xorout[21] >= -4
 R972: sout[5] + sout[13] + sout[29] + sout[37] - sout[53] - sout[61]
   + xorout[20] - xorout[21] >= -2
 R973: sout[5] + sout[13] + sout[29] - sout[37] + sout[53] - sout[61]
   + xorout[20] - xorout[21] >= -2
 R974: sout[5] + sout[13] - sout[29] + sout[37] + sout[53] - sout[61]
   + xorout[20] - xorout[21] >= -2
 R975: sout[5] - sout[13] + sout[29] + sout[37] + sout[53] - sout[61]
   + xorout[20] - xorout[21] >= -2
 R976: - sout[5] + sout[13] + sout[29] + sout[37] + sout[53] - sout[61]
   + xorout[20] - xorout[21] >= -2
 R977: sout[5] + sout[13] + sout[29] - sout[37] - sout[53] + sout[61]
   + xorout[20] - xorout[21] >= -2
 R978: sout[5] + sout[13] - sout[29] + sout[37] - sout[53] + sout[61]
   + xorout[20] - xorout[21] >= -2
 R979: sout[5] - sout[13] + sout[29] + sout[37] - sout[53] + sout[61]
   + xorout[20] - xorout[21] >= -2
 R980: - sout[5] + sout[13] + sout[29] + sout[37] - sout[53] + sout[61]
   + xorout[20] - xorout[21] >= -2
 R981: sout[5] + sout[13] - sout[29] - sout[37] + sout[53] + sout[61]
   + xorout[20] - xorout[21] >= -2
 R982: sout[5] - sout[13] + sout[29] - sout[37] + sout[53] + sout[61]
   + xorout[20] - xorout[21] >= -2
 R983: - sout[5] + sout[13] + sout[29] - sout[37] + sout[53] + sout[61]
   + xorout[20] - xorout[21] >= -2
 R984: sout[5] - sout[13] - sout[29] + sout[37] + sout[53] + sout[61]
   + xorout[20] - xorout[21] >= -2
 R985: - sout[5] + sout[13] - sout[29] + sout[37] + sout[53] + sout[61]
   + xorout[20] - xorout[21] >= -2
 R986: - sout[5] - sout[13] + sout[29] + sout[37] + sout[53] + sout[61]
   + xorout[20] - xorout[21] >= -2
 R987: sout[7] + sout[15] + sout[31] + sout[39] + sout[55] - sout[63]
   + xorout[23] >= 0
 R988: sout[7] + sout[15] + sout[31] + sout[39] - sout[55] + sout[63]
   + xorout[23] >= 0
 R989: sout[7] + sout[15] + sout[31] - sout[39] + sout[55] + sout[63]
   + xorout[23] >= 0
 R990: sout[7] + sout[15] - sout[31] + sout[39] + sout[55] + sout[63]
   + xorout[23] >= 0
 R991: sout[7] - sout[15] + sout[31] + sout[39] + sout[55] + sout[63]
   + xorout[23] >= 0
 R992: - sout[7] + sout[15] + sout[31] + sout[39] + sout[55] + sout[63]
   + xorout[23] >= 0
 R993: sout[7] + sout[15] + sout[31] + sout[39] + sout[55] + sout[63]
   - xorout[23] >= 0
 R994: sout[7] + sout[15] + sout[31] - sout[39] - sout[55] - sout[63]
   + xorout[23] >= -2
 R995: sout[7] + sout[15] - sout[31] + sout[39] - sout[55] - sout[63]
   + xorout[23] >= -2
 R996: sout[7] - sout[15] + sout[31] + sout[39] - sout[55] - sout[63]
   + xorout[23] >= -2
 R997: - sout[7] + sout[15] + sout[31] + sout[39] - sout[55] - sout[63]
   + xorout[23] >= -2
 R998: sout[7] + sout[15] - sout[31] - sout[39] + sout[55] - sout[63]
   + xorout[23] >= -2
 R999: sout[7] - sout[15] + sout[31] - sout[39] + sout[55] - sout[63]
   + xorout[23] >= -2
 R1000: - sout[7] + sout[15] + sout[31] - sout[39] + sout[55] - sout[63]
   + xorout[23] >= -2
 R1001: sout[7] - sout[15] - sout[31] + sout[39] + sout[55] - sout[63]
   + xorout[23] >= -2
 R1002: - sout[7] + sout[15] - sout[31] + sout[39] + sout[55] - sout[63]
   + xorout[23] >= -2
 R1003: - sout[7] - sout[15] + sout[31] + sout[39] + sout[55] - sout[63]
   + xorout[23] >= -2
 R1004: sout[7] + sout[15] - sout[31] - sout[39] - sout[55] + sout[63]
   + xorout[23] >= -2
 R1005: sout[7] - sout[15] + sout[31] - sout[39] - sout[55] + sout[63]
   + xorout[23] >= -2
 R1006: - sout[7] + sout[15] + sout[31] - sout[39] - sout[55] + sout[63]
   + xorout[23] >= -2
 R1007: sout[7] - sout[15] - sout[31] + sout[39] - sout[55] + sout[63]
   + xorout[23] >= -2
 R1008: - sout[7] + sout[15] - sout[31] + sout[39] - sout[55] + sout[63]
   + xorout[23] >= -2
 R1009: - sout[7] - sout[15] + sout[31] + sout[39] - sout[55] + sout[63]
   + xorout[23] >= -2
 R1010: sout[7] - sout[15] - sout[31] - sout[39] + sout[55] + sout[63]
   + xorout[23] >= -2
 R1011: - sout[7] + sout[15] - sout[31] - sout[39] + sout[55] + sout[63]
   + xorout[23] >= -2
 R1012: - sout[7] - sout[15] + sout[31] - sout[39] + sout[55] + sout[63]
   + xorout[23] >= -2
 R1013: - sout[7] - sout[15] - sout[31] + sout[39] + sout[55] + sout[63]
   + xorout[23] >= -2
 R1014: sout[7] - sout[15] - sout[31] - sout[39] - sout[55] - sout[63]
   + xorout[23] >= -4
 R1015: - sout[7] + sout[15] - sout[31] - sout[39] - sout[55] - sout[63]
   + xorout[23] >= -4
 R1016: - sout[7] - sout[15] + sout[31] - sout[39] - sout[55] - sout[63]
   + xorout[23] >= -4
 R1017: - sout[7] - sout[15] - sout[31] + sout[39] - sout[55] - sout[63]
   + xorout[23] >= -4
 R1018: - sout[7] - sout[15] - sout[31] - sout[39] + sout[55] - sout[63]
   + xorout[23] >= -4
 R1019: - sout[7] - sout[15] - sout[31] - sout[39] - sout[55] + sout[63]
   + xorout[23] >= -4
 R1020: - sout[62] + xorout[22] >= 0
 R1021: - sout[54] + xorout[22] >= 0
 R1022: - sout[38] + xorout[22] >= 0
 R1023: - sout[30] + xorout[22] >= 0
 R1024: - sout[14] + xorout[22] >= 0
 R1025: - sout[6] + xorout[22] >= 0
 R1026: - sout[6] + sout[7] >= 0
 R1027: - sout[14] + sout[15] >= 0
 R1028: - sout[30] + sout[31] >= 0
 R1029: - sout[38] + sout[39] >= 0
 R1030: - sout[54] + sout[55] >= 0
 R1031: - sout[62] + sout[63] >= 0
 R1032: sout[6] + sout[14] + sout[30] + sout[38] + sout[54] + sout[62]
   - xorout[22] >= 0
 R1033: - xorout[22] + xorout[23] >= 0
 R1034: - sout[7] - sout[15] - sout[31] - sout[39] - sout[55] - sout[63]
   + xorout[22] - xorout[23] >= -6
 R1035: sout[7] + sout[15] - sout[31] - sout[39] - sout[55] - sout[63]
   + xorout[22] - xorout[23] >= -4
 R1036: sout[7] - sout[15] + sout[31] - sout[39] - sout[55] - sout[63]
   + xorout[22] - xorout[23] >= -4
 R1037: - sout[7] + sout[15] + sout[31] - sout[39] - sout[55] - sout[63]
   + xorout[22] - xorout[23] >= -4
 R1038: sout[7] - sout[15] - sout[31] + sout[39] - sout[55] - sout[63]
   + xorout[22] - xorout[23] >= -4
 R1039: - sout[7] + sout[15] - sout[31] + sout[39] - sout[55] - sout[63]
   + xorout[22] - xorout[23] >= -4
 R1040: - sout[7] - sout[15] + sout[31] + sout[39] - sout[55] - sout[63]
   + xorout[22] - xorout[23] >= -4
 R1041: sout[7] - sout[15] - sout[31] - sout[39] + sout[55] - sout[63]
   + xorout[22] - xorout[23] >= -4
 R1042: - sout[7] + sout[15] - sout[31] - sout[39] + sout[55] - sout[63]
   + xorout[22] - xorout[23] >= -4
 R1043: - sout[7] - sout[15] + sout[31] - sout[39] + sout[55] - sout[63]
   + xorout[22] - xorout[23] >= -4
 R1044: - sout[7] - sout[15] - sout[31] + sout[39] + sout[55] - sout[63]
   + xorout[22] - xorout[23] >= -4
 R1045: sout[7] - sout[15] - sout[31] - sout[39] - sout[55] + sout[63]
   + xorout[22] - xorout[23] >= -4
 R1046: - sout[7] + sout[15] - sout[31] - sout[39] - sout[55] + sout[63]
   + xorout[22] - xorout[23] >= -4
 R1047: - sout[7] - sout[15] + sout[31] - sout[39] - sout[55] + sout[63]
   + xorout[22] - xorout[23] >= -4
 R1048: - sout[7] - sout[15] - sout[31] + sout[39] - sout[55] + sout[63]
   + xorout[22] - xorout[23] >= -4
 R1049: - sout[7] - sout[15] - sout[31] - sout[39] + sout[55] + sout[63]
   + xorout[22] - xorout[23] >= -4
 R1050: sout[7] + sout[15] + sout[31] + sout[39] - sout[55] - sout[63]
   + xorout[22] - xorout[23] >= -2
 R1051: sout[7] + sout[15] + sout[31] - sout[39] + sout[55] - sout[63]
   + xorout[22] - xorout[23] >= -2
 R1052: sout[7] + sout[15] - sout[31] + sout[39] + sout[55] - sout[63]
   + xorout[22] - xorout[23] >= -2
 R1053: sout[7] - sout[15] + sout[31] + sout[39] + sout[55] - sout[63]
   + xorout[22] - xorout[23] >= -2
 R1054: - sout[7] + sout[15] + sout[31] + sout[39] + sout[55] - sout[63]
   + xorout[22] - xorout[23] >= -2
 R1055: sout[7] + sout[15] + sout[31] - sout[39] - sout[55] + sout[63]
   + xorout[22] - xorout[23] >= -2
 R1056: sout[7] + sout[15] - sout[31] + sout[39] - sout[55] + sout[63]
   + xorout[22] - xorout[23] >= -2
 R1057: sout[7] - sout[15] + sout[31] + sout[39] - sout[55] + sout[63]
   + xorout[22] - xorout[23] >= -2
 R1058: - sout[7] + sout[15] + sout[31] + sout[39] - sout[55] + sout[63]
   + xorout[22] - xorout[23] >= -2
 R1059: sout[7] + sout[15] - sout[31] - sout[39] + sout[55] + sout[63]
   + xorout[22] - xorout[23] >= -2
 R1060: sout[7] - sout[15] + sout[31] - sout[39] + sout[55] + sout[63]
   + xorout[22] - xorout[23] >= -2
 R1061: - sout[7] + sout[15] + sout[31] - sout[39] + sout[55] + sout[63]
   + xorout[22] - xorout[23] >= -2
 R1062: sout[7] - sout[15] - sout[31] + sout[39] + sout[55] + sout[63]
   + xorout[22] - xorout[23] >= -2
 R1063: - sout[7] + sout[15] - sout[31] + sout[39] + sout[55] + sout[63]
   + xorout[22] - xorout[23] >= -2
 R1064: - sout[7] - sout[15] + sout[31] + sout[39] + sout[55] + sout[63]
   + xorout[22] - xorout[23] >= -2
 R1065: sout[1] + sout[9] + sout[17] + sout[33] + sout[41] - sout[57]
   + xorout[25] >= 0
 R1066: sout[1] + sout[9] + sout[17] + sout[33] - sout[41] + sout[57]
   + xorout[25] >= 0
 R1067: sout[1] + sout[9] + sout[17] - sout[33] + sout[41] + sout[57]
   + xorout[25] >= 0
 R1068: sout[1] + sout[9] - sout[17] + sout[33] + sout[41] + sout[57]
   + xorout[25] >= 0
 R1069: sout[1] - sout[9] + sout[17] + sout[33] + sout[41] + sout[57]
   + xorout[25] >= 0
 R1070: - sout[1] + sout[9] + sout[17] + sout[33] + sout[41] + sout[57]
   + xorout[25] >= 0
 R1071: sout[1] + sout[9] + sout[17] + sout[33] + sout[41] + sout[57]
   - xorout[25] >= 0
 R1072: sout[1] + sout[9] + sout[17] - sout[33] - sout[41] - sout[57]
   + xorout[25] >= -2
 R1073: sout[1] + sout[9] - sout[17] + sout[33] - sout[41] - sout[57]
   + xorout[25] >= -2
 R1074: sout[1] - sout[9] + sout[17] + sout[33] - sout[41] - sout[57]
   + xorout[25] >= -2
 R1075: - sout[1] + sout[9] + sout[17] + sout[33] - sout[41] - sout[57]
   + xorout[25] >= -2
 R1076: sout[1] + sout[9] - sout[17] - sout[33] + sout[41] - sout[57]
   + xorout[25] >= -2
 R1077: sout[1] - sout[9] + sout[17] - sout[33] + sout[41] - sout[57]
   + xorout[25] >= -2
 R1078: - sout[1] + sout[9] + sout[17] - sout[33] + sout[41] - sout[57]
   + xorout[25] >= -2
 R1079: sout[1] - sout[9] - sout[17] + sout[33] + sout[41] - sout[57]
   + xorout[25] >= -2
 R1080: - sout[1] + sout[9] - sout[17] + sout[33] + sout[41] - sout[57]
   + xorout[25] >= -2
 R1081: - sout[1] - sout[9] + sout[17] + sout[33] + sout[41] - sout[57]
   + xorout[25] >= -2
 R1082: sout[1] + sout[9] - sout[17] - sout[33] - sout[41] + sout[57]
   + xorout[25] >= -2
 R1083: sout[1] - sout[9] + sout[17] - sout[33] - sout[41] + sout[57]
   + xorout[25] >= -2
 R1084: - sout[1] + sout[9] + sout[17] - sout[33] - sout[41] + sout[57]
   + xorout[25] >= -2
 R1085: sout[1] - sout[9] - sout[17] + sout[33] - sout[41] + sout[57]
   + xorout[25] >= -2
 R1086: - sout[1] + sout[9] - sout[17] + sout[33] - sout[41] + sout[57]
   + xorout[25] >= -2
 R1087: - sout[1] - sout[9] + sout[17] + sout[33] - sout[41] + sout[57]
   + xorout[25] >= -2
 R1088: sout[1] - sout[9] - sout[17] - sout[33] + sout[41] + sout[57]
   + xorout[25] >= -2
 R1089: - sout[1] + sout[9] - sout[17] - sout[33] + sout[41] + sout[57]
   + xorout[25] >= -2
 R1090: - sout[1] - sout[9] + sout[17] - sout[33] + sout[41] + sout[57]
   + xorout[25] >= -2
 R1091: - sout[1] - sout[9] - sout[17] + sout[33] + sout[41] + sout[57]
   + xorout[25] >= -2
 R1092: sout[1] - sout[9] - sout[17] - sout[33] - sout[41] - sout[57]
   + xorout[25] >= -4
 R1093: - sout[1] + sout[9] - sout[17] - sout[33] - sout[41] - sout[57]
   + xorout[25] >= -4
 R1094: - sout[1] - sout[9] + sout[17] - sout[33] - sout[41] - sout[57]
   + xorout[25] >= -4
 R1095: - sout[1] - sout[9] - sout[17] + sout[33] - sout[41] - sout[57]
   + xorout[25] >= -4
 R1096: - sout[1] - sout[9] - sout[17] - sout[33] + sout[41] - sout[57]
   + xorout[25] >= -4
 R1097: - sout[1] - sout[9] - sout[17] - sout[33] - sout[41] + sout[57]
   + xorout[25] >= -4
 R1098: - sout[56] + xorout[24] >= 0
 R1099: - sout[40] + xorout[24] >= 0
 R1100: - sout[32] + xorout[24] >= 0
 R1101: - sout[16] + xorout[24] >= 0
 R1102: - sout[8] + xorout[24] >= 0
 R1103: - sout[0] + xorout[24] >= 0
 R1104: - sout[0] + sout[1] >= 0
 R1105: - sout[8] + sout[9] >= 0
 R1106: - sout[16] + sout[17] >= 0
 R1107: - sout[32] + sout[33] >= 0
 R1108: - sout[40] + sout[41] >= 0
 R1109: - sout[56] + sout[57] >= 0
 R1110: sout[0] + sout[8] + sout[16] + sout[32] + sout[40] + sout[56]
   - xorout[24] >= 0
 R1111: - xorout[24] + xorout[25] >= 0
 R1112: - sout[1] - sout[9] - sout[17] - sout[33] - sout[41] - sout[57]
   + xorout[24] - xorout[25] >= -6
 R1113: sout[1] + sout[9] - sout[17] - sout[33] - sout[41] - sout[57]
   + xorout[24] - xorout[25] >= -4
 R1114: sout[1] - sout[9] + sout[17] - sout[33] - sout[41] - sout[57]
   + xorout[24] - xorout[25] >= -4
 R1115: - sout[1] + sout[9] + sout[17] - sout[33] - sout[41] - sout[57]
   + xorout[24] - xorout[25] >= -4
 R1116: sout[1] - sout[9] - sout[17] + sout[33] - sout[41] - sout[57]
   + xorout[24] - xorout[25] >= -4
 R1117: - sout[1] + sout[9] - sout[17] + sout[33] - sout[41] - sout[57]
   + xorout[24] - xorout[25] >= -4
 R1118: - sout[1] - sout[9] + sout[17] + sout[33] - sout[41] - sout[57]
   + xorout[24] - xorout[25] >= -4
 R1119: sout[1] - sout[9] - sout[17] - sout[33] + sout[41] - sout[57]
   + xorout[24] - xorout[25] >= -4
 R1120: - sout[1] + sout[9] - sout[17] - sout[33] + sout[41] - sout[57]
   + xorout[24] - xorout[25] >= -4
 R1121: - sout[1] - sout[9] + sout[17] - sout[33] + sout[41] - sout[57]
   + xorout[24] - xorout[25] >= -4
 R1122: - sout[1] - sout[9] - sout[17] + sout[33] + sout[41] - sout[57]
   + xorout[24] - xorout[25] >= -4
 R1123: sout[1] - sout[9] - sout[17] - sout[33] - sout[41] + sout[57]
   + xorout[24] - xorout[25] >= -4
 R1124: - sout[1] + sout[9] - sout[17] - sout[33] - sout[41] + sout[57]
   + xorout[24] - xorout[25] >= -4
 R1125: - sout[1] - sout[9] + sout[17] - sout[33] - sout[41] + sout[57]
   + xorout[24] - xorout[25] >= -4
 R1126: - sout[1] - sout[9] - sout[17] + sout[33] - sout[41] + sout[57]
   + xorout[24] - xorout[25] >= -4
 R1127: - sout[1] - sout[9] - sout[17] - sout[33] + sout[41] + sout[57]
   + xorout[24] - xorout[25] >= -4
 R1128: sout[1] + sout[9] + sout[17] + sout[33] - sout[41] - sout[57]
   + xorout[24] - xorout[25] >= -2
 R1129: sout[1] + sout[9] + sout[17] - sout[33] + sout[41] - sout[57]
   + xorout[24] - xorout[25] >= -2
 R1130: sout[1] + sout[9] - sout[17] + sout[33] + sout[41] - sout[57]
   + xorout[24] - xorout[25] >= -2
 R1131: sout[1] - sout[9] + sout[17] + sout[33] + sout[41] - sout[57]
   + xorout[24] - xorout[25] >= -2
 R1132: - sout[1] + sout[9] + sout[17] + sout[33] + sout[41] - sout[57]
   + xorout[24] - xorout[25] >= -2
 R1133: sout[1] + sout[9] + sout[17] - sout[33] - sout[41] + sout[57]
   + xorout[24] - xorout[25] >= -2
 R1134: sout[1] + sout[9] - sout[17] + sout[33] - sout[41] + sout[57]
   + xorout[24] - xorout[25] >= -2
 R1135: sout[1] - sout[9] + sout[17] + sout[33] - sout[41] + sout[57]
   + xorout[24] - xorout[25] >= -2
 R1136: - sout[1] + sout[9] + sout[17] + sout[33] - sout[41] + sout[57]
   + xorout[24] - xorout[25] >= -2
 R1137: sout[1] + sout[9] - sout[17] - sout[33] + sout[41] + sout[57]
   + xorout[24] - xorout[25] >= -2
 R1138: sout[1] - sout[9] + sout[17] - sout[33] + sout[41] + sout[57]
   + xorout[24] - xorout[25] >= -2
 R1139: - sout[1] + sout[9] + sout[17] - sout[33] + sout[41] + sout[57]
   + xorout[24] - xorout[25] >= -2
 R1140: sout[1] - sout[9] - sout[17] + sout[33] + sout[41] + sout[57]
   + xorout[24] - xorout[25] >= -2
 R1141: - sout[1] + sout[9] - sout[17] + sout[33] + sout[41] + sout[57]
   + xorout[24] - xorout[25] >= -2
 R1142: - sout[1] - sout[9] + sout[17] + sout[33] + sout[41] + sout[57]
   + xorout[24] - xorout[25] >= -2
 R1143: sout[3] + sout[11] + sout[19] + sout[35] + sout[43] - sout[59]
   + xorout[27] >= 0
 R1144: sout[3] + sout[11] + sout[19] + sout[35] - sout[43] + sout[59]
   + xorout[27] >= 0
 R1145: sout[3] + sout[11] + sout[19] - sout[35] + sout[43] + sout[59]
   + xorout[27] >= 0
 R1146: sout[3] + sout[11] - sout[19] + sout[35] + sout[43] + sout[59]
   + xorout[27] >= 0
 R1147: sout[3] - sout[11] + sout[19] + sout[35] + sout[43] + sout[59]
   + xorout[27] >= 0
 R1148: - sout[3] + sout[11] + sout[19] + sout[35] + sout[43] + sout[59]
   + xorout[27] >= 0
 R1149: sout[3] + sout[11] + sout[19] + sout[35] + sout[43] + sout[59]
   - xorout[27] >= 0
 R1150: sout[3] + sout[11] + sout[19] - sout[35] - sout[43] - sout[59]
   + xorout[27] >= -2
 R1151: sout[3] + sout[11] - sout[19] + sout[35] - sout[43] - sout[59]
   + xorout[27] >= -2
 R1152: sout[3] - sout[11] + sout[19] + sout[35] - sout[43] - sout[59]
   + xorout[27] >= -2
 R1153: - sout[3] + sout[11] + sout[19] + sout[35] - sout[43] - sout[59]
   + xorout[27] >= -2
 R1154: sout[3] + sout[11] - sout[19] - sout[35] + sout[43] - sout[59]
   + xorout[27] >= -2
 R1155: sout[3] - sout[11] + sout[19] - sout[35] + sout[43] - sout[59]
   + xorout[27] >= -2
 R1156: - sout[3] + sout[11] + sout[19] - sout[35] + sout[43] - sout[59]
   + xorout[27] >= -2
 R1157: sout[3] - sout[11] - sout[19] + sout[35] + sout[43] - sout[59]
   + xorout[27] >= -2
 R1158: - sout[3] + sout[11] - sout[19] + sout[35] + sout[43] - sout[59]
   + xorout[27] >= -2
 R1159: - sout[3] - sout[11] + sout[19] + sout[35] + sout[43] - sout[59]
   + xorout[27] >= -2
 R1160: sout[3] + sout[11] - sout[19] - sout[35] - sout[43] + sout[59]
   + xorout[27] >= -2
 R1161: sout[3] - sout[11] + sout[19] - sout[35] - sout[43] + sout[59]
   + xorout[27] >= -2
 R1162: - sout[3] + sout[11] + sout[19] - sout[35] - sout[43] + sout[59]
   + xorout[27] >= -2
 R1163: sout[3] - sout[11] - sout[19] + sout[35] - sout[43] + sout[59]
   + xorout[27] >= -2
 R1164: - sout[3] + sout[11] - sout[19] + sout[35] - sout[43] + sout[59]
   + xorout[27] >= -2
 R1165: - sout[3] - sout[11] + sout[19] + sout[35] - sout[43] + sout[59]
   + xorout[27] >= -2
 R1166: sout[3] - sout[11] - sout[19] - sout[35] + sout[43] + sout[59]
   + xorout[27] >= -2
 R1167: - sout[3] + sout[11] - sout[19] - sout[35] + sout[43] + sout[59]
   + xorout[27] >= -2
 R1168: - sout[3] - sout[11] + sout[19] - sout[35] + sout[43] + sout[59]
   + xorout[27] >= -2
 R1169: - sout[3] - sout[11] - sout[19] + sout[35] + sout[43] + sout[59]
   + xorout[27] >= -2
 R1170: sout[3] - sout[11] - sout[19] - sout[35] - sout[43] - sout[59]
   + xorout[27] >= -4
 R1171: - sout[3] + sout[11] - sout[19] - sout[35] - sout[43] - sout[59]
   + xorout[27] >= -4
 R1172: - sout[3] - sout[11] + sout[19] - sout[35] - sout[43] - sout[59]
   + xorout[27] >= -4
 R1173: - sout[3] - sout[11] - sout[19] + sout[35] - sout[43] - sout[59]
   + xorout[27] >= -4
 R1174: - sout[3] - sout[11] - sout[19] - sout[35] + sout[43] - sout[59]
   + xorout[27] >= -4
 R1175: - sout[3] - sout[11] - sout[19] - sout[35] - sout[43] + sout[59]
   + xorout[27] >= -4
 R1176: - sout[58] + xorout[26] >= 0
 R1177: - sout[42] + xorout[26] >= 0
 R1178: - sout[34] + xorout[26] >= 0
 R1179: - sout[18] + xorout[26] >= 0
 R1180: - sout[10] + xorout[26] >= 0
 R1181: - sout[2] + xorout[26] >= 0
 R1182: - sout[2] + sout[3] >= 0
 R1183: - sout[10] + sout[11] >= 0
 R1184: - sout[18] + sout[19] >= 0
 R1185: - sout[34] + sout[35] >= 0
 R1186: - sout[42] + sout[43] >= 0
 R1187: - sout[58] + sout[59] >= 0
 R1188: sout[2] + sout[10] + sout[18] + sout[34] + sout[42] + sout[58]
   - xorout[26] >= 0
 R1189: - xorout[26] + xorout[27] >= 0
 R1190: - sout[3] - sout[11] - sout[19] - sout[35] - sout[43] - sout[59]
   + xorout[26] - xorout[27] >= -6
 R1191: sout[3] + sout[11] - sout[19] - sout[35] - sout[43] - sout[59]
   + xorout[26] - xorout[27] >= -4
 R1192: sout[3] - sout[11] + sout[19] - sout[35] - sout[43] - sout[59]
   + xorout[26] - xorout[27] >= -4
 R1193: - sout[3] + sout[11] + sout[19] - sout[35] - sout[43] - sout[59]
   + xorout[26] - xorout[27] >= -4
 R1194: sout[3] - sout[11] - sout[19] + sout[35] - sout[43] - sout[59]
   + xorout[26] - xorout[27] >= -4
 R1195: - sout[3] + sout[11] - sout[19] + sout[35] - sout[43] - sout[59]
   + xorout[26] - xorout[27] >= -4
 R1196: - sout[3] - sout[11] + sout[19] + sout[35] - sout[43] - sout[59]
   + xorout[26] - xorout[27] >= -4
 R1197: sout[3] - sout[11] - sout[19] - sout[35] + sout[43] - sout[59]
   + xorout[26] - xorout[27] >= -4
 R1198: - sout[3] + sout[11] - sout[19] - sout[35] + sout[43] - sout[59]
   + xorout[26] - xorout[27] >= -4
 R1199: - sout[3] - sout[11] + sout[19] - sout[35] + sout[43] - sout[59]
   + xorout[26] - xorout[27] >= -4
 R1200: - sout[3] - sout[11] - sout[19] + sout[35] + sout[43] - sout[59]
   + xorout[26] - xorout[27] >= -4
 R1201: sout[3] - sout[11] - sout[19] - sout[35] - sout[43] + sout[59]
   + xorout[26] - xorout[27] >= -4
 R1202: - sout[3] + sout[11] - sout[19] - sout[35] - sout[43] + sout[59]
   + xorout[26] - xorout[27] >= -4
 R1203: - sout[3] - sout[11] + sout[19] - sout[35] - sout[43] + sout[59]
   + xorout[26] - xorout[27] >= -4
 R1204: - sout[3] - sout[11] - sout[19] + sout[35] - sout[43] + sout[59]
   + xorout[26] - xorout[27] >= -4
 R1205: - sout[3] - sout[11] - sout[19] - sout[35] + sout[43] + sout[59]
   + xorout[26] - xorout[27] >= -4
 R1206: sout[3] + sout[11] + sout[19] + sout[35] - sout[43] - sout[59]
   + xorout[26] - xorout[27] >= -2
 R1207: sout[3] + sout[11] + sout[19] - sout[35] + sout[43] - sout[59]
   + xorout[26] - xorout[27] >= -2
 R1208: sout[3] + sout[11] - sout[19] + sout[35] + sout[43] - sout[59]
   + xorout[26] - xorout[27] >= -2
 R1209: sout[3] - sout[11] + sout[19] + sout[35] + sout[43] - sout[59]
   + xorout[26] - xorout[27] >= -2
 R1210: - sout[3] + sout[11] + sout[19] + sout[35] + sout[43] - sout[59]
   + xorout[26] - xorout[27] >= -2
 R1211: sout[3] + sout[11] + sout[19] - sout[35] - sout[43] + sout[59]
   + xorout[26] - xorout[27] >= -2
 R1212: sout[3] + sout[11] - sout[19] + sout[35] - sout[43] + sout[59]
   + xorout[26] - xorout[27] >= -2
 R1213: sout[3] - sout[11] + sout[19] + sout[35] - sout[43] + sout[59]
   + xorout[26] - xorout[27] >= -2
 R1214: - sout[3] + sout[11] + sout[19] + sout[35] - sout[43] + sout[59]
   + xorout[26] - xorout[27] >= -2
 R1215: sout[3] + sout[11] - sout[19] - sout[35] + sout[43] + sout[59]
   + xorout[26] - xorout[27] >= -2
 R1216: sout[3] - sout[11] + sout[19] - sout[35] + sout[43] + sout[59]
   + xorout[26] - xorout[27] >= -2
 R1217: - sout[3] + sout[11] + sout[19] - sout[35] + sout[43] + sout[59]
   + xorout[26] - xorout[27] >= -2
 R1218: sout[3] - sout[11] - sout[19] + sout[35] + sout[43] + sout[59]
   + xorout[26] - xorout[27] >= -2
 R1219: - sout[3] + sout[11] - sout[19] + sout[35] + sout[43] + sout[59]
   + xorout[26] - xorout[27] >= -2
 R1220: - sout[3] - sout[11] + sout[19] + sout[35] + sout[43] + sout[59]
   + xorout[26] - xorout[27] >= -2
 R1221: sout[5] + sout[13] + sout[21] + sout[37] + sout[45] - sout[61]
   + xorout[29] >= 0
 R1222: sout[5] + sout[13] + sout[21] + sout[37] - sout[45] + sout[61]
   + xorout[29] >= 0
 R1223: sout[5] + sout[13] + sout[21] - sout[37] + sout[45] + sout[61]
   + xorout[29] >= 0
 R1224: sout[5] + sout[13] - sout[21] + sout[37] + sout[45] + sout[61]
   + xorout[29] >= 0
 R1225: sout[5] - sout[13] + sout[21] + sout[37] + sout[45] + sout[61]
   + xorout[29] >= 0
 R1226: - sout[5] + sout[13] + sout[21] + sout[37] + sout[45] + sout[61]
   + xorout[29] >= 0
 R1227: sout[5] + sout[13] + sout[21] + sout[37] + sout[45] + sout[61]
   - xorout[29] >= 0
 R1228: sout[5] + sout[13] + sout[21] - sout[37] - sout[45] - sout[61]
   + xorout[29] >= -2
 R1229: sout[5] + sout[13] - sout[21] + sout[37] - sout[45] - sout[61]
   + xorout[29] >= -2
 R1230: sout[5] - sout[13] + sout[21] + sout[37] - sout[45] - sout[61]
   + xorout[29] >= -2
 R1231: - sout[5] + sout[13] + sout[21] + sout[37] - sout[45] - sout[61]
   + xorout[29] >= -2
 R1232: sout[5] + sout[13] - sout[21] - sout[37] + sout[45] - sout[61]
   + xorout[29] >= -2
 R1233: sout[5] - sout[13] + sout[21] - sout[37] + sout[45] - sout[61]
   + xorout[29] >= -2
 R1234: - sout[5] + sout[13] + sout[21] - sout[37] + sout[45] - sout[61]
   + xorout[29] >= -2
 R1235: sout[5] - sout[13] - sout[21] + sout[37] + sout[45] - sout[61]
   + xorout[29] >= -2
 R1236: - sout[5] + sout[13] - sout[21] + sout[37] + sout[45] - sout[61]
   + xorout[29] >= -2
 R1237: - sout[5] - sout[13] + sout[21] + sout[37] + sout[45] - sout[61]
   + xorout[29] >= -2
 R1238: sout[5] + sout[13] - sout[21] - sout[37] - sout[45] + sout[61]
   + xorout[29] >= -2
 R1239: sout[5] - sout[13] + sout[21] - sout[37] - sout[45] + sout[61]
   + xorout[29] >= -2
 R1240: - sout[5] + sout[13] + sout[21] - sout[37] - sout[45] + sout[61]
   + xorout[29] >= -2
 R1241: sout[5] - sout[13] - sout[21] + sout[37] - sout[45] + sout[61]
   + xorout[29] >= -2
 R1242: - sout[5] + sout[13] - sout[21] + sout[37] - sout[45] + sout[61]
   + xorout[29] >= -2
 R1243: - sout[5] - sout[13] + sout[21] + sout[37] - sout[45] + sout[61]
   + xorout[29] >= -2
 R1244: sout[5] - sout[13] - sout[21] - sout[37] + sout[45] + sout[61]
   + xorout[29] >= -2
 R1245: - sout[5] + sout[13] - sout[21] - sout[37] + sout[45] + sout[61]
   + xorout[29] >= -2
 R1246: - sout[5] - sout[13] + sout[21] - sout[37] + sout[45] + sout[61]
   + xorout[29] >= -2
 R1247: - sout[5] - sout[13] - sout[21] + sout[37] + sout[45] + sout[61]
   + xorout[29] >= -2
 R1248: sout[5] - sout[13] - sout[21] - sout[37] - sout[45] - sout[61]
   + xorout[29] >= -4
 R1249: - sout[5] + sout[13] - sout[21] - sout[37] - sout[45] - sout[61]
   + xorout[29] >= -4
 R1250: - sout[5] - sout[13] + sout[21] - sout[37] - sout[45] - sout[61]
   + xorout[29] >= -4
 R1251: - sout[5] - sout[13] - sout[21] + sout[37] - sout[45] - sout[61]
   + xorout[29] >= -4
 R1252: - sout[5] - sout[13] - sout[21] - sout[37] + sout[45] - sout[61]
   + xorout[29] >= -4
 R1253: - sout[5] - sout[13] - sout[21] - sout[37] - sout[45] + sout[61]
   + xorout[29] >= -4
 R1254: - sout[60] + xorout[28] >= 0
 R1255: - sout[44] + xorout[28] >= 0
 R1256: - sout[36] + xorout[28] >= 0
 R1257: - sout[20] + xorout[28] >= 0
 R1258: - sout[12] + xorout[28] >= 0
 R1259: - sout[4] + xorout[28] >= 0
 R1260: - sout[4] + sout[5] >= 0
 R1261: - sout[12] + sout[13] >= 0
 R1262: - sout[20] + sout[21] >= 0
 R1263: - sout[36] + sout[37] >= 0
 R1264: - sout[44] + sout[45] >= 0
 R1265: - sout[60] + sout[61] >= 0
 R1266: sout[4] + sout[12] + sout[20] + sout[36] + sout[44] + sout[60]
   - xorout[28] >= 0
 R1267: - xorout[28] + xorout[29] >= 0
 R1268: - sout[5] - sout[13] - sout[21] - sout[37] - sout[45] - sout[61]
   + xorout[28] - xorout[29] >= -6
 R1269: sout[5] + sout[13] - sout[21] - sout[37] - sout[45] - sout[61]
   + xorout[28] - xorout[29] >= -4
 R1270: sout[5] - sout[13] + sout[21] - sout[37] - sout[45] - sout[61]
   + xorout[28] - xorout[29] >= -4
 R1271: - sout[5] + sout[13] + sout[21] - sout[37] - sout[45] - sout[61]
   + xorout[28] - xorout[29] >= -4
 R1272: sout[5] - sout[13] - sout[21] + sout[37] - sout[45] - sout[61]
   + xorout[28] - xorout[29] >= -4
 R1273: - sout[5] + sout[13] - sout[21] + sout[37] - sout[45] - sout[61]
   + xorout[28] - xorout[29] >= -4
 R1274: - sout[5] - sout[13] + sout[21] + sout[37] - sout[45] - sout[61]
   + xorout[28] - xorout[29] >= -4
 R1275: sout[5] - sout[13] - sout[21] - sout[37] + sout[45] - sout[61]
   + xorout[28] - xorout[29] >= -4
 R1276: - sout[5] + sout[13] - sout[21] - sout[37] + sout[45] - sout[61]
   + xorout[28] - xorout[29] >= -4
 R1277: - sout[5] - sout[13] + sout[21] - sout[37] + sout[45] - sout[61]
   + xorout[28] - xorout[29] >= -4
 R1278: - sout[5] - sout[13] - sout[21] + sout[37] + sout[45] - sout[61]
   + xorout[28] - xorout[29] >= -4
 R1279: sout[5] - sout[13] - sout[21] - sout[37] - sout[45] + sout[61]
   + xorout[28] - xorout[29] >= -4
 R1280: - sout[5] + sout[13] - sout[21] - sout[37] - sout[45] + sout[61]
   + xorout[28] - xorout[29] >= -4
 R1281: - sout[5] - sout[13] + sout[21] - sout[37] - sout[45] + sout[61]
   + xorout[28] - xorout[29] >= -4
 R1282: - sout[5] - sout[13] - sout[21] + sout[37] - sout[45] + sout[61]
   + xorout[28] - xorout[29] >= -4
 R1283: - sout[5] - sout[13] - sout[21] - sout[37] + sout[45] + sout[61]
   + xorout[28] - xorout[29] >= -4
 R1284: sout[5] + sout[13] + sout[21] + sout[37] - sout[45] - sout[61]
   + xorout[28] - xorout[29] >= -2
 R1285: sout[5] + sout[13] + sout[21] - sout[37] + sout[45] - sout[61]
   + xorout[28] - xorout[29] >= -2
 R1286: sout[5] + sout[13] - sout[21] + sout[37] + sout[45] - sout[61]
   + xorout[28] - xorout[29] >= -2
 R1287: sout[5] - sout[13] + sout[21] + sout[37] + sout[45] - sout[61]
   + xorout[28] - xorout[29] >= -2
 R1288: - sout[5] + sout[13] + sout[21] + sout[37] + sout[45] - sout[61]
   + xorout[28] - xorout[29] >= -2
 R1289: sout[5] + sout[13] + sout[21] - sout[37] - sout[45] + sout[61]
   + xorout[28] - xorout[29] >= -2
 R1290: sout[5] + sout[13] - sout[21] + sout[37] - sout[45] + sout[61]
   + xorout[28] - xorout[29] >= -2
 R1291: sout[5] - sout[13] + sout[21] + sout[37] - sout[45] + sout[61]
   + xorout[28] - xorout[29] >= -2
 R1292: - sout[5] + sout[13] + sout[21] + sout[37] - sout[45] + sout[61]
   + xorout[28] - xorout[29] >= -2
 R1293: sout[5] + sout[13] - sout[21] - sout[37] + sout[45] + sout[61]
   + xorout[28] - xorout[29] >= -2
 R1294: sout[5] - sout[13] + sout[21] - sout[37] + sout[45] + sout[61]
   + xorout[28] - xorout[29] >= -2
 R1295: - sout[5] + sout[13] + sout[21] - sout[37] + sout[45] + sout[61]
   + xorout[28] - xorout[29] >= -2
 R1296: sout[5] - sout[13] - sout[21] + sout[37] + sout[45] + sout[61]
   + xorout[28] - xorout[29] >= -2
 R1297: - sout[5] + sout[13] - sout[21] + sout[37] + sout[45] + sout[61]
   + xorout[28] - xorout[29] >= -2
 R1298: - sout[5] - sout[13] + sout[21] + sout[37] + sout[45] + sout[61]
   + xorout[28] - xorout[29] >= -2
 R1299: sout[7] + sout[15] + sout[23] + sout[39] + sout[47] - sout[63]
   + xorout[31] >= 0
 R1300: sout[7] + sout[15] + sout[23] + sout[39] - sout[47] + sout[63]
   + xorout[31] >= 0
 R1301: sout[7] + sout[15] + sout[23] - sout[39] + sout[47] + sout[63]
   + xorout[31] >= 0
 R1302: sout[7] + sout[15] - sout[23] + sout[39] + sout[47] + sout[63]
   + xorout[31] >= 0
 R1303: sout[7] - sout[15] + sout[23] + sout[39] + sout[47] + sout[63]
   + xorout[31] >= 0
 R1304: - sout[7] + sout[15] + sout[23] + sout[39] + sout[47] + sout[63]
   + xorout[31] >= 0
 R1305: sout[7] + sout[15] + sout[23] + sout[39] + sout[47] + sout[63]
   - xorout[31] >= 0
 R1306: sout[7] + sout[15] + sout[23] - sout[39] - sout[47] - sout[63]
   + xorout[31] >= -2
 R1307: sout[7] + sout[15] - sout[23] + sout[39] - sout[47] - sout[63]
   + xorout[31] >= -2
 R1308: sout[7] - sout[15] + sout[23] + sout[39] - sout[47] - sout[63]
   + xorout[31] >= -2
 R1309: - sout[7] + sout[15] + sout[23] + sout[39] - sout[47] - sout[63]
   + xorout[31] >= -2
 R1310: sout[7] + sout[15] - sout[23] - sout[39] + sout[47] - sout[63]
   + xorout[31] >= -2
 R1311: sout[7] - sout[15] + sout[23] - sout[39] + sout[47] - sout[63]
   + xorout[31] >= -2
 R1312: - sout[7] + sout[15] + sout[23] - sout[39] + sout[47] - sout[63]
   + xorout[31] >= -2
 R1313: sout[7] - sout[15] - sout[23] + sout[39] + sout[47] - sout[63]
   + xorout[31] >= -2
 R1314: - sout[7] + sout[15] - sout[23] + sout[39] + sout[47] - sout[63]
   + xorout[31] >= -2
 R1315: - sout[7] - sout[15] + sout[23] + sout[39] + sout[47] - sout[63]
   + xorout[31] >= -2
 R1316: sout[7] + sout[15] - sout[23] - sout[39] - sout[47] + sout[63]
   + xorout[31] >= -2
 R1317: sout[7] - sout[15] + sout[23] - sout[39] - sout[47] + sout[63]
   + xorout[31] >= -2
 R1318: - sout[7] + sout[15] + sout[23] - sout[39] - sout[47] + sout[63]
   + xorout[31] >= -2
 R1319: sout[7] - sout[15] - sout[23] + sout[39] - sout[47] + sout[63]
   + xorout[31] >= -2
 R1320: - sout[7] + sout[15] - sout[23] + sout[39] - sout[47] + sout[63]
   + xorout[31] >= -2
 R1321: - sout[7] - sout[15] + sout[23] + sout[39] - sout[47] + sout[63]
   + xorout[31] >= -2
 R1322: sout[7] - sout[15] - sout[23] - sout[39] + sout[47] + sout[63]
   + xorout[31] >= -2
 R1323: - sout[7] + sout[15] - sout[23] - sout[39] + sout[47] + sout[63]
   + xorout[31] >= -2
 R1324: - sout[7] - sout[15] + sout[23] - sout[39] + sout[47] + sout[63]
   + xorout[31] >= -2
 R1325: - sout[7] - sout[15] - sout[23] + sout[39] + sout[47] + sout[63]
   + xorout[31] >= -2
 R1326: sout[7] - sout[15] - sout[23] - sout[39] - sout[47] - sout[63]
   + xorout[31] >= -4
 R1327: - sout[7] + sout[15] - sout[23] - sout[39] - sout[47] - sout[63]
   + xorout[31] >= -4
 R1328: - sout[7] - sout[15] + sout[23] - sout[39] - sout[47] - sout[63]
   + xorout[31] >= -4
 R1329: - sout[7] - sout[15] - sout[23] + sout[39] - sout[47] - sout[63]
   + xorout[31] >= -4
 R1330: - sout[7] - sout[15] - sout[23] - sout[39] + sout[47] - sout[63]
   + xorout[31] >= -4
 R1331: - sout[7] - sout[15] - sout[23] - sout[39] - sout[47] + sout[63]
   + xorout[31] >= -4
 R1332: - sout[62] + xorout[30] >= 0
 R1333: - sout[46] + xorout[30] >= 0
 R1334: - sout[38] + xorout[30] >= 0
 R1335: - sout[22] + xorout[30] >= 0
 R1336: - sout[14] + xorout[30] >= 0
 R1337: - sout[6] + xorout[30] >= 0
 R1338: - sout[6] + sout[7] >= 0
 R1339: - sout[14] + sout[15] >= 0
 R1340: - sout[22] + sout[23] >= 0
 R1341: - sout[38] + sout[39] >= 0
 R1342: - sout[46] + sout[47] >= 0
 R1343: - sout[62] + sout[63] >= 0
 R1344: sout[6] + sout[14] + sout[22] + sout[38] + sout[46] + sout[62]
   - xorout[30] >= 0
 R1345: - xorout[30] + xorout[31] >= 0
 R1346: - sout[7] - sout[15] - sout[23] - sout[39] - sout[47] - sout[63]
   + xorout[30] - xorout[31] >= -6
 R1347: sout[7] + sout[15] - sout[23] - sout[39] - sout[47] - sout[63]
   + xorout[30] - xorout[31] >= -4
 R1348: sout[7] - sout[15] + sout[23] - sout[39] - sout[47] - sout[63]
   + xorout[30] - xorout[31] >= -4
 R1349: - sout[7] + sout[15] + sout[23] - sout[39] - sout[47] - sout[63]
   + xorout[30] - xorout[31] >= -4
 R1350: sout[7] - sout[15] - sout[23] + sout[39] - sout[47] - sout[63]
   + xorout[30] - xorout[31] >= -4
 R1351: - sout[7] + sout[15] - sout[23] + sout[39] - sout[47] - sout[63]
   + xorout[30] - xorout[31] >= -4
 R1352: - sout[7] - sout[15] + sout[23] + sout[39] - sout[47] - sout[63]
   + xorout[30] - xorout[31] >= -4
 R1353: sout[7] - sout[15] - sout[23] - sout[39] + sout[47] - sout[63]
   + xorout[30] - xorout[31] >= -4
 R1354: - sout[7] + sout[15] - sout[23] - sout[39] + sout[47] - sout[63]
   + xorout[30] - xorout[31] >= -4
 R1355: - sout[7] - sout[15] + sout[23] - sout[39] + sout[47] - sout[63]
   + xorout[30] - xorout[31] >= -4
 R1356: - sout[7] - sout[15] - sout[23] + sout[39] + sout[47] - sout[63]
   + xorout[30] - xorout[31] >= -4
 R1357: sout[7] - sout[15] - sout[23] - sout[39] - sout[47] + sout[63]
   + xorout[30] - xorout[31] >= -4
 R1358: - sout[7] + sout[15] - sout[23] - sout[39] - sout[47] + sout[63]
   + xorout[30] - xorout[31] >= -4
 R1359: - sout[7] - sout[15] + sout[23] - sout[39] - sout[47] + sout[63]
   + xorout[30] - xorout[31] >= -4
 R1360: - sout[7] - sout[15] - sout[23] + sout[39] - sout[47] + sout[63]
   + xorout[30] - xorout[31] >= -4
 R1361: - sout[7] - sout[15] - sout[23] - sout[39] + sout[47] + sout[63]
   + xorout[30] - xorout[31] >= -4
 R1362: sout[7] + sout[15] + sout[23] + sout[39] - sout[47] - sout[63]
   + xorout[30] - xorout[31] >= -2
 R1363: sout[7] + sout[15] + sout[23] - sout[39] + sout[47] - sout[63]
   + xorout[30] - xorout[31] >= -2
 R1364: sout[7] + sout[15] - sout[23] + sout[39] + sout[47] - sout[63]
   + xorout[30] - xorout[31] >= -2
 R1365: sout[7] - sout[15] + sout[23] + sout[39] + sout[47] - sout[63]
   + xorout[30] - xorout[31] >= -2
 R1366: - sout[7] + sout[15] + sout[23] + sout[39] + sout[47] - sout[63]
   + xorout[30] - xorout[31] >= -2
 R1367: sout[7] + sout[15] + sout[23] - sout[39] - sout[47] + sout[63]
   + xorout[30] - xorout[31] >= -2
 R1368: sout[7] + sout[15] - sout[23] + sout[39] - sout[47] + sout[63]
   + xorout[30] - xorout[31] >= -2
 R1369: sout[7] - sout[15] + sout[23] + sout[39] - sout[47] + sout[63]
   + xorout[30] - xorout[31] >= -2
 R1370: - sout[7] + sout[15] + sout[23] + sout[39] - sout[47] + sout[63]
   + xorout[30] - xorout[31] >= -2
 R1371: sout[7] + sout[15] - sout[23] - sout[39] + sout[47] + sout[63]
   + xorout[30] - xorout[31] >= -2
 R1372: sout[7] - sout[15] + sout[23] - sout[39] + sout[47] + sout[63]
   + xorout[30] - xorout[31] >= -2
 R1373: - sout[7] + sout[15] + sout[23] - sout[39] + sout[47] + sout[63]
   + xorout[30] - xorout[31] >= -2
 R1374: sout[7] - sout[15] - sout[23] + sout[39] + sout[47] + sout[63]
   + xorout[30] - xorout[31] >= -2
 R1375: - sout[7] + sout[15] - sout[23] + sout[39] + sout[47] + sout[63]
   + xorout[30] - xorout[31] >= -2
 R1376: - sout[7] - sout[15] + sout[23] + sout[39] + sout[47] + sout[63]
   + xorout[30] - xorout[31] >= -2
 R1377: sout[1] + sout[9] + sout[25] + sout[33] - sout[41] + xorout[33]
   >= 0
 R1378: sout[1] + sout[9] + sout[25] - sout[33] + sout[41] + xorout[33]
   >= 0
 R1379: sout[1] + sout[9] - sout[25] + sout[33] + sout[41] + xorout[33]
   >= 0
 R1380: sout[1] - sout[9] + sout[25] + sout[33] + sout[41] + xorout[33]
   >= 0
 R1381: - sout[1] + sout[9] + sout[25] + sout[33] + sout[41] + xorout[33]
   >= 0
 R1382: sout[1] + sout[9] + sout[25] + sout[33] + sout[41] - xorout[33]
   >= 0
 R1383: sout[1] + sout[9] - sout[25] - sout[33] - sout[41] + xorout[33]
   >= -2
 R1384: sout[1] - sout[9] + sout[25] - sout[33] - sout[41] + xorout[33]
   >= -2
 R1385: - sout[1] + sout[9] + sout[25] - sout[33] - sout[41] + xorout[33]
   >= -2
 R1386: sout[1] - sout[9] - sout[25] + sout[33] - sout[41] + xorout[33]
   >= -2
 R1387: - sout[1] + sout[9] - sout[25] + sout[33] - sout[41] + xorout[33]
   >= -2
 R1388: - sout[1] - sout[9] + sout[25] + sout[33] - sout[41] + xorout[33]
   >= -2
 R1389: sout[1] - sout[9] - sout[25] - sout[33] + sout[41] + xorout[33]
   >= -2
 R1390: - sout[1] + sout[9] - sout[25] - sout[33] + sout[41] + xorout[33]
   >= -2
 R1391: - sout[1] - sout[9] + sout[25] - sout[33] + sout[41] + xorout[33]
   >= -2
 R1392: - sout[1] - sout[9] - sout[25] + sout[33] + sout[41] + xorout[33]
   >= -2
 R1393: - sout[1] - sout[9] - sout[25] - sout[33] - sout[41] + xorout[33]
   >= -4
 R1394: - sout[40] + xorout[32] >= 0
 R1395: - sout[32] + xorout[32] >= 0
 R1396: - sout[24] + xorout[32] >= 0
 R1397: - sout[8] + xorout[32] >= 0
 R1398: - sout[0] + xorout[32] >= 0
 R1399: - sout[0] + sout[1] >= 0
 R1400: - sout[8] + sout[9] >= 0
 R1401: - sout[24] + sout[25] >= 0
 R1402: - sout[32] + sout[33] >= 0
 R1403: - sout[40] + sout[41] >= 0
 R1404: sout[0] + sout[8] + sout[24] + sout[32] + sout[40] - xorout[32]
   >= 0
 R1405: - xorout[32] + xorout[33] >= 0
 R1406: sout[1] - sout[9] - sout[25] - sout[33] - sout[41] + xorout[32]
   - xorout[33] >= -4
 R1407: - sout[1] + sout[9] - sout[25] - sout[33] - sout[41] + xorout[32]
   - xorout[33] >= -4
 R1408: - sout[1] - sout[9] + sout[25] - sout[33] - sout[41] + xorout[32]
   - xorout[33] >= -4
 R1409: - sout[1] - sout[9] - sout[25] + sout[33] - sout[41] + xorout[32]
   - xorout[33] >= -4
 R1410: - sout[1] - sout[9] - sout[25] - sout[33] + sout[41] + xorout[32]
   - xorout[33] >= -4
 R1411: sout[1] + sout[9] + sout[25] - sout[33] - sout[41] + xorout[32]
   - xorout[33] >= -2
 R1412: sout[1] + sout[9] - sout[25] + sout[33] - sout[41] + xorout[32]
   - xorout[33] >= -2
 R1413: sout[1] - sout[9] + sout[25] + sout[33] - sout[41] + xorout[32]
   - xorout[33] >= -2
 R1414: - sout[1] + sout[9] + sout[25] + sout[33] - sout[41] + xorout[32]
   - xorout[33] >= -2
 R1415: sout[1] + sout[9] - sout[25] - sout[33] + sout[41] + xorout[32]
   - xorout[33] >= -2
 R1416: sout[1] - sout[9] + sout[25] - sout[33] + sout[41] + xorout[32]
   - xorout[33] >= -2
 R1417: - sout[1] + sout[9] + sout[25] - sout[33] + sout[41] + xorout[32]
   - xorout[33] >= -2
 R1418: sout[1] - sout[9] - sout[25] + sout[33] + sout[41] + xorout[32]
   - xorout[33] >= -2
 R1419: - sout[1] + sout[9] - sout[25] + sout[33] + sout[41] + xorout[32]
   - xorout[33] >= -2
 R1420: - sout[1] - sout[9] + sout[25] + sout[33] + sout[41] + xorout[32]
   - xorout[33] >= -2
 R1421: sout[3] + sout[11] + sout[27] + sout[35] - sout[43] + xorout[35]
   >= 0
 R1422: sout[3] + sout[11] + sout[27] - sout[35] + sout[43] + xorout[35]
   >= 0
 R1423: sout[3] + sout[11] - sout[27] + sout[35] + sout[43] + xorout[35]
   >= 0
 R1424: sout[3] - sout[11] + sout[27] + sout[35] + sout[43] + xorout[35]
   >= 0
 R1425: - sout[3] + sout[11] + sout[27] + sout[35] + sout[43] + xorout[35]
   >= 0
 R1426: sout[3] + sout[11] + sout[27] + sout[35] + sout[43] - xorout[35]
   >= 0
 R1427: sout[3] + sout[11] - sout[27] - sout[35] - sout[43] + xorout[35]
   >= -2
 R1428: sout[3] - sout[11] + sout[27] - sout[35] - sout[43] + xorout[35]
   >= -2
 R1429: - sout[3] + sout[11] + sout[27] - sout[35] - sout[43] + xorout[35]
   >= -2
 R1430: sout[3] - sout[11] - sout[27] + sout[35] - sout[43] + xorout[35]
   >= -2
 R1431: - sout[3] + sout[11] - sout[27] + sout[35] - sout[43] + xorout[35]
   >= -2
 R1432: - sout[3] - sout[11] + sout[27] + sout[35] - sout[43] + xorout[35]
   >= -2
 R1433: sout[3] - sout[11] - sout[27] - sout[35] + sout[43] + xorout[35]
   >= -2
 R1434: - sout[3] + sout[11] - sout[27] - sout[35] + sout[43] + xorout[35]
   >= -2
 R1435: - sout[3] - sout[11] + sout[27] - sout[35] + sout[43] + xorout[35]
   >= -2
 R1436: - sout[3] - sout[11] - sout[27] + sout[35] + sout[43] + xorout[35]
   >= -2
 R1437: - sout[3] - sout[11] - sout[27] - sout[35] - sout[43] + xorout[35]
   >= -4
 R1438: - sout[42] + xorout[34] >= 0
 R1439: - sout[34] + xorout[34] >= 0
 R1440: - sout[26] + xorout[34] >= 0
 R1441: - sout[10] + xorout[34] >= 0
 R1442: - sout[2] + xorout[34] >= 0
 R1443: - sout[2] + sout[3] >= 0
 R1444: - sout[10] + sout[11] >= 0
 R1445: - sout[26] + sout[27] >= 0
 R1446: - sout[34] + sout[35] >= 0
 R1447: - sout[42] + sout[43] >= 0
 R1448: sout[2] + sout[10] + sout[26] + sout[34] + sout[42] - xorout[34]
   >= 0
 R1449: - xorout[34] + xorout[35] >= 0
 R1450: sout[3] - sout[11] - sout[27] - sout[35] - sout[43] + xorout[34]
   - xorout[35] >= -4
 R1451: - sout[3] + sout[11] - sout[27] - sout[35] - sout[43] + xorout[34]
   - xorout[35] >= -4
 R1452: - sout[3] - sout[11] + sout[27] - sout[35] - sout[43] + xorout[34]
   - xorout[35] >= -4
 R1453: - sout[3] - sout[11] - sout[27] + sout[35] - sout[43] + xorout[34]
   - xorout[35] >= -4
 R1454: - sout[3] - sout[11] - sout[27] - sout[35] + sout[43] + xorout[34]
   - xorout[35] >= -4
 R1455: sout[3] + sout[11] + sout[27] - sout[35] - sout[43] + xorout[34]
   - xorout[35] >= -2
 R1456: sout[3] + sout[11] - sout[27] + sout[35] - sout[43] + xorout[34]
   - xorout[35] >= -2
 R1457: sout[3] - sout[11] + sout[27] + sout[35] - sout[43] + xorout[34]
   - xorout[35] >= -2
 R1458: - sout[3] + sout[11] + sout[27] + sout[35] - sout[43] + xorout[34]
   - xorout[35] >= -2
 R1459: sout[3] + sout[11] - sout[27] - sout[35] + sout[43] + xorout[34]
   - xorout[35] >= -2
 R1460: sout[3] - sout[11] + sout[27] - sout[35] + sout[43] + xorout[34]
   - xorout[35] >= -2
 R1461: - sout[3] + sout[11] + sout[27] - sout[35] + sout[43] + xorout[34]
   - xorout[35] >= -2
 R1462: sout[3] - sout[11] - sout[27] + sout[35] + sout[43] + xorout[34]
   - xorout[35] >= -2
 R1463: - sout[3] + sout[11] - sout[27] + sout[35] + sout[43] + xorout[34]
   - xorout[35] >= -2
 R1464: - sout[3] - sout[11] + sout[27] + sout[35] + sout[43] + xorout[34]
   - xorout[35] >= -2
 R1465: sout[5] + sout[13] + sout[29] + sout[37] - sout[45] + xorout[37]
   >= 0
 R1466: sout[5] + sout[13] + sout[29] - sout[37] + sout[45] + xorout[37]
   >= 0
 R1467: sout[5] + sout[13] - sout[29] + sout[37] + sout[45] + xorout[37]
   >= 0
 R1468: sout[5] - sout[13] + sout[29] + sout[37] + sout[45] + xorout[37]
   >= 0
 R1469: - sout[5] + sout[13] + sout[29] + sout[37] + sout[45] + xorout[37]
   >= 0
 R1470: sout[5] + sout[13] + sout[29] + sout[37] + sout[45] - xorout[37]
   >= 0
 R1471: sout[5] + sout[13] - sout[29] - sout[37] - sout[45] + xorout[37]
   >= -2
 R1472: sout[5] - sout[13] + sout[29] - sout[37] - sout[45] + xorout[37]
   >= -2
 R1473: - sout[5] + sout[13] + sout[29] - sout[37] - sout[45] + xorout[37]
   >= -2
 R1474: sout[5] - sout[13] - sout[29] + sout[37] - sout[45] + xorout[37]
   >= -2
 R1475: - sout[5] + sout[13] - sout[29] + sout[37] - sout[45] + xorout[37]
   >= -2
 R1476: - sout[5] - sout[13] + sout[29] + sout[37] - sout[45] + xorout[37]
   >= -2
 R1477: sout[5] - sout[13] - sout[29] - sout[37] + sout[45] + xorout[37]
   >= -2
 R1478: - sout[5] + sout[13] - sout[29] - sout[37] + sout[45] + xorout[37]
   >= -2
 R1479: - sout[5] - sout[13] + sout[29] - sout[37] + sout[45] + xorout[37]
   >= -2
 R1480: - sout[5] - sout[13] - sout[29] + sout[37] + sout[45] + xorout[37]
   >= -2
 R1481: - sout[5] - sout[13] - sout[29] - sout[37] - sout[45] + xorout[37]
   >= -4
 R1482: - sout[44] + xorout[36] >= 0
 R1483: - sout[36] + xorout[36] >= 0
 R1484: - sout[28] + xorout[36] >= 0
 R1485: - sout[12] + xorout[36] >= 0
 R1486: - sout[4] + xorout[36] >= 0
 R1487: - sout[4] + sout[5] >= 0
 R1488: - sout[12] + sout[13] >= 0
 R1489: - sout[28] + sout[29] >= 0
 R1490: - sout[36] + sout[37] >= 0
 R1491: - sout[44] + sout[45] >= 0
 R1492: sout[4] + sout[12] + sout[28] + sout[36] + sout[44] - xorout[36]
   >= 0
 R1493: - xorout[36] + xorout[37] >= 0
 R1494: sout[5] - sout[13] - sout[29] - sout[37] - sout[45] + xorout[36]
   - xorout[37] >= -4
 R1495: - sout[5] + sout[13] - sout[29] - sout[37] - sout[45] + xorout[36]
   - xorout[37] >= -4
 R1496: - sout[5] - sout[13] + sout[29] - sout[37] - sout[45] + xorout[36]
   - xorout[37] >= -4
 R1497: - sout[5] - sout[13] - sout[29] + sout[37] - sout[45] + xorout[36]
   - xorout[37] >= -4
 R1498: - sout[5] - sout[13] - sout[29] - sout[37] + sout[45] + xorout[36]
   - xorout[37] >= -4
 R1499: sout[5] + sout[13] + sout[29] - sout[37] - sout[45] + xorout[36]
   - xorout[37] >= -2
 R1500: sout[5] + sout[13] - sout[29] + sout[37] - sout[45] + xorout[36]
   - xorout[37] >= -2
 R1501: sout[5] - sout[13] + sout[29] + sout[37] - sout[45] + xorout[36]
   - xorout[37] >= -2
 R1502: - sout[5] + sout[13] + sout[29] + sout[37] - sout[45] + xorout[36]
   - xorout[37] >= -2
 R1503: sout[5] + sout[13] - sout[29] - sout[37] + sout[45] + xorout[36]
   - xorout[37] >= -2
 R1504: sout[5] - sout[13] + sout[29] - sout[37] + sout[45] + xorout[36]
   - xorout[37] >= -2
 R1505: - sout[5] + sout[13] + sout[29] - sout[37] + sout[45] + xorout[36]
   - xorout[37] >= -2
 R1506: sout[5] - sout[13] - sout[29] + sout[37] + sout[45] + xorout[36]
   - xorout[37] >= -2
 R1507: - sout[5] + sout[13] - sout[29] + sout[37] + sout[45] + xorout[36]
   - xorout[37] >= -2
 R1508: - sout[5] - sout[13] + sout[29] + sout[37] + sout[45] + xorout[36]
   - xorout[37] >= -2
 R1509: sout[7] + sout[15] + sout[31] + sout[39] - sout[47] + xorout[39]
   >= 0
 R1510: sout[7] + sout[15] + sout[31] - sout[39] + sout[47] + xorout[39]
   >= 0
 R1511: sout[7] + sout[15] - sout[31] + sout[39] + sout[47] + xorout[39]
   >= 0
 R1512: sout[7] - sout[15] + sout[31] + sout[39] + sout[47] + xorout[39]
   >= 0
 R1513: - sout[7] + sout[15] + sout[31] + sout[39] + sout[47] + xorout[39]
   >= 0
 R1514: sout[7] + sout[15] + sout[31] + sout[39] + sout[47] - xorout[39]
   >= 0
 R1515: sout[7] + sout[15] - sout[31] - sout[39] - sout[47] + xorout[39]
   >= -2
 R1516: sout[7] - sout[15] + sout[31] - sout[39] - sout[47] + xorout[39]
   >= -2
 R1517: - sout[7] + sout[15] + sout[31] - sout[39] - sout[47] + xorout[39]
   >= -2
 R1518: sout[7] - sout[15] - sout[31] + sout[39] - sout[47] + xorout[39]
   >= -2
 R1519: - sout[7] + sout[15] - sout[31] + sout[39] - sout[47] + xorout[39]
   >= -2
 R1520: - sout[7] - sout[15] + sout[31] + sout[39] - sout[47] + xorout[39]
   >= -2
 R1521: sout[7] - sout[15] - sout[31] - sout[39] + sout[47] + xorout[39]
   >= -2
 R1522: - sout[7] + sout[15] - sout[31] - sout[39] + sout[47] + xorout[39]
   >= -2
 R1523: - sout[7] - sout[15] + sout[31] - sout[39] + sout[47] + xorout[39]
   >= -2
 R1524: - sout[7] - sout[15] - sout[31] + sout[39] + sout[47] + xorout[39]
   >= -2
 R1525: - sout[7] - sout[15] - sout[31] - sout[39] - sout[47] + xorout[39]
   >= -4
 R1526: - sout[46] + xorout[38] >= 0
 R1527: - sout[38] + xorout[38] >= 0
 R1528: - sout[30] + xorout[38] >= 0
 R1529: - sout[14] + xorout[38] >= 0
 R1530: - sout[6] + xorout[38] >= 0
 R1531: - sout[6] + sout[7] >= 0
 R1532: - sout[14] + sout[15] >= 0
 R1533: - sout[30] + sout[31] >= 0
 R1534: - sout[38] + sout[39] >= 0
 R1535: - sout[46] + sout[47] >= 0
 R1536: sout[6] + sout[14] + sout[30] + sout[38] + sout[46] - xorout[38]
   >= 0
 R1537: - xorout[38] + xorout[39] >= 0
 R1538: sout[7] - sout[15] - sout[31] - sout[39] - sout[47] + xorout[38]
   - xorout[39] >= -4
 R1539: - sout[7] + sout[15] - sout[31] - sout[39] - sout[47] + xorout[38]
   - xorout[39] >= -4
 R1540: - sout[7] - sout[15] + sout[31] - sout[39] - sout[47] + xorout[38]
   - xorout[39] >= -4
 R1541: - sout[7] - sout[15] - sout[31] + sout[39] - sout[47] + xorout[38]
   - xorout[39] >= -4
 R1542: - sout[7] - sout[15] - sout[31] - sout[39] + sout[47] + xorout[38]
   - xorout[39] >= -4
 R1543: sout[7] + sout[15] + sout[31] - sout[39] - sout[47] + xorout[38]
   - xorout[39] >= -2
 R1544: sout[7] + sout[15] - sout[31] + sout[39] - sout[47] + xorout[38]
   - xorout[39] >= -2
 R1545: sout[7] - sout[15] + sout[31] + sout[39] - sout[47] + xorout[38]
   - xorout[39] >= -2
 R1546: - sout[7] + sout[15] + sout[31] + sout[39] - sout[47] + xorout[38]
   - xorout[39] >= -2
 R1547: sout[7] + sout[15] - sout[31] - sout[39] + sout[47] + xorout[38]
   - xorout[39] >= -2
 R1548: sout[7] - sout[15] + sout[31] - sout[39] + sout[47] + xorout[38]
   - xorout[39] >= -2
 R1549: - sout[7] + sout[15] + sout[31] - sout[39] + sout[47] + xorout[38]
   - xorout[39] >= -2
 R1550: sout[7] - sout[15] - sout[31] + sout[39] + sout[47] + xorout[38]
   - xorout[39] >= -2
 R1551: - sout[7] + sout[15] - sout[31] + sout[39] + sout[47] + xorout[38]
   - xorout[39] >= -2
 R1552: - sout[7] - sout[15] + sout[31] + sout[39] + sout[47] + xorout[38]
   - xorout[39] >= -2
 R1553: sout[1] + sout[9] + sout[17] + sout[41] - sout[49] + xorout[41]
   >= 0
 R1554: sout[1] + sout[9] + sout[17] - sout[41] + sout[49] + xorout[41]
   >= 0
 R1555: sout[1] + sout[9] - sout[17] + sout[41] + sout[49] + xorout[41]
   >= 0
 R1556: sout[1] - sout[9] + sout[17] + sout[41] + sout[49] + xorout[41]
   >= 0
 R1557: - sout[1] + sout[9] + sout[17] + sout[41] + sout[49] + xorout[41]
   >= 0
 R1558: sout[1] + sout[9] + sout[17] + sout[41] + sout[49] - xorout[41]
   >= 0
 R1559: sout[1] + sout[9] - sout[17] - sout[41] - sout[49] + xorout[41]
   >= -2
 R1560: sout[1] - sout[9] + sout[17] - sout[41] - sout[49] + xorout[41]
   >= -2
 R1561: - sout[1] + sout[9] + sout[17] - sout[41] - sout[49] + xorout[41]
   >= -2
 R1562: sout[1] - sout[9] - sout[17] + sout[41] - sout[49] + xorout[41]
   >= -2
 R1563: - sout[1] + sout[9] - sout[17] + sout[41] - sout[49] + xorout[41]
   >= -2
 R1564: - sout[1] - sout[9] + sout[17] + sout[41] - sout[49] + xorout[41]
   >= -2
 R1565: sout[1] - sout[9] - sout[17] - sout[41] + sout[49] + xorout[41]
   >= -2
 R1566: - sout[1] + sout[9] - sout[17] - sout[41] + sout[49] + xorout[41]
   >= -2
 R1567: - sout[1] - sout[9] + sout[17] - sout[41] + sout[49] + xorout[41]
   >= -2
 R1568: - sout[1] - sout[9] - sout[17] + sout[41] + sout[49] + xorout[41]
   >= -2
 R1569: - sout[1] - sout[9] - sout[17] - sout[41] - sout[49] + xorout[41]
   >= -4
 R1570: - sout[48] + xorout[40] >= 0
 R1571: - sout[40] + xorout[40] >= 0
 R1572: - sout[16] + xorout[40] >= 0
 R1573: - sout[8] + xorout[40] >= 0
 R1574: - sout[0] + xorout[40] >= 0
 R1575: - sout[0] + sout[1] >= 0
 R1576: - sout[8] + sout[9] >= 0
 R1577: - sout[16] + sout[17] >= 0
 R1578: - sout[40] + sout[41] >= 0
 R1579: - sout[48] + sout[49] >= 0
 R1580: sout[0] + sout[8] + sout[16] + sout[40] + sout[48] - xorout[40]
   >= 0
 R1581: - xorout[40] + xorout[41] >= 0
 R1582: sout[1] - sout[9] - sout[17] - sout[41] - sout[49] + xorout[40]
   - xorout[41] >= -4
 R1583: - sout[1] + sout[9] - sout[17] - sout[41] - sout[49] + xorout[40]
   - xorout[41] >= -4
 R1584: - sout[1] - sout[9] + sout[17] - sout[41] - sout[49] + xorout[40]
   - xorout[41] >= -4
 R1585: - sout[1] - sout[9] - sout[17] + sout[41] - sout[49] + xorout[40]
   - xorout[41] >= -4
 R1586: - sout[1] - sout[9] - sout[17] - sout[41] + sout[49] + xorout[40]
   - xorout[41] >= -4
 R1587: sout[1] + sout[9] + sout[17] - sout[41] - sout[49] + xorout[40]
   - xorout[41] >= -2
 R1588: sout[1] + sout[9] - sout[17] + sout[41] - sout[49] + xorout[40]
   - xorout[41] >= -2
 R1589: sout[1] - sout[9] + sout[17] + sout[41] - sout[49] + xorout[40]
   - xorout[41] >= -2
 R1590: - sout[1] + sout[9] + sout[17] + sout[41] - sout[49] + xorout[40]
   - xorout[41] >= -2
 R1591: sout[1] + sout[9] - sout[17] - sout[41] + sout[49] + xorout[40]
   - xorout[41] >= -2
 R1592: sout[1] - sout[9] + sout[17] - sout[41] + sout[49] + xorout[40]
   - xorout[41] >= -2
 R1593: - sout[1] + sout[9] + sout[17] - sout[41] + sout[49] + xorout[40]
   - xorout[41] >= -2
 R1594: sout[1] - sout[9] - sout[17] + sout[41] + sout[49] + xorout[40]
   - xorout[41] >= -2
 R1595: - sout[1] + sout[9] - sout[17] + sout[41] + sout[49] + xorout[40]
   - xorout[41] >= -2
 R1596: - sout[1] - sout[9] + sout[17] + sout[41] + sout[49] + xorout[40]
   - xorout[41] >= -2
 R1597: sout[3] + sout[11] + sout[19] + sout[43] - sout[51] + xorout[43]
   >= 0
 R1598: sout[3] + sout[11] + sout[19] - sout[43] + sout[51] + xorout[43]
   >= 0
 R1599: sout[3] + sout[11] - sout[19] + sout[43] + sout[51] + xorout[43]
   >= 0
 R1600: sout[3] - sout[11] + sout[19] + sout[43] + sout[51] + xorout[43]
   >= 0
 R1601: - sout[3] + sout[11] + sout[19] + sout[43] + sout[51] + xorout[43]
   >= 0
 R1602: sout[3] + sout[11] + sout[19] + sout[43] + sout[51] - xorout[43]
   >= 0
 R1603: sout[3] + sout[11] - sout[19] - sout[43] - sout[51] + xorout[43]
   >= -2
 R1604: sout[3] - sout[11] + sout[19] - sout[43] - sout[51] + xorout[43]
   >= -2
 R1605: - sout[3] + sout[11] + sout[19] - sout[43] - sout[51] + xorout[43]
   >= -2
 R1606: sout[3] - sout[11] - sout[19] + sout[43] - sout[51] + xorout[43]
   >= -2
 R1607: - sout[3] + sout[11] - sout[19] + sout[43] - sout[51] + xorout[43]
   >= -2
 R1608: - sout[3] - sout[11] + sout[19] + sout[43] - sout[51] + xorout[43]
   >= -2
 R1609: sout[3] - sout[11] - sout[19] - sout[43] + sout[51] + xorout[43]
   >= -2
 R1610: - sout[3] + sout[11] - sout[19] - sout[43] + sout[51] + xorout[43]
   >= -2
 R1611: - sout[3] - sout[11] + sout[19] - sout[43] + sout[51] + xorout[43]
   >= -2
 R1612: - sout[3] - sout[11] - sout[19] + sout[43] + sout[51] + xorout[43]
   >= -2
 R1613: - sout[3] - sout[11] - sout[19] - sout[43] - sout[51] + xorout[43]
   >= -4
 R1614: - sout[50] + xorout[42] >= 0
 R1615: - sout[42] + xorout[42] >= 0
 R1616: - sout[18] + xorout[42] >= 0
 R1617: - sout[10] + xorout[42] >= 0
 R1618: - sout[2] + xorout[42] >= 0
 R1619: - sout[2] + sout[3] >= 0
 R1620: - sout[10] + sout[11] >= 0
 R1621: - sout[18] + sout[19] >= 0
 R1622: - sout[42] + sout[43] >= 0
 R1623: - sout[50] + sout[51] >= 0
 R1624: sout[2] + sout[10] + sout[18] + sout[42] + sout[50] - xorout[42]
   >= 0
 R1625: - xorout[42] + xorout[43] >= 0
 R1626: sout[3] - sout[11] - sout[19] - sout[43] - sout[51] + xorout[42]
   - xorout[43] >= -4
 R1627: - sout[3] + sout[11] - sout[19] - sout[43] - sout[51] + xorout[42]
   - xorout[43] >= -4
 R1628: - sout[3] - sout[11] + sout[19] - sout[43] - sout[51] + xorout[42]
   - xorout[43] >= -4
 R1629: - sout[3] - sout[11] - sout[19] + sout[43] - sout[51] + xorout[42]
   - xorout[43] >= -4
 R1630: - sout[3] - sout[11] - sout[19] - sout[43] + sout[51] + xorout[42]
   - xorout[43] >= -4
 R1631: sout[3] + sout[11] + sout[19] - sout[43] - sout[51] + xorout[42]
   - xorout[43] >= -2
 R1632: sout[3] + sout[11] - sout[19] + sout[43] - sout[51] + xorout[42]
   - xorout[43] >= -2
 R1633: sout[3] - sout[11] + sout[19] + sout[43] - sout[51] + xorout[42]
   - xorout[43] >= -2
 R1634: - sout[3] + sout[11] + sout[19] + sout[43] - sout[51] + xorout[42]
   - xorout[43] >= -2
 R1635: sout[3] + sout[11] - sout[19] - sout[43] + sout[51] + xorout[42]
   - xorout[43] >= -2
 R1636: sout[3] - sout[11] + sout[19] - sout[43] + sout[51] + xorout[42]
   - xorout[43] >= -2
 R1637: - sout[3] + sout[11] + sout[19] - sout[43] + sout[51] + xorout[42]
   - xorout[43] >= -2
 R1638: sout[3] - sout[11] - sout[19] + sout[43] + sout[51] + xorout[42]
   - xorout[43] >= -2
 R1639: - sout[3] + sout[11] - sout[19] + sout[43] + sout[51] + xorout[42]
   - xorout[43] >= -2
 R1640: - sout[3] - sout[11] + sout[19] + sout[43] + sout[51] + xorout[42]
   - xorout[43] >= -2
 R1641: sout[5] + sout[13] + sout[21] + sout[45] - sout[53] + xorout[45]
   >= 0
 R1642: sout[5] + sout[13] + sout[21] - sout[45] + sout[53] + xorout[45]
   >= 0
 R1643: sout[5] + sout[13] - sout[21] + sout[45] + sout[53] + xorout[45]
   >= 0
 R1644: sout[5] - sout[13] + sout[21] + sout[45] + sout[53] + xorout[45]
   >= 0
 R1645: - sout[5] + sout[13] + sout[21] + sout[45] + sout[53] + xorout[45]
   >= 0
 R1646: sout[5] + sout[13] + sout[21] + sout[45] + sout[53] - xorout[45]
   >= 0
 R1647: sout[5] + sout[13] - sout[21] - sout[45] - sout[53] + xorout[45]
   >= -2
 R1648: sout[5] - sout[13] + sout[21] - sout[45] - sout[53] + xorout[45]
   >= -2
 R1649: - sout[5] + sout[13] + sout[21] - sout[45] - sout[53] + xorout[45]
   >= -2
 R1650: sout[5] - sout[13] - sout[21] + sout[45] - sout[53] + xorout[45]
   >= -2
 R1651: - sout[5] + sout[13] - sout[21] + sout[45] - sout[53] + xorout[45]
   >= -2
 R1652: - sout[5] - sout[13] + sout[21] + sout[45] - sout[53] + xorout[45]
   >= -2
 R1653: sout[5] - sout[13] - sout[21] - sout[45] + sout[53] + xorout[45]
   >= -2
 R1654: - sout[5] + sout[13] - sout[21] - sout[45] + sout[53] + xorout[45]
   >= -2
 R1655: - sout[5] - sout[13] + sout[21] - sout[45] + sout[53] + xorout[45]
   >= -2
 R1656: - sout[5] - sout[13] - sout[21] + sout[45] + sout[53] + xorout[45]
   >= -2
 R1657: - sout[5] - sout[13] - sout[21] - sout[45] - sout[53] + xorout[45]
   >= -4
 R1658: - sout[52] + xorout[44] >= 0
 R1659: - sout[44] + xorout[44] >= 0
 R1660: - sout[20] + xorout[44] >= 0
 R1661: - sout[12] + xorout[44] >= 0
 R1662: - sout[4] + xorout[44] >= 0
 R1663: - sout[4] + sout[5] >= 0
 R1664: - sout[12] + sout[13] >= 0
 R1665: - sout[20] + sout[21] >= 0
 R1666: - sout[44] + sout[45] >= 0
 R1667: - sout[52] + sout[53] >= 0
 R1668: sout[4] + sout[12] + sout[20] + sout[44] + sout[52] - xorout[44]
   >= 0
 R1669: - xorout[44] + xorout[45] >= 0
 R1670: sout[5] - sout[13] - sout[21] - sout[45] - sout[53] + xorout[44]
   - xorout[45] >= -4
 R1671: - sout[5] + sout[13] - sout[21] - sout[45] - sout[53] + xorout[44]
   - xorout[45] >= -4
 R1672: - sout[5] - sout[13] + sout[21] - sout[45] - sout[53] + xorout[44]
   - xorout[45] >= -4
 R1673: - sout[5] - sout[13] - sout[21] + sout[45] - sout[53] + xorout[44]
   - xorout[45] >= -4
 R1674: - sout[5] - sout[13] - sout[21] - sout[45] + sout[53] + xorout[44]
   - xorout[45] >= -4
 R1675: sout[5] + sout[13] + sout[21] - sout[45] - sout[53] + xorout[44]
   - xorout[45] >= -2
 R1676: sout[5] + sout[13] - sout[21] + sout[45] - sout[53] + xorout[44]
   - xorout[45] >= -2
 R1677: sout[5] - sout[13] + sout[21] + sout[45] - sout[53] + xorout[44]
   - xorout[45] >= -2
 R1678: - sout[5] + sout[13] + sout[21] + sout[45] - sout[53] + xorout[44]
   - xorout[45] >= -2
 R1679: sout[5] + sout[13] - sout[21] - sout[45] + sout[53] + xorout[44]
   - xorout[45] >= -2
 R1680: sout[5] - sout[13] + sout[21] - sout[45] + sout[53] + xorout[44]
   - xorout[45] >= -2
 R1681: - sout[5] + sout[13] + sout[21] - sout[45] + sout[53] + xorout[44]
   - xorout[45] >= -2
 R1682: sout[5] - sout[13] - sout[21] + sout[45] + sout[53] + xorout[44]
   - xorout[45] >= -2
 R1683: - sout[5] + sout[13] - sout[21] + sout[45] + sout[53] + xorout[44]
   - xorout[45] >= -2
 R1684: - sout[5] - sout[13] + sout[21] + sout[45] + sout[53] + xorout[44]
   - xorout[45] >= -2
 R1685: sout[7] + sout[15] + sout[23] + sout[47] - sout[55] + xorout[47]
   >= 0
 R1686: sout[7] + sout[15] + sout[23] - sout[47] + sout[55] + xorout[47]
   >= 0
 R1687: sout[7] + sout[15] - sout[23] + sout[47] + sout[55] + xorout[47]
   >= 0
 R1688: sout[7] - sout[15] + sout[23] + sout[47] + sout[55] + xorout[47]
   >= 0
 R1689: - sout[7] + sout[15] + sout[23] + sout[47] + sout[55] + xorout[47]
   >= 0
 R1690: sout[7] + sout[15] + sout[23] + sout[47] + sout[55] - xorout[47]
   >= 0
 R1691: sout[7] + sout[15] - sout[23] - sout[47] - sout[55] + xorout[47]
   >= -2
 R1692: sout[7] - sout[15] + sout[23] - sout[47] - sout[55] + xorout[47]
   >= -2
 R1693: - sout[7] + sout[15] + sout[23] - sout[47] - sout[55] + xorout[47]
   >= -2
 R1694: sout[7] - sout[15] - sout[23] + sout[47] - sout[55] + xorout[47]
   >= -2
 R1695: - sout[7] + sout[15] - sout[23] + sout[47] - sout[55] + xorout[47]
   >= -2
 R1696: - sout[7] - sout[15] + sout[23] + sout[47] - sout[55] + xorout[47]
   >= -2
 R1697: sout[7] - sout[15] - sout[23] - sout[47] + sout[55] + xorout[47]
   >= -2
 R1698: - sout[7] + sout[15] - sout[23] - sout[47] + sout[55] + xorout[47]
   >= -2
 R1699: - sout[7] - sout[15] + sout[23] - sout[47] + sout[55] + xorout[47]
   >= -2
 R1700: - sout[7] - sout[15] - sout[23] + sout[47] + sout[55] + xorout[47]
   >= -2
 R1701: - sout[7] - sout[15] - sout[23] - sout[47] - sout[55] + xorout[47]
   >= -4
 R1702: - sout[54] + xorout[46] >= 0
 R1703: - sout[46] + xorout[46] >= 0
 R1704: - sout[22] + xorout[46] >= 0
 R1705: - sout[14] + xorout[46] >= 0
 R1706: - sout[6] + xorout[46] >= 0
 R1707: - sout[6] + sout[7] >= 0
 R1708: - sout[14] + sout[15] >= 0
 R1709: - sout[22] + sout[23] >= 0
 R1710: - sout[46] + sout[47] >= 0
 R1711: - sout[54] + sout[55] >= 0
 R1712: sout[6] + sout[14] + sout[22] + sout[46] + sout[54] - xorout[46]
   >= 0
 R1713: - xorout[46] + xorout[47] >= 0
 R1714: sout[7] - sout[15] - sout[23] - sout[47] - sout[55] + xorout[46]
   - xorout[47] >= -4
 R1715: - sout[7] + sout[15] - sout[23] - sout[47] - sout[55] + xorout[46]
   - xorout[47] >= -4
 R1716: - sout[7] - sout[15] + sout[23] - sout[47] - sout[55] + xorout[46]
   - xorout[47] >= -4
 R1717: - sout[7] - sout[15] - sout[23] + sout[47] - sout[55] + xorout[46]
   - xorout[47] >= -4
 R1718: - sout[7] - sout[15] - sout[23] - sout[47] + sout[55] + xorout[46]
   - xorout[47] >= -4
 R1719: sout[7] + sout[15] + sout[23] - sout[47] - sout[55] + xorout[46]
   - xorout[47] >= -2
 R1720: sout[7] + sout[15] - sout[23] + sout[47] - sout[55] + xorout[46]
   - xorout[47] >= -2
 R1721: sout[7] - sout[15] + sout[23] + sout[47] - sout[55] + xorout[46]
   - xorout[47] >= -2
 R1722: - sout[7] + sout[15] + sout[23] + sout[47] - sout[55] + xorout[46]
   - xorout[47] >= -2
 R1723: sout[7] + sout[15] - sout[23] - sout[47] + sout[55] + xorout[46]
   - xorout[47] >= -2
 R1724: sout[7] - sout[15] + sout[23] - sout[47] + sout[55] + xorout[46]
   - xorout[47] >= -2
 R1725: - sout[7] + sout[15] + sout[23] - sout[47] + sout[55] + xorout[46]
   - xorout[47] >= -2
 R1726: sout[7] - sout[15] - sout[23] + sout[47] + sout[55] + xorout[46]
   - xorout[47] >= -2
 R1727: - sout[7] + sout[15] - sout[23] + sout[47] + sout[55] + xorout[46]
   - xorout[47] >= -2
 R1728: - sout[7] - sout[15] + sout[23] + sout[47] + sout[55] + xorout[46]
   - xorout[47] >= -2
 R1729: sout[9] + sout[17] + sout[25] + sout[49] - sout[57] + xorout[49]
   >= 0
 R1730: sout[9] + sout[17] + sout[25] - sout[49] + sout[57] + xorout[49]
   >= 0
 R1731: sout[9] + sout[17] - sout[25] + sout[49] + sout[57] + xorout[49]
   >= 0
 R1732: sout[9] - sout[17] + sout[25] + sout[49] + sout[57] + xorout[49]
   >= 0
 R1733: - sout[9] + sout[17] + sout[25] + sout[49] + sout[57] + xorout[49]
   >= 0
 R1734: sout[9] + sout[17] + sout[25] + sout[49] + sout[57] - xorout[49]
   >= 0
 R1735: sout[9] + sout[17] - sout[25] - sout[49] - sout[57] + xorout[49]
   >= -2
 R1736: sout[9] - sout[17] + sout[25] - sout[49] - sout[57] + xorout[49]
   >= -2
 R1737: - sout[9] + sout[17] + sout[25] - sout[49] - sout[57] + xorout[49]
   >= -2
 R1738: sout[9] - sout[17] - sout[25] + sout[49] - sout[57] + xorout[49]
   >= -2
 R1739: - sout[9] + sout[17] - sout[25] + sout[49] - sout[57] + xorout[49]
   >= -2
 R1740: - sout[9] - sout[17] + sout[25] + sout[49] - sout[57] + xorout[49]
   >= -2
 R1741: sout[9] - sout[17] - sout[25] - sout[49] + sout[57] + xorout[49]
   >= -2
 R1742: - sout[9] + sout[17] - sout[25] - sout[49] + sout[57] + xorout[49]
   >= -2
 R1743: - sout[9] - sout[17] + sout[25] - sout[49] + sout[57] + xorout[49]
   >= -2
 R1744: - sout[9] - sout[17] - sout[25] + sout[49] + sout[57] + xorout[49]
   >= -2
 R1745: - sout[9] - sout[17] - sout[25] - sout[49] - sout[57] + xorout[49]
   >= -4
 R1746: - sout[56] + xorout[48] >= 0
 R1747: - sout[48] + xorout[48] >= 0
 R1748: - sout[24] + xorout[48] >= 0
 R1749: - sout[16] + xorout[48] >= 0
 R1750: - sout[8] + xorout[48] >= 0
 R1751: - sout[8] + sout[9] >= 0
 R1752: - sout[16] + sout[17] >= 0
 R1753: - sout[24] + sout[25] >= 0
 R1754: - sout[48] + sout[49] >= 0
 R1755: - sout[56] + sout[57] >= 0
 R1756: sout[8] + sout[16] + sout[24] + sout[48] + sout[56] - xorout[48]
   >= 0
 R1757: - xorout[48] + xorout[49] >= 0
 R1758: sout[9] - sout[17] - sout[25] - sout[49] - sout[57] + xorout[48]
   - xorout[49] >= -4
 R1759: - sout[9] + sout[17] - sout[25] - sout[49] - sout[57] + xorout[48]
   - xorout[49] >= -4
 R1760: - sout[9] - sout[17] + sout[25] - sout[49] - sout[57] + xorout[48]
   - xorout[49] >= -4
 R1761: - sout[9] - sout[17] - sout[25] + sout[49] - sout[57] + xorout[48]
   - xorout[49] >= -4
 R1762: - sout[9] - sout[17] - sout[25] - sout[49] + sout[57] + xorout[48]
   - xorout[49] >= -4
 R1763: sout[9] + sout[17] + sout[25] - sout[49] - sout[57] + xorout[48]
   - xorout[49] >= -2
 R1764: sout[9] + sout[17] - sout[25] + sout[49] - sout[57] + xorout[48]
   - xorout[49] >= -2
 R1765: sout[9] - sout[17] + sout[25] + sout[49] - sout[57] + xorout[48]
   - xorout[49] >= -2
 R1766: - sout[9] + sout[17] + sout[25] + sout[49] - sout[57] + xorout[48]
   - xorout[49] >= -2
 R1767: sout[9] + sout[17] - sout[25] - sout[49] + sout[57] + xorout[48]
   - xorout[49] >= -2
 R1768: sout[9] - sout[17] + sout[25] - sout[49] + sout[57] + xorout[48]
   - xorout[49] >= -2
 R1769: - sout[9] + sout[17] + sout[25] - sout[49] + sout[57] + xorout[48]
   - xorout[49] >= -2
 R1770: sout[9] - sout[17] - sout[25] + sout[49] + sout[57] + xorout[48]
   - xorout[49] >= -2
 R1771: - sout[9] + sout[17] - sout[25] + sout[49] + sout[57] + xorout[48]
   - xorout[49] >= -2
 R1772: - sout[9] - sout[17] + sout[25] + sout[49] + sout[57] + xorout[48]
   - xorout[49] >= -2
 R1773: sout[11] + sout[19] + sout[27] + sout[51] - sout[59] + xorout[51]
   >= 0
 R1774: sout[11] + sout[19] + sout[27] - sout[51] + sout[59] + xorout[51]
   >= 0
 R1775: sout[11] + sout[19] - sout[27] + sout[51] + sout[59] + xorout[51]
   >= 0
 R1776: sout[11] - sout[19] + sout[27] + sout[51] + sout[59] + xorout[51]
   >= 0
 R1777: - sout[11] + sout[19] + sout[27] + sout[51] + sout[59] + xorout[51]
   >= 0
 R1778: sout[11] + sout[19] + sout[27] + sout[51] + sout[59] - xorout[51]
   >= 0
 R1779: sout[11] + sout[19] - sout[27] - sout[51] - sout[59] + xorout[51]
   >= -2
 R1780: sout[11] - sout[19] + sout[27] - sout[51] - sout[59] + xorout[51]
   >= -2
 R1781: - sout[11] + sout[19] + sout[27] - sout[51] - sout[59] + xorout[51]
   >= -2
 R1782: sout[11] - sout[19] - sout[27] + sout[51] - sout[59] + xorout[51]
   >= -2
 R1783: - sout[11] + sout[19] - sout[27] + sout[51] - sout[59] + xorout[51]
   >= -2
 R1784: - sout[11] - sout[19] + sout[27] + sout[51] - sout[59] + xorout[51]
   >= -2
 R1785: sout[11] - sout[19] - sout[27] - sout[51] + sout[59] + xorout[51]
   >= -2
 R1786: - sout[11] + sout[19] - sout[27] - sout[51] + sout[59] + xorout[51]
   >= -2
 R1787: - sout[11] - sout[19] + sout[27] - sout[51] + sout[59] + xorout[51]
   >= -2
 R1788: - sout[11] - sout[19] - sout[27] + sout[51] + sout[59] + xorout[51]
   >= -2
 R1789: - sout[11] - sout[19] - sout[27] - sout[51] - sout[59] + xorout[51]
   >= -4
 R1790: - sout[58] + xorout[50] >= 0
 R1791: - sout[50] + xorout[50] >= 0
 R1792: - sout[26] + xorout[50] >= 0
 R1793: - sout[18] + xorout[50] >= 0
 R1794: - sout[10] + xorout[50] >= 0
 R1795: - sout[10] + sout[11] >= 0
 R1796: - sout[18] + sout[19] >= 0
 R1797: - sout[26] + sout[27] >= 0
 R1798: - sout[50] + sout[51] >= 0
 R1799: - sout[58] + sout[59] >= 0
 R1800: sout[10] + sout[18] + sout[26] + sout[50] + sout[58] - xorout[50]
   >= 0
 R1801: - xorout[50] + xorout[51] >= 0
 R1802: sout[11] - sout[19] - sout[27] - sout[51] - sout[59] + xorout[50]
   - xorout[51] >= -4
 R1803: - sout[11] + sout[19] - sout[27] - sout[51] - sout[59] + xorout[50]
   - xorout[51] >= -4
 R1804: - sout[11] - sout[19] + sout[27] - sout[51] - sout[59] + xorout[50]
   - xorout[51] >= -4
 R1805: - sout[11] - sout[19] - sout[27] + sout[51] - sout[59] + xorout[50]
   - xorout[51] >= -4
 R1806: - sout[11] - sout[19] - sout[27] - sout[51] + sout[59] + xorout[50]
   - xorout[51] >= -4
 R1807: sout[11] + sout[19] + sout[27] - sout[51] - sout[59] + xorout[50]
   - xorout[51] >= -2
 R1808: sout[11] + sout[19] - sout[27] + sout[51] - sout[59] + xorout[50]
   - xorout[51] >= -2
 R1809: sout[11] - sout[19] + sout[27] + sout[51] - sout[59] + xorout[50]
   - xorout[51] >= -2
 R1810: - sout[11] + sout[19] + sout[27] + sout[51] - sout[59] + xorout[50]
   - xorout[51] >= -2
 R1811: sout[11] + sout[19] - sout[27] - sout[51] + sout[59] + xorout[50]
   - xorout[51] >= -2
 R1812: sout[11] - sout[19] + sout[27] - sout[51] + sout[59] + xorout[50]
   - xorout[51] >= -2
 R1813: - sout[11] + sout[19] + sout[27] - sout[51] + sout[59] + xorout[50]
   - xorout[51] >= -2
 R1814: sout[11] - sout[19] - sout[27] + sout[51] + sout[59] + xorout[50]
   - xorout[51] >= -2
 R1815: - sout[11] + sout[19] - sout[27] + sout[51] + sout[59] + xorout[50]
   - xorout[51] >= -2
 R1816: - sout[11] - sout[19] + sout[27] + sout[51] + sout[59] + xorout[50]
   - xorout[51] >= -2
 R1817: sout[13] + sout[21] + sout[29] + sout[53] - sout[61] + xorout[53]
   >= 0
 R1818: sout[13] + sout[21] + sout[29] - sout[53] + sout[61] + xorout[53]
   >= 0
 R1819: sout[13] + sout[21] - sout[29] + sout[53] + sout[61] + xorout[53]
   >= 0
 R1820: sout[13] - sout[21] + sout[29] + sout[53] + sout[61] + xorout[53]
   >= 0
 R1821: - sout[13] + sout[21] + sout[29] + sout[53] + sout[61] + xorout[53]
   >= 0
 R1822: sout[13] + sout[21] + sout[29] + sout[53] + sout[61] - xorout[53]
   >= 0
 R1823: sout[13] + sout[21] - sout[29] - sout[53] - sout[61] + xorout[53]
   >= -2
 R1824: sout[13] - sout[21] + sout[29] - sout[53] - sout[61] + xorout[53]
   >= -2
 R1825: - sout[13] + sout[21] + sout[29] - sout[53] - sout[61] + xorout[53]
   >= -2
 R1826: sout[13] - sout[21] - sout[29] + sout[53] - sout[61] + xorout[53]
   >= -2
 R1827: - sout[13] + sout[21] - sout[29] + sout[53] - sout[61] + xorout[53]
   >= -2
 R1828: - sout[13] - sout[21] + sout[29] + sout[53] - sout[61] + xorout[53]
   >= -2
 R1829: sout[13] - sout[21] - sout[29] - sout[53] + sout[61] + xorout[53]
   >= -2
 R1830: - sout[13] + sout[21] - sout[29] - sout[53] + sout[61] + xorout[53]
   >= -2
 R1831: - sout[13] - sout[21] + sout[29] - sout[53] + sout[61] + xorout[53]
   >= -2
 R1832: - sout[13] - sout[21] - sout[29] + sout[53] + sout[61] + xorout[53]
   >= -2
 R1833: - sout[13] - sout[21] - sout[29] - sout[53] - sout[61] + xorout[53]
   >= -4
 R1834: - sout[60] + xorout[52] >= 0
 R1835: - sout[52] + xorout[52] >= 0
 R1836: - sout[28] + xorout[52] >= 0
 R1837: - sout[20] + xorout[52] >= 0
 R1838: - sout[12] + xorout[52] >= 0
 R1839: - sout[12] + sout[13] >= 0
 R1840: - sout[20] + sout[21] >= 0
 R1841: - sout[28] + sout[29] >= 0
 R1842: - sout[52] + sout[53] >= 0
 R1843: - sout[60] + sout[61] >= 0
 R1844: sout[12] + sout[20] + sout[28] + sout[52] + sout[60] - xorout[52]
   >= 0
 R1845: - xorout[52] + xorout[53] >= 0
 R1846: sout[13] - sout[21] - sout[29] - sout[53] - sout[61] + xorout[52]
   - xorout[53] >= -4
 R1847: - sout[13] + sout[21] - sout[29] - sout[53] - sout[61] + xorout[52]
   - xorout[53] >= -4
 R1848: - sout[13] - sout[21] + sout[29] - sout[53] - sout[61] + xorout[52]
   - xorout[53] >= -4
 R1849: - sout[13] - sout[21] - sout[29] + sout[53] - sout[61] + xorout[52]
   - xorout[53] >= -4
 R1850: - sout[13] - sout[21] - sout[29] - sout[53] + sout[61] + xorout[52]
   - xorout[53] >= -4
 R1851: sout[13] + sout[21] + sout[29] - sout[53] - sout[61] + xorout[52]
   - xorout[53] >= -2
 R1852: sout[13] + sout[21] - sout[29] + sout[53] - sout[61] + xorout[52]
   - xorout[53] >= -2
 R1853: sout[13] - sout[21] + sout[29] + sout[53] - sout[61] + xorout[52]
   - xorout[53] >= -2
 R1854: - sout[13] + sout[21] + sout[29] + sout[53] - sout[61] + xorout[52]
   - xorout[53] >= -2
 R1855: sout[13] + sout[21] - sout[29] - sout[53] + sout[61] + xorout[52]
   - xorout[53] >= -2
 R1856: sout[13] - sout[21] + sout[29] - sout[53] + sout[61] + xorout[52]
   - xorout[53] >= -2
 R1857: - sout[13] + sout[21] + sout[29] - sout[53] + sout[61] + xorout[52]
   - xorout[53] >= -2
 R1858: sout[13] - sout[21] - sout[29] + sout[53] + sout[61] + xorout[52]
   - xorout[53] >= -2
 R1859: - sout[13] + sout[21] - sout[29] + sout[53] + sout[61] + xorout[52]
   - xorout[53] >= -2
 R1860: - sout[13] - sout[21] + sout[29] + sout[53] + sout[61] + xorout[52]
   - xorout[53] >= -2
 R1861: sout[15] + sout[23] + sout[31] + sout[55] - sout[63] + xorout[55]
   >= 0
 R1862: sout[15] + sout[23] + sout[31] - sout[55] + sout[63] + xorout[55]
   >= 0
 R1863: sout[15] + sout[23] - sout[31] + sout[55] + sout[63] + xorout[55]
   >= 0
 R1864: sout[15] - sout[23] + sout[31] + sout[55] + sout[63] + xorout[55]
   >= 0
 R1865: - sout[15] + sout[23] + sout[31] + sout[55] + sout[63] + xorout[55]
   >= 0
 R1866: sout[15] + sout[23] + sout[31] + sout[55] + sout[63] - xorout[55]
   >= 0
 R1867: sout[15] + sout[23] - sout[31] - sout[55] - sout[63] + xorout[55]
   >= -2
 R1868: sout[15] - sout[23] + sout[31] - sout[55] - sout[63] + xorout[55]
   >= -2
 R1869: - sout[15] + sout[23] + sout[31] - sout[55] - sout[63] + xorout[55]
   >= -2
 R1870: sout[15] - sout[23] - sout[31] + sout[55] - sout[63] + xorout[55]
   >= -2
 R1871: - sout[15] + sout[23] - sout[31] + sout[55] - sout[63] + xorout[55]
   >= -2
 R1872: - sout[15] - sout[23] + sout[31] + sout[55] - sout[63] + xorout[55]
   >= -2
 R1873: sout[15] - sout[23] - sout[31] - sout[55] + sout[63] + xorout[55]
   >= -2
 R1874: - sout[15] + sout[23] - sout[31] - sout[55] + sout[63] + xorout[55]
   >= -2
 R1875: - sout[15] - sout[23] + sout[31] - sout[55] + sout[63] + xorout[55]
   >= -2
 R1876: - sout[15] - sout[23] - sout[31] + sout[55] + sout[63] + xorout[55]
   >= -2
 R1877: - sout[15] - sout[23] - sout[31] - sout[55] - sout[63] + xorout[55]
   >= -4
 R1878: - sout[62] + xorout[54] >= 0
 R1879: - sout[54] + xorout[54] >= 0
 R1880: - sout[30] + xorout[54] >= 0
 R1881: - sout[22] + xorout[54] >= 0
 R1882: - sout[14] + xorout[54] >= 0
 R1883: - sout[14] + sout[15] >= 0
 R1884: - sout[22] + sout[23] >= 0
 R1885: - sout[30] + sout[31] >= 0
 R1886: - sout[54] + sout[55] >= 0
 R1887: - sout[62] + sout[63] >= 0
 R1888: sout[14] + sout[22] + sout[30] + sout[54] + sout[62] - xorout[54]
   >= 0
 R1889: - xorout[54] + xorout[55] >= 0
 R1890: sout[15] - sout[23] - sout[31] - sout[55] - sout[63] + xorout[54]
   - xorout[55] >= -4
 R1891: - sout[15] + sout[23] - sout[31] - sout[55] - sout[63] + xorout[54]
   - xorout[55] >= -4
 R1892: - sout[15] - sout[23] + sout[31] - sout[55] - sout[63] + xorout[54]
   - xorout[55] >= -4
 R1893: - sout[15] - sout[23] - sout[31] + sout[55] - sout[63] + xorout[54]
   - xorout[55] >= -4
 R1894: - sout[15] - sout[23] - sout[31] - sout[55] + sout[63] + xorout[54]
   - xorout[55] >= -4
 R1895: sout[15] + sout[23] + sout[31] - sout[55] - sout[63] + xorout[54]
   - xorout[55] >= -2
 R1896: sout[15] + sout[23] - sout[31] + sout[55] - sout[63] + xorout[54]
   - xorout[55] >= -2
 R1897: sout[15] - sout[23] + sout[31] + sout[55] - sout[63] + xorout[54]
   - xorout[55] >= -2
 R1898: - sout[15] + sout[23] + sout[31] + sout[55] - sout[63] + xorout[54]
   - xorout[55] >= -2
 R1899: sout[15] + sout[23] - sout[31] - sout[55] + sout[63] + xorout[54]
   - xorout[55] >= -2
 R1900: sout[15] - sout[23] + sout[31] - sout[55] + sout[63] + xorout[54]
   - xorout[55] >= -2
 R1901: - sout[15] + sout[23] + sout[31] - sout[55] + sout[63] + xorout[54]
   - xorout[55] >= -2
 R1902: sout[15] - sout[23] - sout[31] + sout[55] + sout[63] + xorout[54]
   - xorout[55] >= -2
 R1903: - sout[15] + sout[23] - sout[31] + sout[55] + sout[63] + xorout[54]
   - xorout[55] >= -2
 R1904: - sout[15] - sout[23] + sout[31] + sout[55] + sout[63] + xorout[54]
   - xorout[55] >= -2
 R1905: sout[1] + sout[17] + sout[25] + sout[33] - sout[57] + xorout[57]
   >= 0
 R1906: sout[1] + sout[17] + sout[25] - sout[33] + sout[57] + xorout[57]
   >= 0
 R1907: sout[1] + sout[17] - sout[25] + sout[33] + sout[57] + xorout[57]
   >= 0
 R1908: sout[1] - sout[17] + sout[25] + sout[33] + sout[57] + xorout[57]
   >= 0
 R1909: - sout[1] + sout[17] + sout[25] + sout[33] + sout[57] + xorout[57]
   >= 0
 R1910: sout[1] + sout[17] + sout[25] + sout[33] + sout[57] - xorout[57]
   >= 0
 R1911: sout[1] + sout[17] - sout[25] - sout[33] - sout[57] + xorout[57]
   >= -2
 R1912: sout[1] - sout[17] + sout[25] - sout[33] - sout[57] + xorout[57]
   >= -2
 R1913: - sout[1] + sout[17] + sout[25] - sout[33] - sout[57] + xorout[57]
   >= -2
 R1914: sout[1] - sout[17] - sout[25] + sout[33] - sout[57] + xorout[57]
   >= -2
 R1915: - sout[1] + sout[17] - sout[25] + sout[33] - sout[57] + xorout[57]
   >= -2
 R1916: - sout[1] - sout[17] + sout[25] + sout[33] - sout[57] + xorout[57]
   >= -2
 R1917: sout[1] - sout[17] - sout[25] - sout[33] + sout[57] + xorout[57]
   >= -2
 R1918: - sout[1] + sout[17] - sout[25] - sout[33] + sout[57] + xorout[57]
   >= -2
 R1919: - sout[1] - sout[17] + sout[25] - sout[33] + sout[57] + xorout[57]
   >= -2
 R1920: - sout[1] - sout[17] - sout[25] + sout[33] + sout[57] + xorout[57]
   >= -2
 R1921: - sout[1] - sout[17] - sout[25] - sout[33] - sout[57] + xorout[57]
   >= -4
 R1922: - sout[56] + xorout[56] >= 0
 R1923: - sout[32] + xorout[56] >= 0
 R1924: - sout[24] + xorout[56] >= 0
 R1925: - sout[16] + xorout[56] >= 0
 R1926: - sout[0] + xorout[56] >= 0
 R1927: - sout[0] + sout[1] >= 0
 R1928: - sout[16] + sout[17] >= 0
 R1929: - sout[24] + sout[25] >= 0
 R1930: - sout[32] + sout[33] >= 0
 R1931: - sout[56] + sout[57] >= 0
 R1932: sout[0] + sout[16] + sout[24] + sout[32] + sout[56] - xorout[56]
   >= 0
 R1933: - xorout[56] + xorout[57] >= 0
 R1934: sout[1] - sout[17] - sout[25] - sout[33] - sout[57] + xorout[56]
   - xorout[57] >= -4
 R1935: - sout[1] + sout[17] - sout[25] - sout[33] - sout[57] + xorout[56]
   - xorout[57] >= -4
 R1936: - sout[1] - sout[17] + sout[25] - sout[33] - sout[57] + xorout[56]
   - xorout[57] >= -4
 R1937: - sout[1] - sout[17] - sout[25] + sout[33] - sout[57] + xorout[56]
   - xorout[57] >= -4
 R1938: - sout[1] - sout[17] - sout[25] - sout[33] + sout[57] + xorout[56]
   - xorout[57] >= -4
 R1939: sout[1] + sout[17] + sout[25] - sout[33] - sout[57] + xorout[56]
   - xorout[57] >= -2
 R1940: sout[1] + sout[17] - sout[25] + sout[33] - sout[57] + xorout[56]
   - xorout[57] >= -2
 R1941: sout[1] - sout[17] + sout[25] + sout[33] - sout[57] + xorout[56]
   - xorout[57] >= -2
 R1942: - sout[1] + sout[17] + sout[25] + sout[33] - sout[57] + xorout[56]
   - xorout[57] >= -2
 R1943: sout[1] + sout[17] - sout[25] - sout[33] + sout[57] + xorout[56]
   - xorout[57] >= -2
 R1944: sout[1] - sout[17] + sout[25] - sout[33] + sout[57] + xorout[56]
   - xorout[57] >= -2
 R1945: - sout[1] + sout[17] + sout[25] - sout[33] + sout[57] + xorout[56]
   - xorout[57] >= -2
 R1946: sout[1] - sout[17] - sout[25] + sout[33] + sout[57] + xorout[56]
   - xorout[57] >= -2
 R1947: - sout[1] + sout[17] - sout[25] + sout[33] + sout[57] + xorout[56]
   - xorout[57] >= -2
 R1948: - sout[1] - sout[17] + sout[25] + sout[33] + sout[57] + xorout[56]
   - xorout[57] >= -2
 R1949: sout[3] + sout[19] + sout[27] + sout[35] - sout[59] + xorout[59]
   >= 0
 R1950: sout[3] + sout[19] + sout[27] - sout[35] + sout[59] + xorout[59]
   >= 0
 R1951: sout[3] + sout[19] - sout[27] + sout[35] + sout[59] + xorout[59]
   >= 0
 R1952: sout[3] - sout[19] + sout[27] + sout[35] + sout[59] + xorout[59]
   >= 0
 R1953: - sout[3] + sout[19] + sout[27] + sout[35] + sout[59] + xorout[59]
   >= 0
 R1954: sout[3] + sout[19] + sout[27] + sout[35] + sout[59] - xorout[59]
   >= 0
 R1955: sout[3] + sout[19] - sout[27] - sout[35] - sout[59] + xorout[59]
   >= -2
 R1956: sout[3] - sout[19] + sout[27] - sout[35] - sout[59] + xorout[59]
   >= -2
 R1957: - sout[3] + sout[19] + sout[27] - sout[35] - sout[59] + xorout[59]
   >= -2
 R1958: sout[3] - sout[19] - sout[27] + sout[35] - sout[59] + xorout[59]
   >= -2
 R1959: - sout[3] + sout[19] - sout[27] + sout[35] - sout[59] + xorout[59]
   >= -2
 R1960: - sout[3] - sout[19] + sout[27] + sout[35] - sout[59] + xorout[59]
   >= -2
 R1961: sout[3] - sout[19] - sout[27] - sout[35] + sout[59] + xorout[59]
   >= -2
 R1962: - sout[3] + sout[19] - sout[27] - sout[35] + sout[59] + xorout[59]
   >= -2
 R1963: - sout[3] - sout[19] + sout[27] - sout[35] + sout[59] + xorout[59]
   >= -2
 R1964: - sout[3] - sout[19] - sout[27] + sout[35] + sout[59] + xorout[59]
   >= -2
 R1965: - sout[3] - sout[19] - sout[27] - sout[35] - sout[59] + xorout[59]
   >= -4
 R1966: - sout[58] + xorout[58] >= 0
 R1967: - sout[34] + xorout[58] >= 0
 R1968: - sout[26] + xorout[58] >= 0
 R1969: - sout[18] + xorout[58] >= 0
 R1970: - sout[2] + xorout[58] >= 0
 R1971: - sout[2] + sout[3] >= 0
 R1972: - sout[18] + sout[19] >= 0
 R1973: - sout[26] + sout[27] >= 0
 R1974: - sout[34] + sout[35] >= 0
 R1975: - sout[58] + sout[59] >= 0
 R1976: sout[2] + sout[18] + sout[26] + sout[34] + sout[58] - xorout[58]
   >= 0
 R1977: - xorout[58] + xorout[59] >= 0
 R1978: sout[3] - sout[19] - sout[27] - sout[35] - sout[59] + xorout[58]
   - xorout[59] >= -4
 R1979: - sout[3] + sout[19] - sout[27] - sout[35] - sout[59] + xorout[58]
   - xorout[59] >= -4
 R1980: - sout[3] - sout[19] + sout[27] - sout[35] - sout[59] + xorout[58]
   - xorout[59] >= -4
 R1981: - sout[3] - sout[19] - sout[27] + sout[35] - sout[59] + xorout[58]
   - xorout[59] >= -4
 R1982: - sout[3] - sout[19] - sout[27] - sout[35] + sout[59] + xorout[58]
   - xorout[59] >= -4
 R1983: sout[3] + sout[19] + sout[27] - sout[35] - sout[59] + xorout[58]
   - xorout[59] >= -2
 R1984: sout[3] + sout[19] - sout[27] + sout[35] - sout[59] + xorout[58]
   - xorout[59] >= -2
 R1985: sout[3] - sout[19] + sout[27] + sout[35] - sout[59] + xorout[58]
   - xorout[59] >= -2
 R1986: - sout[3] + sout[19] + sout[27] + sout[35] - sout[59] + xorout[58]
   - xorout[59] >= -2
 R1987: sout[3] + sout[19] - sout[27] - sout[35] + sout[59] + xorout[58]
   - xorout[59] >= -2
 R1988: sout[3] - sout[19] + sout[27] - sout[35] + sout[59] + xorout[58]
   - xorout[59] >= -2
 R1989: - sout[3] + sout[19] + sout[27] - sout[35] + sout[59] + xorout[58]
   - xorout[59] >= -2
 R1990: sout[3] - sout[19] - sout[27] + sout[35] + sout[59] + xorout[58]
   - xorout[59] >= -2
 R1991: - sout[3] + sout[19] - sout[27] + sout[35] + sout[59] + xorout[58]
   - xorout[59] >= -2
 R1992: - sout[3] - sout[19] + sout[27] + sout[35] + sout[59] + xorout[58]
   - xorout[59] >= -2
 R1993: sout[5] + sout[21] + sout[29] + sout[37] - sout[61] + xorout[61]
   >= 0
 R1994: sout[5] + sout[21] + sout[29] - sout[37] + sout[61] + xorout[61]
   >= 0
 R1995: sout[5] + sout[21] - sout[29] + sout[37] + sout[61] + xorout[61]
   >= 0
 R1996: sout[5] - sout[21] + sout[29] + sout[37] + sout[61] + xorout[61]
   >= 0
 R1997: - sout[5] + sout[21] + sout[29] + sout[37] + sout[61] + xorout[61]
   >= 0
 R1998: sout[5] + sout[21] + sout[29] + sout[37] + sout[61] - xorout[61]
   >= 0
 R1999: sout[5] + sout[21] - sout[29] - sout[37] - sout[61] + xorout[61]
   >= -2
 R2000: sout[5] - sout[21] + sout[29] - sout[37] - sout[61] + xorout[61]
   >= -2
 R2001: - sout[5] + sout[21] + sout[29] - sout[37] - sout[61] + xorout[61]
   >= -2
 R2002: sout[5] - sout[21] - sout[29] + sout[37] - sout[61] + xorout[61]
   >= -2
 R2003: - sout[5] + sout[21] - sout[29] + sout[37] - sout[61] + xorout[61]
   >= -2
 R2004: - sout[5] - sout[21] + sout[29] + sout[37] - sout[61] + xorout[61]
   >= -2
 R2005: sout[5] - sout[21] - sout[29] - sout[37] + sout[61] + xorout[61]
   >= -2
 R2006: - sout[5] + sout[21] - sout[29] - sout[37] + sout[61] + xorout[61]
   >= -2
 R2007: - sout[5] - sout[21] + sout[29] - sout[37] + sout[61] + xorout[61]
   >= -2
 R2008: - sout[5] - sout[21] - sout[29] + sout[37] + sout[61] + xorout[61]
   >= -2
 R2009: - sout[5] - sout[21] - sout[29] - sout[37] - sout[61] + xorout[61]
   >= -4
 R2010: - sout[60] + xorout[60] >= 0
 R2011: - sout[36] + xorout[60] >= 0
 R2012: - sout[28] + xorout[60] >= 0
 R2013: - sout[20] + xorout[60] >= 0
 R2014: - sout[4] + xorout[60] >= 0
 R2015: - sout[4] + sout[5] >= 0
 R2016: - sout[20] + sout[21] >= 0
 R2017: - sout[28] + sout[29] >= 0
 R2018: - sout[36] + sout[37] >= 0
 R2019: - sout[60] + sout[61] >= 0
 R2020: sout[4] + sout[20] + sout[28] + sout[36] + sout[60] - xorout[60]
   >= 0
 R2021: - xorout[60] + xorout[61] >= 0
 R2022: sout[5] - sout[21] - sout[29] - sout[37] - sout[61] + xorout[60]
   - xorout[61] >= -4
 R2023: - sout[5] + sout[21] - sout[29] - sout[37] - sout[61] + xorout[60]
   - xorout[61] >= -4
 R2024: - sout[5] - sout[21] + sout[29] - sout[37] - sout[61] + xorout[60]
   - xorout[61] >= -4
 R2025: - sout[5] - sout[21] - sout[29] + sout[37] - sout[61] + xorout[60]
   - xorout[61] >= -4
 R2026: - sout[5] - sout[21] - sout[29] - sout[37] + sout[61] + xorout[60]
   - xorout[61] >= -4
 R2027: sout[5] + sout[21] + sout[29] - sout[37] - sout[61] + xorout[60]
   - xorout[61] >= -2
 R2028: sout[5] + sout[21] - sout[29] + sout[37] - sout[61] + xorout[60]
   - xorout[61] >= -2
 R2029: sout[5] - sout[21] + sout[29] + sout[37] - sout[61] + xorout[60]
   - xorout[61] >= -2
 R2030: - sout[5] + sout[21] + sout[29] + sout[37] - sout[61] + xorout[60]
   - xorout[61] >= -2
 R2031: sout[5] + sout[21] - sout[29] - sout[37] + sout[61] + xorout[60]
   - xorout[61] >= -2
 R2032: sout[5] - sout[21] + sout[29] - sout[37] + sout[61] + xorout[60]
   - xorout[61] >= -2
 R2033: - sout[5] + sout[21] + sout[29] - sout[37] + sout[61] + xorout[60]
   - xorout[61] >= -2
 R2034: sout[5] - sout[21] - sout[29] + sout[37] + sout[61] + xorout[60]
   - xorout[61] >= -2
 R2035: - sout[5] + sout[21] - sout[29] + sout[37] + sout[61] + xorout[60]
   - xorout[61] >= -2
 R2036: - sout[5] - sout[21] + sout[29] + sout[37] + sout[61] + xorout[60]
   - xorout[61] >= -2
 R2037: sout[7] + sout[23] + sout[31] + sout[39] - sout[63] + xorout[63]
   >= 0
 R2038: sout[7] + sout[23] + sout[31] - sout[39] + sout[63] + xorout[63]
   >= 0
 R2039: sout[7] + sout[23] - sout[31] + sout[39] + sout[63] + xorout[63]
   >= 0
 R2040: sout[7] - sout[23] + sout[31] + sout[39] + sout[63] + xorout[63]
   >= 0
 R2041: - sout[7] + sout[23] + sout[31] + sout[39] + sout[63] + xorout[63]
   >= 0
 R2042: sout[7] + sout[23] + sout[31] + sout[39] + sout[63] - xorout[63]
   >= 0
 R2043: sout[7] + sout[23] - sout[31] - sout[39] - sout[63] + xorout[63]
   >= -2
 R2044: sout[7] - sout[23] + sout[31] - sout[39] - sout[63] + xorout[63]
   >= -2
 R2045: - sout[7] + sout[23] + sout[31] - sout[39] - sout[63] + xorout[63]
   >= -2
 R2046: sout[7] - sout[23] - sout[31] + sout[39] - sout[63] + xorout[63]
   >= -2
 R2047: - sout[7] + sout[23] - sout[31] + sout[39] - sout[63] + xorout[63]
   >= -2
 R2048: - sout[7] - sout[23] + sout[31] + sout[39] - sout[63] + xorout[63]
   >= -2
 R2049: sout[7] - sout[23] - sout[31] - sout[39] + sout[63] + xorout[63]
   >= -2
 R2050: - sout[7] + sout[23] - sout[31] - sout[39] + sout[63] + xorout[63]
   >= -2
 R2051: - sout[7] - sout[23] + sout[31] - sout[39] + sout[63] + xorout[63]
   >= -2
 R2052: - sout[7] - sout[23] - sout[31] + sout[39] + sout[63] + xorout[63]
   >= -2
 R2053: - sout[7] - sout[23] - sout[31] - sout[39] - sout[63] + xorout[63]
   >= -4
 R2054: - sout[62] + xorout[62] >= 0
 R2055: - sout[38] + xorout[62] >= 0
 R2056: - sout[30] + xorout[62] >= 0
 R2057: - sout[22] + xorout[62] >= 0
 R2058: - sout[6] + xorout[62] >= 0
 R2059: - sout[6] + sout[7] >= 0
 R2060: - sout[22] + sout[23] >= 0
 R2061: - sout[30] + sout[31] >= 0
 R2062: - sout[38] + sout[39] >= 0
 R2063: - sout[62] + sout[63] >= 0
 R2064: sout[6] + sout[22] + sout[30] + sout[38] + sout[62] - xorout[62]
   >= 0
 R2065: - xorout[62] + xorout[63] >= 0
 R2066: sout[7] - sout[23] - sout[31] - sout[39] - sout[63] + xorout[62]
   - xorout[63] >= -4
 R2067: - sout[7] + sout[23] - sout[31] - sout[39] - sout[63] + xorout[62]
   - xorout[63] >= -4
 R2068: - sout[7] - sout[23] + sout[31] - sout[39] - sout[63] + xorout[62]
   - xorout[63] >= -4
 R2069: - sout[7] - sout[23] - sout[31] + sout[39] - sout[63] + xorout[62]
   - xorout[63] >= -4
 R2070: - sout[7] - sout[23] - sout[31] - sout[39] + sout[63] + xorout[62]
   - xorout[63] >= -4
 R2071: sout[7] + sout[23] + sout[31] - sout[39] - sout[63] + xorout[62]
   - xorout[63] >= -2
 R2072: sout[7] + sout[23] - sout[31] + sout[39] - sout[63] + xorout[62]
   - xorout[63] >= -2
 R2073: sout[7] - sout[23] + sout[31] + sout[39] - sout[63] + xorout[62]
   - xorout[63] >= -2
 R2074: - sout[7] + sout[23] + sout[31] + sout[39] - sout[63] + xorout[62]
   - xorout[63] >= -2
 R2075: sout[7] + sout[23] - sout[31] - sout[39] + sout[63] + xorout[62]
   - xorout[63] >= -2
 R2076: sout[7] - sout[23] + sout[31] - sout[39] + sout[63] + xorout[62]
   - xorout[63] >= -2
 R2077: - sout[7] + sout[23] + sout[31] - sout[39] + sout[63] + xorout[62]
   - xorout[63] >= -2
 R2078: sout[7] - sout[23] - sout[31] + sout[39] + sout[63] + xorout[62]
   - xorout[63] >= -2
 R2079: - sout[7] + sout[23] - sout[31] + sout[39] + sout[63] + xorout[62]
   - xorout[63] >= -2
 R2080: - sout[7] - sout[23] + sout[31] + sout[39] + sout[63] + xorout[62]
   - xorout[63] >= -2
 R2081: - x[65] + x[129] + xorout[1] >= 0
 R2082: x[65] + x[129] - xorout[1] >= 0
 R2083: - x[64] + x[128] >= 0
 R2084: x[128] - xorout[0] >= 0
 R2085: x[65] - x[129] + xorout[1] >= 0
 R2086: - xorout[0] + xorout[1] >= 0
 R2087: - x[64] + x[65] >= 0
 R2088: x[64] - x[128] + xorout[0] >= 0
 R2089: - x[128] + x[129] >= 0
 R2090: - x[65] + x[128] - x[129] - xorout[1] >= -2
 R2091: - x[67] + x[131] + xorout[3] >= 0
 R2092: x[67] + x[131] - xorout[3] >= 0
 R2093: - x[66] + x[130] >= 0
 R2094: x[130] - xorout[2] >= 0
 R2095: x[67] - x[131] + xorout[3] >= 0
 R2096: - xorout[2] + xorout[3] >= 0
 R2097: - x[66] + x[67] >= 0
 R2098: x[66] - x[130] + xorout[2] >= 0
 R2099: - x[130] + x[131] >= 0
 R2100: - x[67] + x[130] - x[131] - xorout[3] >= -2
 R2101: - x[69] + x[133] + xorout[5] >= 0
 R2102: x[69] + x[133] - xorout[5] >= 0
 R2103: - x[68] + x[132] >= 0
 R2104: x[132] - xorout[4] >= 0
 R2105: x[69] - x[133] + xorout[5] >= 0
 R2106: - xorout[4] + xorout[5] >= 0
 R2107: - x[68] + x[69] >= 0
 R2108: x[68] - x[132] + xorout[4] >= 0
 R2109: - x[132] + x[133] >= 0
 R2110: - x[69] + x[132] - x[133] - xorout[5] >= -2
 R2111: - x[71] + x[135] + xorout[7] >= 0
 R2112: x[71] + x[135] - xorout[7] >= 0
 R2113: - x[70] + x[134] >= 0
 R2114: x[134] - xorout[6] >= 0
 R2115: x[71] - x[135] + xorout[7] >= 0
 R2116: - xorout[6] + xorout[7] >= 0
 R2117: - x[70] + x[71] >= 0
 R2118: x[70] - x[134] + xorout[6] >= 0
 R2119: - x[134] + x[135] >= 0
 R2120: - x[71] + x[134] - x[135] - xorout[7] >= -2
 R2121: - x[73] + x[137] + xorout[9] >= 0
 R2122: x[73] + x[137] - xorout[9] >= 0
 R2123: - x[72] + x[136] >= 0
 R2124: x[136] - xorout[8] >= 0
 R2125: x[73] - x[137] + xorout[9] >= 0
 R2126: - xorout[8] + xorout[9] >= 0
 R2127: - x[72] + x[73] >= 0
 R2128: x[72] - x[136] + xorout[8] >= 0
 R2129: - x[136] + x[137] >= 0
 R2130: - x[73] + x[136] - x[137] - xorout[9] >= -2
 R2131: - x[75] + x[139] + xorout[11] >= 0
 R2132: x[75] + x[139] - xorout[11] >= 0
 R2133: - x[74] + x[138] >= 0
 R2134: x[138] - xorout[10] >= 0
 R2135: x[75] - x[139] + xorout[11] >= 0
 R2136: - xorout[10] + xorout[11] >= 0
 R2137: - x[74] + x[75] >= 0
 R2138: x[74] - x[138] + xorout[10] >= 0
 R2139: - x[138] + x[139] >= 0
 R2140: - x[75] + x[138] - x[139] - xorout[11] >= -2
 R2141: - x[77] + x[141] + xorout[13] >= 0
 R2142: x[77] + x[141] - xorout[13] >= 0
 R2143: - x[76] + x[140] >= 0
 R2144: x[140] - xorout[12] >= 0
 R2145: x[77] - x[141] + xorout[13] >= 0
 R2146: - xorout[12] + xorout[13] >= 0
 R2147: - x[76] + x[77] >= 0
 R2148: x[76] - x[140] + xorout[12] >= 0
 R2149: - x[140] + x[141] >= 0
 R2150: - x[77] + x[140] - x[141] - xorout[13] >= -2
 R2151: - x[79] + x[143] + xorout[15] >= 0
 R2152: x[79] + x[143] - xorout[15] >= 0
 R2153: - x[78] + x[142] >= 0
 R2154: x[142] - xorout[14] >= 0
 R2155: x[79] - x[143] + xorout[15] >= 0
 R2156: - xorout[14] + xorout[15] >= 0
 R2157: - x[78] + x[79] >= 0
 R2158: x[78] - x[142] + xorout[14] >= 0
 R2159: - x[142] + x[143] >= 0
 R2160: - x[79] + x[142] - x[143] - xorout[15] >= -2
 R2161: - x[81] + x[145] + xorout[17] >= 0
 R2162: x[81] + x[145] - xorout[17] >= 0
 R2163: - x[80] + x[144] >= 0
 R2164: x[144] - xorout[16] >= 0
 R2165: x[81] - x[145] + xorout[17] >= 0
 R2166: - xorout[16] + xorout[17] >= 0
 R2167: - x[80] + x[81] >= 0
 R2168: x[80] - x[144] + xorout[16] >= 0
 R2169: - x[144] + x[145] >= 0
 R2170: - x[81] + x[144] - x[145] - xorout[17] >= -2
 R2171: - x[83] + x[147] + xorout[19] >= 0
 R2172: x[83] + x[147] - xorout[19] >= 0
 R2173: - x[82] + x[146] >= 0
 R2174: x[146] - xorout[18] >= 0
 R2175: x[83] - x[147] + xorout[19] >= 0
 R2176: - xorout[18] + xorout[19] >= 0
 R2177: - x[82] + x[83] >= 0
 R2178: x[82] - x[146] + xorout[18] >= 0
 R2179: - x[146] + x[147] >= 0
 R2180: - x[83] + x[146] - x[147] - xorout[19] >= -2
 R2181: - x[85] + x[149] + xorout[21] >= 0
 R2182: x[85] + x[149] - xorout[21] >= 0
 R2183: - x[84] + x[148] >= 0
 R2184: x[148] - xorout[20] >= 0
 R2185: x[85] - x[149] + xorout[21] >= 0
 R2186: - xorout[20] + xorout[21] >= 0
 R2187: - x[84] + x[85] >= 0
 R2188: x[84] - x[148] + xorout[20] >= 0
 R2189: - x[148] + x[149] >= 0
 R2190: - x[85] + x[148] - x[149] - xorout[21] >= -2
 R2191: - x[87] + x[151] + xorout[23] >= 0
 R2192: x[87] + x[151] - xorout[23] >= 0
 R2193: - x[86] + x[150] >= 0
 R2194: x[150] - xorout[22] >= 0
 R2195: x[87] - x[151] + xorout[23] >= 0
 R2196: - xorout[22] + xorout[23] >= 0
 R2197: - x[86] + x[87] >= 0
 R2198: x[86] - x[150] + xorout[22] >= 0
 R2199: - x[150] + x[151] >= 0
 R2200: - x[87] + x[150] - x[151] - xorout[23] >= -2
 R2201: - x[89] + x[153] + xorout[25] >= 0
 R2202: x[89] + x[153] - xorout[25] >= 0
 R2203: - x[88] + x[152] >= 0
 R2204: x[152] - xorout[24] >= 0
 R2205: x[89] - x[153] + xorout[25] >= 0
 R2206: - xorout[24] + xorout[25] >= 0
 R2207: - x[88] + x[89] >= 0
 R2208: x[88] - x[152] + xorout[24] >= 0
 R2209: - x[152] + x[153] >= 0
 R2210: - x[89] + x[152] - x[153] - xorout[25] >= -2
 R2211: - x[91] + x[155] + xorout[27] >= 0
 R2212: x[91] + x[155] - xorout[27] >= 0
 R2213: - x[90] + x[154] >= 0
 R2214: x[154] - xorout[26] >= 0
 R2215: x[91] - x[155] + xorout[27] >= 0
 R2216: - xorout[26] + xorout[27] >= 0
 R2217: - x[90] + x[91] >= 0
 R2218: x[90] - x[154] + xorout[26] >= 0
 R2219: - x[154] + x[155] >= 0
 R2220: - x[91] + x[154] - x[155] - xorout[27] >= -2
 R2221: - x[93] + x[157] + xorout[29] >= 0
 R2222: x[93] + x[157] - xorout[29] >= 0
 R2223: - x[92] + x[156] >= 0
 R2224: x[156] - xorout[28] >= 0
 R2225: x[93] - x[157] + xorout[29] >= 0
 R2226: - xorout[28] + xorout[29] >= 0
 R2227: - x[92] + x[93] >= 0
 R2228: x[92] - x[156] + xorout[28] >= 0
 R2229: - x[156] + x[157] >= 0
 R2230: - x[93] + x[156] - x[157] - xorout[29] >= -2
 R2231: - x[95] + x[159] + xorout[31] >= 0
 R2232: x[95] + x[159] - xorout[31] >= 0
 R2233: - x[94] + x[158] >= 0
 R2234: x[158] - xorout[30] >= 0
 R2235: x[95] - x[159] + xorout[31] >= 0
 R2236: - xorout[30] + xorout[31] >= 0
 R2237: - x[94] + x[95] >= 0
 R2238: x[94] - x[158] + xorout[30] >= 0
 R2239: - x[158] + x[159] >= 0
 R2240: - x[95] + x[158] - x[159] - xorout[31] >= -2
 R2241: - x[97] + x[161] + xorout[33] >= 0
 R2242: x[97] + x[161] - xorout[33] >= 0
 R2243: - x[96] + x[160] >= 0
 R2244: x[160] - xorout[32] >= 0
 R2245: x[97] - x[161] + xorout[33] >= 0
 R2246: - xorout[32] + xorout[33] >= 0
 R2247: - x[96] + x[97] >= 0
 R2248: x[96] - x[160] + xorout[32] >= 0
 R2249: - x[160] + x[161] >= 0
 R2250: - x[97] + x[160] - x[161] - xorout[33] >= -2
 R2251: - x[99] + x[163] + xorout[35] >= 0
 R2252: x[99] + x[163] - xorout[35] >= 0
 R2253: - x[98] + x[162] >= 0
 R2254: x[162] - xorout[34] >= 0
 R2255: x[99] - x[163] + xorout[35] >= 0
 R2256: - xorout[34] + xorout[35] >= 0
 R2257: - x[98] + x[99] >= 0
 R2258: x[98] - x[162] + xorout[34] >= 0
 R2259: - x[162] + x[163] >= 0
 R2260: - x[99] + x[162] - x[163] - xorout[35] >= -2
 R2261: - x[101] + x[165] + xorout[37] >= 0
 R2262: x[101] + x[165] - xorout[37] >= 0
 R2263: - x[100] + x[164] >= 0
 R2264: x[164] - xorout[36] >= 0
 R2265: x[101] - x[165] + xorout[37] >= 0
 R2266: - xorout[36] + xorout[37] >= 0
 R2267: - x[100] + x[101] >= 0
 R2268: x[100] - x[164] + xorout[36] >= 0
 R2269: - x[164] + x[165] >= 0
 R2270: - x[101] + x[164] - x[165] - xorout[37] >= -2
 R2271: - x[103] + x[167] + xorout[39] >= 0
 R2272: x[103] + x[167] - xorout[39] >= 0
 R2273: - x[102] + x[166] >= 0
 R2274: x[166] - xorout[38] >= 0
 R2275: x[103] - x[167] + xorout[39] >= 0
 R2276: - xorout[38] + xorout[39] >= 0
 R2277: - x[102] + x[103] >= 0
 R2278: x[102] - x[166] + xorout[38] >= 0
 R2279: - x[166] + x[167] >= 0
 R2280: - x[103] + x[166] - x[167] - xorout[39] >= -2
 R2281: - x[105] + x[169] + xorout[41] >= 0
 R2282: x[105] + x[169] - xorout[41] >= 0
 R2283: - x[104] + x[168] >= 0
 R2284: x[168] - xorout[40] >= 0
 R2285: x[105] - x[169] + xorout[41] >= 0
 R2286: - xorout[40] + xorout[41] >= 0
 R2287: - x[104] + x[105] >= 0
 R2288: x[104] - x[168] + xorout[40] >= 0
 R2289: - x[168] + x[169] >= 0
 R2290: - x[105] + x[168] - x[169] - xorout[41] >= -2
 R2291: - x[107] + x[171] + xorout[43] >= 0
 R2292: x[107] + x[171] - xorout[43] >= 0
 R2293: - x[106] + x[170] >= 0
 R2294: x[170] - xorout[42] >= 0
 R2295: x[107] - x[171] + xorout[43] >= 0
 R2296: - xorout[42] + xorout[43] >= 0
 R2297: - x[106] + x[107] >= 0
 R2298: x[106] - x[170] + xorout[42] >= 0
 R2299: - x[170] + x[171] >= 0
 R2300: - x[107] + x[170] - x[171] - xorout[43] >= -2
 R2301: - x[109] + x[173] + xorout[45] >= 0
 R2302: x[109] + x[173] - xorout[45] >= 0
 R2303: - x[108] + x[172] >= 0
 R2304: x[172] - xorout[44] >= 0
 R2305: x[109] - x[173] + xorout[45] >= 0
 R2306: - xorout[44] + xorout[45] >= 0
 R2307: - x[108] + x[109] >= 0
 R2308: x[108] - x[172] + xorout[44] >= 0
 R2309: - x[172] + x[173] >= 0
 R2310: - x[109] + x[172] - x[173] - xorout[45] >= -2
 R2311: - x[111] + x[175] + xorout[47] >= 0
 R2312: x[111] + x[175] - xorout[47] >= 0
 R2313: - x[110] + x[174] >= 0
 R2314: x[174] - xorout[46] >= 0
 R2315: x[111] - x[175] + xorout[47] >= 0
 R2316: - xorout[46] + xorout[47] >= 0
 R2317: - x[110] + x[111] >= 0
 R2318: x[110] - x[174] + xorout[46] >= 0
 R2319: - x[174] + x[175] >= 0
 R2320: - x[111] + x[174] - x[175] - xorout[47] >= -2
 R2321: - x[113] + x[177] + xorout[49] >= 0
 R2322: x[113] + x[177] - xorout[49] >= 0
 R2323: - x[112] + x[176] >= 0
 R2324: x[176] - xorout[48] >= 0
 R2325: x[113] - x[177] + xorout[49] >= 0
 R2326: - xorout[48] + xorout[49] >= 0
 R2327: - x[112] + x[113] >= 0
 R2328: x[112] - x[176] + xorout[48] >= 0
 R2329: - x[176] + x[177] >= 0
 R2330: - x[113] + x[176] - x[177] - xorout[49] >= -2
 R2331: - x[115] + x[179] + xorout[51] >= 0
 R2332: x[115] + x[179] - xorout[51] >= 0
 R2333: - x[114] + x[178] >= 0
 R2334: x[178] - xorout[50] >= 0
 R2335: x[115] - x[179] + xorout[51] >= 0
 R2336: - xorout[50] + xorout[51] >= 0
 R2337: - x[114] + x[115] >= 0
 R2338: x[114] - x[178] + xorout[50] >= 0
 R2339: - x[178] + x[179] >= 0
 R2340: - x[115] + x[178] - x[179] - xorout[51] >= -2
 R2341: - x[117] + x[181] + xorout[53] >= 0
 R2342: x[117] + x[181] - xorout[53] >= 0
 R2343: - x[116] + x[180] >= 0
 R2344: x[180] - xorout[52] >= 0
 R2345: x[117] - x[181] + xorout[53] >= 0
 R2346: - xorout[52] + xorout[53] >= 0
 R2347: - x[116] + x[117] >= 0
 R2348: x[116] - x[180] + xorout[52] >= 0
 R2349: - x[180] + x[181] >= 0
 R2350: - x[117] + x[180] - x[181] - xorout[53] >= -2
 R2351: - x[119] + x[183] + xorout[55] >= 0
 R2352: x[119] + x[183] - xorout[55] >= 0
 R2353: - x[118] + x[182] >= 0
 R2354: x[182] - xorout[54] >= 0
 R2355: x[119] - x[183] + xorout[55] >= 0
 R2356: - xorout[54] + xorout[55] >= 0
 R2357: - x[118] + x[119] >= 0
 R2358: x[118] - x[182] + xorout[54] >= 0
 R2359: - x[182] + x[183] >= 0
 R2360: - x[119] + x[182] - x[183] - xorout[55] >= -2
 R2361: - x[121] + x[185] + xorout[57] >= 0
 R2362: x[121] + x[185] - xorout[57] >= 0
 R2363: - x[120] + x[184] >= 0
 R2364: x[184] - xorout[56] >= 0
 R2365: x[121] - x[185] + xorout[57] >= 0
 R2366: - xorout[56] + xorout[57] >= 0
 R2367: - x[120] + x[121] >= 0
 R2368: x[120] - x[184] + xorout[56] >= 0
 R2369: - x[184] + x[185] >= 0
 R2370: - x[121] + x[184] - x[185] - xorout[57] >= -2
 R2371: - x[123] + x[187] + xorout[59] >= 0
 R2372: x[123] + x[187] - xorout[59] >= 0
 R2373: - x[122] + x[186] >= 0
 R2374: x[186] - xorout[58] >= 0
 R2375: x[123] - x[187] + xorout[59] >= 0
 R2376: - xorout[58] + xorout[59] >= 0
 R2377: - x[122] + x[123] >= 0
 R2378: x[122] - x[186] + xorout[58] >= 0
 R2379: - x[186] + x[187] >= 0
 R2380: - x[123] + x[186] - x[187] - xorout[59] >= -2
 R2381: - x[125] + x[189] + xorout[61] >= 0
 R2382: x[125] + x[189] - xorout[61] >= 0
 R2383: - x[124] + x[188] >= 0
 R2384: x[188] - xorout[60] >= 0
 R2385: x[125] - x[189] + xorout[61] >= 0
 R2386: - xorout[60] + xorout[61] >= 0
 R2387: - x[124] + x[125] >= 0
 R2388: x[124] - x[188] + xorout[60] >= 0
 R2389: - x[188] + x[189] >= 0
 R2390: - x[125] + x[188] - x[189] - xorout[61] >= -2
 R2391: - x[127] + x[191] + xorout[63] >= 0
 R2392: x[127] + x[191] - xorout[63] >= 0
 R2393: - x[126] + x[190] >= 0
 R2394: x[190] - xorout[62] >= 0
 R2395: x[127] - x[191] + xorout[63] >= 0
 R2396: - xorout[62] + xorout[63] >= 0
 R2397: - x[126] + x[127] >= 0
 R2398: x[126] - x[190] + xorout[62] >= 0
 R2399: - x[190] + x[191] >= 0
 R2400: - x[127] + x[190] - x[191] - xorout[63] >= -2
 R2401: x[0] - x[192] = 0
 R2402: x[1] - x[193] = 0
 R2403: x[2] - x[194] = 0
 R2404: x[3] - x[195] = 0
 R2405: x[4] - x[196] = 0
 R2406: x[5] - x[197] = 0
 R2407: x[6] - x[198] = 0
 R2408: x[7] - x[199] = 0
 R2409: x[8] - x[200] = 0
 R2410: x[9] - x[201] = 0
 R2411: x[10] - x[202] = 0
 R2412: x[11] - x[203] = 0
 R2413: x[12] - x[204] = 0
 R2414: x[13] - x[205] = 0
 R2415: x[14] - x[206] = 0
 R2416: x[15] - x[207] = 0
 R2417: x[16] - x[208] = 0
 R2418: x[17] - x[209] = 0
 R2419: x[18] - x[210] = 0
 R2420: x[19] - x[211] = 0
 R2421: x[20] - x[212] = 0
 R2422: x[21] - x[213] = 0
 R2423: x[22] - x[214] = 0
 R2424: x[23] - x[215] = 0
 R2425: x[24] - x[216] = 0
 R2426: x[25] - x[217] = 0
 R2427: x[26] - x[218] = 0
 R2428: x[27] - x[219] = 0
 R2429: x[28] - x[220] = 0
 R2430: x[29] - x[221] = 0
 R2431: x[30] - x[222] = 0
 R2432: x[31] - x[223] = 0
 R2433: x[32] - x[224] = 0
 R2434: x[33] - x[225] = 0
 R2435: x[34] - x[226] = 0
 R2436: x[35] - x[227] = 0
 R2437: x[36] - x[228] = 0
 R2438: x[37] - x[229] = 0
 R2439: x[38] - x[230] = 0
 R2440: x[39] - x[231] = 0
 R2441: x[40] - x[232] = 0
 R2442: x[41] - x[233] = 0
 R2443: x[42] - x[234] = 0
 R2444: x[43] - x[235] = 0
 R2445: x[44] - x[236] = 0
 R2446: x[45] - x[237] = 0
 R2447: x[46] - x[238] = 0
 R2448: x[47] - x[239] = 0
 R2449: x[48] - x[240] = 0
 R2450: x[49] - x[241] = 0
 R2451: x[50] - x[242] = 0
 R2452: x[51] - x[243] = 0
 R2453: x[52] - x[244] = 0
 R2454: x[53] - x[245] = 0
 R2455: x[54] - x[246] = 0
 R2456: x[55] - x[247] = 0
 R2457: x[56] - x[248] = 0
 R2458: x[57] - x[249] = 0
 R2459: x[58] - x[250] = 0
 R2460: x[59] - x[251] = 0
 R2461: x[60] - x[252] = 0
 R2462: x[61] - x[253] = 0
 R2463: x[62] - x[254] = 0
 R2464: x[63] - x[255] = 0
 R2465: - x[128] + x[129] >= 0
 R2466: - x[130] + x[131] >= 0
 R2467: - x[132] + x[133] >= 0
 R2468: - x[134] + x[135] >= 0
 R2469: - x[135] + sout[71] >= 0
 R2470: sout[64] - sout[65] >= 0
 R2471: sout[65] - sout[66] >= 0
 R2472: sout[66] - sout[67] >= 0
 R2473: sout[67] - sout[68] >= 0
 R2474: sout[68] - sout[69] >= 0
 R2475: sout[69] - sout[70] >= 0
 R2476: sout[70] - sout[71] >= 0
 R2477: - x[133] + sout[71] >= 0
 R2478: - x[131] + sout[71] >= 0
 R2479: - x[129] + sout[71] >= 0
 R2480: x[129] + x[131] + x[133] + x[135] - sout[64] >= 0
 R2481: - x[136] + x[137] >= 0
 R2482: - x[138] + x[139] >= 0
 R2483: - x[140] + x[141] >= 0
 R2484: - x[142] + x[143] >= 0
 R2485: - x[143] + sout[79] >= 0
 R2486: sout[72] - sout[73] >= 0
 R2487: sout[73] - sout[74] >= 0
 R2488: sout[74] - sout[75] >= 0
 R2489: sout[75] - sout[76] >= 0
 R2490: sout[76] - sout[77] >= 0
 R2491: sout[77] - sout[78] >= 0
 R2492: sout[78] - sout[79] >= 0
 R2493: - x[141] + sout[79] >= 0
 R2494: - x[139] + sout[79] >= 0
 R2495: - x[137] + sout[79] >= 0
 R2496: x[137] + x[139] + x[141] + x[143] - sout[72] >= 0
 R2497: - x[144] + x[145] >= 0
 R2498: - x[146] + x[147] >= 0
 R2499: - x[148] + x[149] >= 0
 R2500: - x[150] + x[151] >= 0
 R2501: - x[151] + sout[87] >= 0
 R2502: sout[80] - sout[81] >= 0
 R2503: sout[81] - sout[82] >= 0
 R2504: sout[82] - sout[83] >= 0
 R2505: sout[83] - sout[84] >= 0
 R2506: sout[84] - sout[85] >= 0
 R2507: sout[85] - sout[86] >= 0
 R2508: sout[86] - sout[87] >= 0
 R2509: - x[149] + sout[87] >= 0
 R2510: - x[147] + sout[87] >= 0
 R2511: - x[145] + sout[87] >= 0
 R2512: x[145] + x[147] + x[149] + x[151] - sout[80] >= 0
 R2513: - x[152] + x[153] >= 0
 R2514: - x[154] + x[155] >= 0
 R2515: - x[156] + x[157] >= 0
 R2516: - x[158] + x[159] >= 0
 R2517: - x[159] + sout[95] >= 0
 R2518: sout[88] - sout[89] >= 0
 R2519: sout[89] - sout[90] >= 0
 R2520: sout[90] - sout[91] >= 0
 R2521: sout[91] - sout[92] >= 0
 R2522: sout[92] - sout[93] >= 0
 R2523: sout[93] - sout[94] >= 0
 R2524: sout[94] - sout[95] >= 0
 R2525: - x[157] + sout[95] >= 0
 R2526: - x[155] + sout[95] >= 0
 R2527: - x[153] + sout[95] >= 0
 R2528: x[153] + x[155] + x[157] + x[159] - sout[88] >= 0
 R2529: - x[160] + x[161] >= 0
 R2530: - x[162] + x[163] >= 0
 R2531: - x[164] + x[165] >= 0
 R2532: - x[166] + x[167] >= 0
 R2533: - x[167] + sout[103] >= 0
 R2534: sout[96] - sout[97] >= 0
 R2535: sout[97] - sout[98] >= 0
 R2536: sout[98] - sout[99] >= 0
 R2537: sout[99] - sout[100] >= 0
 R2538: sout[100] - sout[101] >= 0
 R2539: sout[101] - sout[102] >= 0
 R2540: sout[102] - sout[103] >= 0
 R2541: - x[165] + sout[103] >= 0
 R2542: - x[163] + sout[103] >= 0
 R2543: - x[161] + sout[103] >= 0
 R2544: x[161] + x[163] + x[165] + x[167] - sout[96] >= 0
 R2545: - x[168] + x[169] >= 0
 R2546: - x[170] + x[171] >= 0
 R2547: - x[172] + x[173] >= 0
 R2548: - x[174] + x[175] >= 0
 R2549: - x[175] + sout[111] >= 0
 R2550: sout[104] - sout[105] >= 0
 R2551: sout[105] - sout[106] >= 0
 R2552: sout[106] - sout[107] >= 0
 R2553: sout[107] - sout[108] >= 0
 R2554: sout[108] - sout[109] >= 0
 R2555: sout[109] - sout[110] >= 0
 R2556: sout[110] - sout[111] >= 0
 R2557: - x[173] + sout[111] >= 0
 R2558: - x[171] + sout[111] >= 0
 R2559: - x[169] + sout[111] >= 0
 R2560: x[169] + x[171] + x[173] + x[175] - sout[104] >= 0
 R2561: - x[176] + x[177] >= 0
 R2562: - x[178] + x[179] >= 0
 R2563: - x[180] + x[181] >= 0
 R2564: - x[182] + x[183] >= 0
 R2565: - x[183] + sout[119] >= 0
 R2566: sout[112] - sout[113] >= 0
 R2567: sout[113] - sout[114] >= 0
 R2568: sout[114] - sout[115] >= 0
 R2569: sout[115] - sout[116] >= 0
 R2570: sout[116] - sout[117] >= 0
 R2571: sout[117] - sout[118] >= 0
 R2572: sout[118] - sout[119] >= 0
 R2573: - x[181] + sout[119] >= 0
 R2574: - x[179] + sout[119] >= 0
 R2575: - x[177] + sout[119] >= 0
 R2576: x[177] + x[179] + x[181] + x[183] - sout[112] >= 0
 R2577: - x[184] + x[185] >= 0
 R2578: - x[186] + x[187] >= 0
 R2579: - x[188] + x[189] >= 0
 R2580: - x[190] + x[191] >= 0
 R2581: - x[191] + sout[127] >= 0
 R2582: sout[120] - sout[121] >= 0
 R2583: sout[121] - sout[122] >= 0
 R2584: sout[122] - sout[123] >= 0
 R2585: sout[123] - sout[124] >= 0
 R2586: sout[124] - sout[125] >= 0
 R2587: sout[125] - sout[126] >= 0
 R2588: sout[126] - sout[127] >= 0
 R2589: - x[189] + sout[127] >= 0
 R2590: - x[187] + sout[127] >= 0
 R2591: - x[185] + sout[127] >= 0
 R2592: x[185] + x[187] + x[189] + x[191] - sout[120] >= 0
 R2593: sout[73] + sout[81] + sout[89] + sout[97] + sout[105] - sout[113]
   + xorout[65] >= 0
 R2594: sout[73] + sout[81] + sout[89] + sout[97] - sout[105] + sout[113]
   + xorout[65] >= 0
 R2595: sout[73] + sout[81] + sout[89] - sout[97] + sout[105] + sout[113]
   + xorout[65] >= 0
 R2596: sout[73] + sout[81] - sout[89] + sout[97] + sout[105] + sout[113]
   + xorout[65] >= 0
 R2597: sout[73] - sout[81] + sout[89] + sout[97] + sout[105] + sout[113]
   + xorout[65] >= 0
 R2598: - sout[73] + sout[81] + sout[89] + sout[97] + sout[105] + sout[113]
   + xorout[65] >= 0
 R2599: sout[73] + sout[81] + sout[89] + sout[97] + sout[105] + sout[113]
   - xorout[65] >= 0
 R2600: sout[73] + sout[81] + sout[89] - sout[97] - sout[105] - sout[113]
   + xorout[65] >= -2
 R2601: sout[73] + sout[81] - sout[89] + sout[97] - sout[105] - sout[113]
   + xorout[65] >= -2
 R2602: sout[73] - sout[81] + sout[89] + sout[97] - sout[105] - sout[113]
   + xorout[65] >= -2
 R2603: - sout[73] + sout[81] + sout[89] + sout[97] - sout[105] - sout[113]
   + xorout[65] >= -2
 R2604: sout[73] + sout[81] - sout[89] - sout[97] + sout[105] - sout[113]
   + xorout[65] >= -2
 R2605: sout[73] - sout[81] + sout[89] - sout[97] + sout[105] - sout[113]
   + xorout[65] >= -2
 R2606: - sout[73] + sout[81] + sout[89] - sout[97] + sout[105] - sout[113]
   + xorout[65] >= -2
 R2607: sout[73] - sout[81] - sout[89] + sout[97] + sout[105] - sout[113]
   + xorout[65] >= -2
 R2608: - sout[73] + sout[81] - sout[89] + sout[97] + sout[105] - sout[113]
   + xorout[65] >= -2
 R2609: - sout[73] - sout[81] + sout[89] + sout[97] + sout[105] - sout[113]
   + xorout[65] >= -2
 R2610: sout[73] + sout[81] - sout[89] - sout[97] - sout[105] + sout[113]
   + xorout[65] >= -2
 R2611: sout[73] - sout[81] + sout[89] - sout[97] - sout[105] + sout[113]
   + xorout[65] >= -2
 R2612: - sout[73] + sout[81] + sout[89] - sout[97] - sout[105] + sout[113]
   + xorout[65] >= -2
 R2613: sout[73] - sout[81] - sout[89] + sout[97] - sout[105] + sout[113]
   + xorout[65] >= -2
 R2614: - sout[73] + sout[81] - sout[89] + sout[97] - sout[105] + sout[113]
   + xorout[65] >= -2
 R2615: - sout[73] - sout[81] + sout[89] + sout[97] - sout[105] + sout[113]
   + xorout[65] >= -2
 R2616: sout[73] - sout[81] - sout[89] - sout[97] + sout[105] + sout[113]
   + xorout[65] >= -2
 R2617: - sout[73] + sout[81] - sout[89] - sout[97] + sout[105] + sout[113]
   + xorout[65] >= -2
 R2618: - sout[73] - sout[81] + sout[89] - sout[97] + sout[105] + sout[113]
   + xorout[65] >= -2
 R2619: - sout[73] - sout[81] - sout[89] + sout[97] + sout[105] + sout[113]
   + xorout[65] >= -2
 R2620: sout[73] - sout[81] - sout[89] - sout[97] - sout[105] - sout[113]
   + xorout[65] >= -4
 R2621: - sout[73] + sout[81] - sout[89] - sout[97] - sout[105] - sout[113]
   + xorout[65] >= -4
 R2622: - sout[73] - sout[81] + sout[89] - sout[97] - sout[105] - sout[113]
   + xorout[65] >= -4
 R2623: - sout[73] - sout[81] - sout[89] + sout[97] - sout[105] - sout[113]
   + xorout[65] >= -4
 R2624: - sout[73] - sout[81] - sout[89] - sout[97] + sout[105] - sout[113]
   + xorout[65] >= -4
 R2625: - sout[73] - sout[81] - sout[89] - sout[97] - sout[105] + sout[113]
   + xorout[65] >= -4
 R2626: - sout[112] + xorout[64] >= 0
 R2627: - sout[104] + xorout[64] >= 0
 R2628: - sout[96] + xorout[64] >= 0
 R2629: - sout[88] + xorout[64] >= 0
 R2630: - sout[80] + xorout[64] >= 0
 R2631: - sout[72] + xorout[64] >= 0
 R2632: - sout[72] + sout[73] >= 0
 R2633: - sout[80] + sout[81] >= 0
 R2634: - sout[88] + sout[89] >= 0
 R2635: - sout[96] + sout[97] >= 0
 R2636: - sout[104] + sout[105] >= 0
 R2637: - sout[112] + sout[113] >= 0
 R2638: sout[72] + sout[80] + sout[88] + sout[96] + sout[104] + sout[112]
   - xorout[64] >= 0
 R2639: - xorout[64] + xorout[65] >= 0
 R2640: - sout[73] - sout[81] - sout[89] - sout[97] - sout[105] - sout[113]
   + xorout[64] - xorout[65] >= -6
 R2641: sout[73] + sout[81] - sout[89] - sout[97] - sout[105] - sout[113]
   + xorout[64] - xorout[65] >= -4
 R2642: sout[73] - sout[81] + sout[89] - sout[97] - sout[105] - sout[113]
   + xorout[64] - xorout[65] >= -4
 R2643: - sout[73] + sout[81] + sout[89] - sout[97] - sout[105] - sout[113]
   + xorout[64] - xorout[65] >= -4
 R2644: sout[73] - sout[81] - sout[89] + sout[97] - sout[105] - sout[113]
   + xorout[64] - xorout[65] >= -4
 R2645: - sout[73] + sout[81] - sout[89] + sout[97] - sout[105] - sout[113]
   + xorout[64] - xorout[65] >= -4
 R2646: - sout[73] - sout[81] + sout[89] + sout[97] - sout[105] - sout[113]
   + xorout[64] - xorout[65] >= -4
 R2647: sout[73] - sout[81] - sout[89] - sout[97] + sout[105] - sout[113]
   + xorout[64] - xorout[65] >= -4
 R2648: - sout[73] + sout[81] - sout[89] - sout[97] + sout[105] - sout[113]
   + xorout[64] - xorout[65] >= -4
 R2649: - sout[73] - sout[81] + sout[89] - sout[97] + sout[105] - sout[113]
   + xorout[64] - xorout[65] >= -4
 R2650: - sout[73] - sout[81] - sout[89] + sout[97] + sout[105] - sout[113]
   + xorout[64] - xorout[65] >= -4
 R2651: sout[73] - sout[81] - sout[89] - sout[97] - sout[105] + sout[113]
   + xorout[64] - xorout[65] >= -4
 R2652: - sout[73] + sout[81] - sout[89] - sout[97] - sout[105] + sout[113]
   + xorout[64] - xorout[65] >= -4
 R2653: - sout[73] - sout[81] + sout[89] - sout[97] - sout[105] + sout[113]
   + xorout[64] - xorout[65] >= -4
 R2654: - sout[73] - sout[81] - sout[89] + sout[97] - sout[105] + sout[113]
   + xorout[64] - xorout[65] >= -4
 R2655: - sout[73] - sout[81] - sout[89] - sout[97] + sout[105] + sout[113]
   + xorout[64] - xorout[65] >= -4
 R2656: sout[73] + sout[81] + sout[89] + sout[97] - sout[105] - sout[113]
   + xorout[64] - xorout[65] >= -2
 R2657: sout[73] + sout[81] + sout[89] - sout[97] + sout[105] - sout[113]
   + xorout[64] - xorout[65] >= -2
 R2658: sout[73] + sout[81] - sout[89] + sout[97] + sout[105] - sout[113]
   + xorout[64] - xorout[65] >= -2
 R2659: sout[73] - sout[81] + sout[89] + sout[97] + sout[105] - sout[113]
   + xorout[64] - xorout[65] >= -2
 R2660: - sout[73] + sout[81] + sout[89] + sout[97] + sout[105] - sout[113]
   + xorout[64] - xorout[65] >= -2
 R2661: sout[73] + sout[81] + sout[89] - sout[97] - sout[105] + sout[113]
   + xorout[64] - xorout[65] >= -2
 R2662: sout[73] + sout[81] - sout[89] + sout[97] - sout[105] + sout[113]
   + xorout[64] - xorout[65] >= -2
 R2663: sout[73] - sout[81] + sout[89] + sout[97] - sout[105] + sout[113]
   + xorout[64] - xorout[65] >= -2
 R2664: - sout[73] + sout[81] + sout[89] + sout[97] - sout[105] + sout[113]
   + xorout[64] - xorout[65] >= -2
 R2665: sout[73] + sout[81] - sout[89] - sout[97] + sout[105] + sout[113]
   + xorout[64] - xorout[65] >= -2
 R2666: sout[73] - sout[81] + sout[89] - sout[97] + sout[105] + sout[113]
   + xorout[64] - xorout[65] >= -2
 R2667: - sout[73] + sout[81] + sout[89] - sout[97] + sout[105] + sout[113]
   + xorout[64] - xorout[65] >= -2
 R2668: sout[73] - sout[81] - sout[89] + sout[97] + sout[105] + sout[113]
   + xorout[64] - xorout[65] >= -2
 R2669: - sout[73] + sout[81] - sout[89] + sout[97] + sout[105] + sout[113]
   + xorout[64] - xorout[65] >= -2
 R2670: - sout[73] - sout[81] + sout[89] + sout[97] + sout[105] + sout[113]
   + xorout[64] - xorout[65] >= -2
 R2671: sout[75] + sout[83] + sout[91] + sout[99] + sout[107] - sout[115]
   + xorout[67] >= 0
 R2672: sout[75] + sout[83] + sout[91] + sout[99] - sout[107] + sout[115]
   + xorout[67] >= 0
 R2673: sout[75] + sout[83] + sout[91] - sout[99] + sout[107] + sout[115]
   + xorout[67] >= 0
 R2674: sout[75] + sout[83] - sout[91] + sout[99] + sout[107] + sout[115]
   + xorout[67] >= 0
 R2675: sout[75] - sout[83] + sout[91] + sout[99] + sout[107] + sout[115]
   + xorout[67] >= 0
 R2676: - sout[75] + sout[83] + sout[91] + sout[99] + sout[107] + sout[115]
   + xorout[67] >= 0
 R2677: sout[75] + sout[83] + sout[91] + sout[99] + sout[107] + sout[115]
   - xorout[67] >= 0
 R2678: sout[75] + sout[83] + sout[91] - sout[99] - sout[107] - sout[115]
   + xorout[67] >= -2
 R2679: sout[75] + sout[83] - sout[91] + sout[99] - sout[107] - sout[115]
   + xorout[67] >= -2
 R2680: sout[75] - sout[83] + sout[91] + sout[99] - sout[107] - sout[115]
   + xorout[67] >= -2
 R2681: - sout[75] + sout[83] + sout[91] + sout[99] - sout[107] - sout[115]
   + xorout[67] >= -2
 R2682: sout[75] + sout[83] - sout[91] - sout[99] + sout[107] - sout[115]
   + xorout[67] >= -2
 R2683: sout[75] - sout[83] + sout[91] - sout[99] + sout[107] - sout[115]
   + xorout[67] >= -2
 R2684: - sout[75] + sout[83] + sout[91] - sout[99] + sout[107] - sout[115]
   + xorout[67] >= -2
 R2685: sout[75] - sout[83] - sout[91] + sout[99] + sout[107] - sout[115]
   + xorout[67] >= -2
 R2686: - sout[75] + sout[83] - sout[91] + sout[99] + sout[107] - sout[115]
   + xorout[67] >= -2
 R2687: - sout[75] - sout[83] + sout[91] + sout[99] + sout[107] - sout[115]
   + xorout[67] >= -2
 R2688: sout[75] + sout[83] - sout[91] - sout[99] - sout[107] + sout[115]
   + xorout[67] >= -2
 R2689: sout[75] - sout[83] + sout[91] - sout[99] - sout[107] + sout[115]
   + xorout[67] >= -2
 R2690: - sout[75] + sout[83] + sout[91] - sout[99] - sout[107] + sout[115]
   + xorout[67] >= -2
 R2691: sout[75] - sout[83] - sout[91] + sout[99] - sout[107] + sout[115]
   + xorout[67] >= -2
 R2692: - sout[75] + sout[83] - sout[91] + sout[99] - sout[107] + sout[115]
   + xorout[67] >= -2
 R2693: - sout[75] - sout[83] + sout[91] + sout[99] - sout[107] + sout[115]
   + xorout[67] >= -2
 R2694: sout[75] - sout[83] - sout[91] - sout[99] + sout[107] + sout[115]
   + xorout[67] >= -2
 R2695: - sout[75] + sout[83] - sout[91] - sout[99] + sout[107] + sout[115]
   + xorout[67] >= -2
 R2696: - sout[75] - sout[83] + sout[91] - sout[99] + sout[107] + sout[115]
   + xorout[67] >= -2
 R2697: - sout[75] - sout[83] - sout[91] + sout[99] + sout[107] + sout[115]
   + xorout[67] >= -2
 R2698: sout[75] - sout[83] - sout[91] - sout[99] - sout[107] - sout[115]
   + xorout[67] >= -4
 R2699: - sout[75] + sout[83] - sout[91] - sout[99] - sout[107] - sout[115]
   + xorout[67] >= -4
 R2700: - sout[75] - sout[83] + sout[91] - sout[99] - sout[107] - sout[115]
   + xorout[67] >= -4
 R2701: - sout[75] - sout[83] - sout[91] + sout[99] - sout[107] - sout[115]
   + xorout[67] >= -4
 R2702: - sout[75] - sout[83] - sout[91] - sout[99] + sout[107] - sout[115]
   + xorout[67] >= -4
 R2703: - sout[75] - sout[83] - sout[91] - sout[99] - sout[107] + sout[115]
   + xorout[67] >= -4
 R2704: - sout[114] + xorout[66] >= 0
 R2705: - sout[106] + xorout[66] >= 0
 R2706: - sout[98] + xorout[66] >= 0
 R2707: - sout[90] + xorout[66] >= 0
 R2708: - sout[82] + xorout[66] >= 0
 R2709: - sout[74] + xorout[66] >= 0
 R2710: - sout[74] + sout[75] >= 0
 R2711: - sout[82] + sout[83] >= 0
 R2712: - sout[90] + sout[91] >= 0
 R2713: - sout[98] + sout[99] >= 0
 R2714: - sout[106] + sout[107] >= 0
 R2715: - sout[114] + sout[115] >= 0
 R2716: sout[74] + sout[82] + sout[90] + sout[98] + sout[106] + sout[114]
   - xorout[66] >= 0
 R2717: - xorout[66] + xorout[67] >= 0
 R2718: - sout[75] - sout[83] - sout[91] - sout[99] - sout[107] - sout[115]
   + xorout[66] - xorout[67] >= -6
 R2719: sout[75] + sout[83] - sout[91] - sout[99] - sout[107] - sout[115]
   + xorout[66] - xorout[67] >= -4
 R2720: sout[75] - sout[83] + sout[91] - sout[99] - sout[107] - sout[115]
   + xorout[66] - xorout[67] >= -4
 R2721: - sout[75] + sout[83] + sout[91] - sout[99] - sout[107] - sout[115]
   + xorout[66] - xorout[67] >= -4
 R2722: sout[75] - sout[83] - sout[91] + sout[99] - sout[107] - sout[115]
   + xorout[66] - xorout[67] >= -4
 R2723: - sout[75] + sout[83] - sout[91] + sout[99] - sout[107] - sout[115]
   + xorout[66] - xorout[67] >= -4
 R2724: - sout[75] - sout[83] + sout[91] + sout[99] - sout[107] - sout[115]
   + xorout[66] - xorout[67] >= -4
 R2725: sout[75] - sout[83] - sout[91] - sout[99] + sout[107] - sout[115]
   + xorout[66] - xorout[67] >= -4
 R2726: - sout[75] + sout[83] - sout[91] - sout[99] + sout[107] - sout[115]
   + xorout[66] - xorout[67] >= -4
 R2727: - sout[75] - sout[83] + sout[91] - sout[99] + sout[107] - sout[115]
   + xorout[66] - xorout[67] >= -4
 R2728: - sout[75] - sout[83] - sout[91] + sout[99] + sout[107] - sout[115]
   + xorout[66] - xorout[67] >= -4
 R2729: sout[75] - sout[83] - sout[91] - sout[99] - sout[107] + sout[115]
   + xorout[66] - xorout[67] >= -4
 R2730: - sout[75] + sout[83] - sout[91] - sout[99] - sout[107] + sout[115]
   + xorout[66] - xorout[67] >= -4
 R2731: - sout[75] - sout[83] + sout[91] - sout[99] - sout[107] + sout[115]
   + xorout[66] - xorout[67] >= -4
 R2732: - sout[75] - sout[83] - sout[91] + sout[99] - sout[107] + sout[115]
   + xorout[66] - xorout[67] >= -4
 R2733: - sout[75] - sout[83] - sout[91] - sout[99] + sout[107] + sout[115]
   + xorout[66] - xorout[67] >= -4
 R2734: sout[75] + sout[83] + sout[91] + sout[99] - sout[107] - sout[115]
   + xorout[66] - xorout[67] >= -2
 R2735: sout[75] + sout[83] + sout[91] - sout[99] + sout[107] - sout[115]
   + xorout[66] - xorout[67] >= -2
 R2736: sout[75] + sout[83] - sout[91] + sout[99] + sout[107] - sout[115]
   + xorout[66] - xorout[67] >= -2
 R2737: sout[75] - sout[83] + sout[91] + sout[99] + sout[107] - sout[115]
   + xorout[66] - xorout[67] >= -2
 R2738: - sout[75] + sout[83] + sout[91] + sout[99] + sout[107] - sout[115]
   + xorout[66] - xorout[67] >= -2
 R2739: sout[75] + sout[83] + sout[91] - sout[99] - sout[107] + sout[115]
   + xorout[66] - xorout[67] >= -2
 R2740: sout[75] + sout[83] - sout[91] + sout[99] - sout[107] + sout[115]
   + xorout[66] - xorout[67] >= -2
 R2741: sout[75] - sout[83] + sout[91] + sout[99] - sout[107] + sout[115]
   + xorout[66] - xorout[67] >= -2
 R2742: - sout[75] + sout[83] + sout[91] + sout[99] - sout[107] + sout[115]
   + xorout[66] - xorout[67] >= -2
 R2743: sout[75] + sout[83] - sout[91] - sout[99] + sout[107] + sout[115]
   + xorout[66] - xorout[67] >= -2
 R2744: sout[75] - sout[83] + sout[91] - sout[99] + sout[107] + sout[115]
   + xorout[66] - xorout[67] >= -2
 R2745: - sout[75] + sout[83] + sout[91] - sout[99] + sout[107] + sout[115]
   + xorout[66] - xorout[67] >= -2
 R2746: sout[75] - sout[83] - sout[91] + sout[99] + sout[107] + sout[115]
   + xorout[66] - xorout[67] >= -2
 R2747: - sout[75] + sout[83] - sout[91] + sout[99] + sout[107] + sout[115]
   + xorout[66] - xorout[67] >= -2
 R2748: - sout[75] - sout[83] + sout[91] + sout[99] + sout[107] + sout[115]
   + xorout[66] - xorout[67] >= -2
 R2749: sout[77] + sout[85] + sout[93] + sout[101] + sout[109] - sout[117]
   + xorout[69] >= 0
 R2750: sout[77] + sout[85] + sout[93] + sout[101] - sout[109] + sout[117]
   + xorout[69] >= 0
 R2751: sout[77] + sout[85] + sout[93] - sout[101] + sout[109] + sout[117]
   + xorout[69] >= 0
 R2752: sout[77] + sout[85] - sout[93] + sout[101] + sout[109] + sout[117]
   + xorout[69] >= 0
 R2753: sout[77] - sout[85] + sout[93] + sout[101] + sout[109] + sout[117]
   + xorout[69] >= 0
 R2754: - sout[77] + sout[85] + sout[93] + sout[101] + sout[109]
   + sout[117] + xorout[69] >= 0
 R2755: sout[77] + sout[85] + sout[93] + sout[101] + sout[109] + sout[117]
   - xorout[69] >= 0
 R2756: sout[77] + sout[85] + sout[93] - sout[101] - sout[109] - sout[117]
   + xorout[69] >= -2
 R2757: sout[77] + sout[85] - sout[93] + sout[101] - sout[109] - sout[117]
   + xorout[69] >= -2
 R2758: sout[77] - sout[85] + sout[93] + sout[101] - sout[109] - sout[117]
   + xorout[69] >= -2
 R2759: - sout[77] + sout[85] + sout[93] + sout[101] - sout[109]
   - sout[117] + xorout[69] >= -2
 R2760: sout[77] + sout[85] - sout[93] - sout[101] + sout[109] - sout[117]
   + xorout[69] >= -2
 R2761: sout[77] - sout[85] + sout[93] - sout[101] + sout[109] - sout[117]
   + xorout[69] >= -2
 R2762: - sout[77] + sout[85] + sout[93] - sout[101] + sout[109]
   - sout[117] + xorout[69] >= -2
 R2763: sout[77] - sout[85] - sout[93] + sout[101] + sout[109] - sout[117]
   + xorout[69] >= -2
 R2764: - sout[77] + sout[85] - sout[93] + sout[101] + sout[109]
   - sout[117] + xorout[69] >= -2
 R2765: - sout[77] - sout[85] + sout[93] + sout[101] + sout[109]
   - sout[117] + xorout[69] >= -2
 R2766: sout[77] + sout[85] - sout[93] - sout[101] - sout[109] + sout[117]
   + xorout[69] >= -2
 R2767: sout[77] - sout[85] + sout[93] - sout[101] - sout[109] + sout[117]
   + xorout[69] >= -2
 R2768: - sout[77] + sout[85] + sout[93] - sout[101] - sout[109]
   + sout[117] + xorout[69] >= -2
 R2769: sout[77] - sout[85] - sout[93] + sout[101] - sout[109] + sout[117]
   + xorout[69] >= -2
 R2770: - sout[77] + sout[85] - sout[93] + sout[101] - sout[109]
   + sout[117] + xorout[69] >= -2
 R2771: - sout[77] - sout[85] + sout[93] + sout[101] - sout[109]
   + sout[117] + xorout[69] >= -2
 R2772: sout[77] - sout[85] - sout[93] - sout[101] + sout[109] + sout[117]
   + xorout[69] >= -2
 R2773: - sout[77] + sout[85] - sout[93] - sout[101] + sout[109]
   + sout[117] + xorout[69] >= -2
 R2774: - sout[77] - sout[85] + sout[93] - sout[101] + sout[109]
   + sout[117] + xorout[69] >= -2
 R2775: - sout[77] - sout[85] - sout[93] + sout[101] + sout[109]
   + sout[117] + xorout[69] >= -2
 R2776: sout[77] - sout[85] - sout[93] - sout[101] - sout[109] - sout[117]
   + xorout[69] >= -4
 R2777: - sout[77] + sout[85] - sout[93] - sout[101] - sout[109]
   - sout[117] + xorout[69] >= -4
 R2778: - sout[77] - sout[85] + sout[93] - sout[101] - sout[109]
   - sout[117] + xorout[69] >= -4
 R2779: - sout[77] - sout[85] - sout[93] + sout[101] - sout[109]
   - sout[117] + xorout[69] >= -4
 R2780: - sout[77] - sout[85] - sout[93] - sout[101] + sout[109]
   - sout[117] + xorout[69] >= -4
 R2781: - sout[77] - sout[85] - sout[93] - sout[101] - sout[109]
   + sout[117] + xorout[69] >= -4
 R2782: - sout[116] + xorout[68] >= 0
 R2783: - sout[108] + xorout[68] >= 0
 R2784: - sout[100] + xorout[68] >= 0
 R2785: - sout[92] + xorout[68] >= 0
 R2786: - sout[84] + xorout[68] >= 0
 R2787: - sout[76] + xorout[68] >= 0
 R2788: - sout[76] + sout[77] >= 0
 R2789: - sout[84] + sout[85] >= 0
 R2790: - sout[92] + sout[93] >= 0
 R2791: - sout[100] + sout[101] >= 0
 R2792: - sout[108] + sout[109] >= 0
 R2793: - sout[116] + sout[117] >= 0
 R2794: sout[76] + sout[84] + sout[92] + sout[100] + sout[108] + sout[116]
   - xorout[68] >= 0
 R2795: - xorout[68] + xorout[69] >= 0
 R2796: - sout[77] - sout[85] - sout[93] - sout[101] - sout[109]
   - sout[117] + xorout[68] - xorout[69] >= -6
 R2797: sout[77] + sout[85] - sout[93] - sout[101] - sout[109] - sout[117]
   + xorout[68] - xorout[69] >= -4
 R2798: sout[77] - sout[85] + sout[93] - sout[101] - sout[109] - sout[117]
   + xorout[68] - xorout[69] >= -4
 R2799: - sout[77] + sout[85] + sout[93] - sout[101] - sout[109]
   - sout[117] + xorout[68] - xorout[69] >= -4
 R2800: sout[77] - sout[85] - sout[93] + sout[101] - sout[109] - sout[117]
   + xorout[68] - xorout[69] >= -4
 R2801: - sout[77] + sout[85] - sout[93] + sout[101] - sout[109]
   - sout[117] + xorout[68] - xorout[69] >= -4
 R2802: - sout[77] - sout[85] + sout[93] + sout[101] - sout[109]
   - sout[117] + xorout[68] - xorout[69] >= -4
 R2803: sout[77] - sout[85] - sout[93] - sout[101] + sout[109] - sout[117]
   + xorout[68] - xorout[69] >= -4
 R2804: - sout[77] + sout[85] - sout[93] - sout[101] + sout[109]
   - sout[117] + xorout[68] - xorout[69] >= -4
 R2805: - sout[77] - sout[85] + sout[93] - sout[101] + sout[109]
   - sout[117] + xorout[68] - xorout[69] >= -4
 R2806: - sout[77] - sout[85] - sout[93] + sout[101] + sout[109]
   - sout[117] + xorout[68] - xorout[69] >= -4
 R2807: sout[77] - sout[85] - sout[93] - sout[101] - sout[109] + sout[117]
   + xorout[68] - xorout[69] >= -4
 R2808: - sout[77] + sout[85] - sout[93] - sout[101] - sout[109]
   + sout[117] + xorout[68] - xorout[69] >= -4
 R2809: - sout[77] - sout[85] + sout[93] - sout[101] - sout[109]
   + sout[117] + xorout[68] - xorout[69] >= -4
 R2810: - sout[77] - sout[85] - sout[93] + sout[101] - sout[109]
   + sout[117] + xorout[68] - xorout[69] >= -4
 R2811: - sout[77] - sout[85] - sout[93] - sout[101] + sout[109]
   + sout[117] + xorout[68] - xorout[69] >= -4
 R2812: sout[77] + sout[85] + sout[93] + sout[101] - sout[109] - sout[117]
   + xorout[68] - xorout[69] >= -2
 R2813: sout[77] + sout[85] + sout[93] - sout[101] + sout[109] - sout[117]
   + xorout[68] - xorout[69] >= -2
 R2814: sout[77] + sout[85] - sout[93] + sout[101] + sout[109] - sout[117]
   + xorout[68] - xorout[69] >= -2
 R2815: sout[77] - sout[85] + sout[93] + sout[101] + sout[109] - sout[117]
   + xorout[68] - xorout[69] >= -2
 R2816: - sout[77] + sout[85] + sout[93] + sout[101] + sout[109]
   - sout[117] + xorout[68] - xorout[69] >= -2
 R2817: sout[77] + sout[85] + sout[93] - sout[101] - sout[109] + sout[117]
   + xorout[68] - xorout[69] >= -2
 R2818: sout[77] + sout[85] - sout[93] + sout[101] - sout[109] + sout[117]
   + xorout[68] - xorout[69] >= -2
 R2819: sout[77] - sout[85] + sout[93] + sout[101] - sout[109] + sout[117]
   + xorout[68] - xorout[69] >= -2
 R2820: - sout[77] + sout[85] + sout[93] + sout[101] - sout[109]
   + sout[117] + xorout[68] - xorout[69] >= -2
 R2821: sout[77] + sout[85] - sout[93] - sout[101] + sout[109] + sout[117]
   + xorout[68] - xorout[69] >= -2
 R2822: sout[77] - sout[85] + sout[93] - sout[101] + sout[109] + sout[117]
   + xorout[68] - xorout[69] >= -2
 R2823: - sout[77] + sout[85] + sout[93] - sout[101] + sout[109]
   + sout[117] + xorout[68] - xorout[69] >= -2
 R2824: sout[77] - sout[85] - sout[93] + sout[101] + sout[109] + sout[117]
   + xorout[68] - xorout[69] >= -2
 R2825: - sout[77] + sout[85] - sout[93] + sout[101] + sout[109]
   + sout[117] + xorout[68] - xorout[69] >= -2
 R2826: - sout[77] - sout[85] + sout[93] + sout[101] + sout[109]
   + sout[117] + xorout[68] - xorout[69] >= -2
 R2827: sout[79] + sout[87] + sout[95] + sout[103] + sout[111] - sout[119]
   + xorout[71] >= 0
 R2828: sout[79] + sout[87] + sout[95] + sout[103] - sout[111] + sout[119]
   + xorout[71] >= 0
 R2829: sout[79] + sout[87] + sout[95] - sout[103] + sout[111] + sout[119]
   + xorout[71] >= 0
 R2830: sout[79] + sout[87] - sout[95] + sout[103] + sout[111] + sout[119]
   + xorout[71] >= 0
 R2831: sout[79] - sout[87] + sout[95] + sout[103] + sout[111] + sout[119]
   + xorout[71] >= 0
 R2832: - sout[79] + sout[87] + sout[95] + sout[103] + sout[111]
   + sout[119] + xorout[71] >= 0
 R2833: sout[79] + sout[87] + sout[95] + sout[103] + sout[111] + sout[119]
   - xorout[71] >= 0
 R2834: sout[79] + sout[87] + sout[95] - sout[103] - sout[111] - sout[119]
   + xorout[71] >= -2
 R2835: sout[79] + sout[87] - sout[95] + sout[103] - sout[111] - sout[119]
   + xorout[71] >= -2
 R2836: sout[79] - sout[87] + sout[95] + sout[103] - sout[111] - sout[119]
   + xorout[71] >= -2
 R2837: - sout[79] + sout[87] + sout[95] + sout[103] - sout[111]
   - sout[119] + xorout[71] >= -2
 R2838: sout[79] + sout[87] - sout[95] - sout[103] + sout[111] - sout[119]
   + xorout[71] >= -2
 R2839: sout[79] - sout[87] + sout[95] - sout[103] + sout[111] - sout[119]
   + xorout[71] >= -2
 R2840: - sout[79] + sout[87] + sout[95] - sout[103] + sout[111]
   - sout[119] + xorout[71] >= -2
 R2841: sout[79] - sout[87] - sout[95] + sout[103] + sout[111] - sout[119]
   + xorout[71] >= -2
 R2842: - sout[79] + sout[87] - sout[95] + sout[103] + sout[111]
   - sout[119] + xorout[71] >= -2
 R2843: - sout[79] - sout[87] + sout[95] + sout[103] + sout[111]
   - sout[119] + xorout[71] >= -2
 R2844: sout[79] + sout[87] - sout[95] - sout[103] - sout[111] + sout[119]
   + xorout[71] >= -2
 R2845: sout[79] - sout[87] + sout[95] - sout[103] - sout[111] + sout[119]
   + xorout[71] >= -2
 R2846: - sout[79] + sout[87] + sout[95] - sout[103] - sout[111]
   + sout[119] + xorout[71] >= -2
 R2847: sout[79] - sout[87] - sout[95] + sout[103] - sout[111] + sout[119]
   + xorout[71] >= -2
 R2848: - sout[79] + sout[87] - sout[95] + sout[103] - sout[111]
   + sout[119] + xorout[71] >= -2
 R2849: - sout[79] - sout[87] + sout[95] + sout[103] - sout[111]
   + sout[119] + xorout[71] >= -2
 R2850: sout[79] - sout[87] - sout[95] - sout[103] + sout[111] + sout[119]
   + xorout[71] >= -2
 R2851: - sout[79] + sout[87] - sout[95] - sout[103] + sout[111]
   + sout[119] + xorout[71] >= -2
 R2852: - sout[79] - sout[87] + sout[95] - sout[103] + sout[111]
   + sout[119] + xorout[71] >= -2
 R2853: - sout[79] - sout[87] - sout[95] + sout[103] + sout[111]
   + sout[119] + xorout[71] >= -2
 R2854: sout[79] - sout[87] - sout[95] - sout[103] - sout[111] - sout[119]
   + xorout[71] >= -4
 R2855: - sout[79] + sout[87] - sout[95] - sout[103] - sout[111]
   - sout[119] + xorout[71] >= -4
 R2856: - sout[79] - sout[87] + sout[95] - sout[103] - sout[111]
   - sout[119] + xorout[71] >= -4
 R2857: - sout[79] - sout[87] - sout[95] + sout[103] - sout[111]
   - sout[119] + xorout[71] >= -4
 R2858: - sout[79] - sout[87] - sout[95] - sout[103] + sout[111]
   - sout[119] + xorout[71] >= -4
 R2859: - sout[79] - sout[87] - sout[95] - sout[103] - sout[111]
   + sout[119] + xorout[71] >= -4
 R2860: - sout[118] + xorout[70] >= 0
 R2861: - sout[110] + xorout[70] >= 0
 R2862: - sout[102] + xorout[70] >= 0
 R2863: - sout[94] + xorout[70] >= 0
 R2864: - sout[86] + xorout[70] >= 0
 R2865: - sout[78] + xorout[70] >= 0
 R2866: - sout[78] + sout[79] >= 0
 R2867: - sout[86] + sout[87] >= 0
 R2868: - sout[94] + sout[95] >= 0
 R2869: - sout[102] + sout[103] >= 0
 R2870: - sout[110] + sout[111] >= 0
 R2871: - sout[118] + sout[119] >= 0
 R2872: sout[78] + sout[86] + sout[94] + sout[102] + sout[110] + sout[118]
   - xorout[70] >= 0
 R2873: - xorout[70] + xorout[71] >= 0
 R2874: - sout[79] - sout[87] - sout[95] - sout[103] - sout[111]
   - sout[119] + xorout[70] - xorout[71] >= -6
 R2875: sout[79] + sout[87] - sout[95] - sout[103] - sout[111] - sout[119]
   + xorout[70] - xorout[71] >= -4
 R2876: sout[79] - sout[87] + sout[95] - sout[103] - sout[111] - sout[119]
   + xorout[70] - xorout[71] >= -4
 R2877: - sout[79] + sout[87] + sout[95] - sout[103] - sout[111]
   - sout[119] + xorout[70] - xorout[71] >= -4
 R2878: sout[79] - sout[87] - sout[95] + sout[103] - sout[111] - sout[119]
   + xorout[70] - xorout[71] >= -4
 R2879: - sout[79] + sout[87] - sout[95] + sout[103] - sout[111]
   - sout[119] + xorout[70] - xorout[71] >= -4
 R2880: - sout[79] - sout[87] + sout[95] + sout[103] - sout[111]
   - sout[119] + xorout[70] - xorout[71] >= -4
 R2881: sout[79] - sout[87] - sout[95] - sout[103] + sout[111] - sout[119]
   + xorout[70] - xorout[71] >= -4
 R2882: - sout[79] + sout[87] - sout[95] - sout[103] + sout[111]
   - sout[119] + xorout[70] - xorout[71] >= -4
 R2883: - sout[79] - sout[87] + sout[95] - sout[103] + sout[111]
   - sout[119] + xorout[70] - xorout[71] >= -4
 R2884: - sout[79] - sout[87] - sout[95] + sout[103] + sout[111]
   - sout[119] + xorout[70] - xorout[71] >= -4
 R2885: sout[79] - sout[87] - sout[95] - sout[103] - sout[111] + sout[119]
   + xorout[70] - xorout[71] >= -4
 R2886: - sout[79] + sout[87] - sout[95] - sout[103] - sout[111]
   + sout[119] + xorout[70] - xorout[71] >= -4
 R2887: - sout[79] - sout[87] + sout[95] - sout[103] - sout[111]
   + sout[119] + xorout[70] - xorout[71] >= -4
 R2888: - sout[79] - sout[87] - sout[95] + sout[103] - sout[111]
   + sout[119] + xorout[70] - xorout[71] >= -4
 R2889: - sout[79] - sout[87] - sout[95] - sout[103] + sout[111]
   + sout[119] + xorout[70] - xorout[71] >= -4
 R2890: sout[79] + sout[87] + sout[95] + sout[103] - sout[111] - sout[119]
   + xorout[70] - xorout[71] >= -2
 R2891: sout[79] + sout[87] + sout[95] - sout[103] + sout[111] - sout[119]
   + xorout[70] - xorout[71] >= -2
 R2892: sout[79] + sout[87] - sout[95] + sout[103] + sout[111] - sout[119]
   + xorout[70] - xorout[71] >= -2
 R2893: sout[79] - sout[87] + sout[95] + sout[103] + sout[111] - sout[119]
   + xorout[70] - xorout[71] >= -2
 R2894: - sout[79] + sout[87] + sout[95] + sout[103] + sout[111]
   - sout[119] + xorout[70] - xorout[71] >= -2
 R2895: sout[79] + sout[87] + sout[95] - sout[103] - sout[111] + sout[119]
   + xorout[70] - xorout[71] >= -2
 R2896: sout[79] + sout[87] - sout[95] + sout[103] - sout[111] + sout[119]
   + xorout[70] - xorout[71] >= -2
 R2897: sout[79] - sout[87] + sout[95] + sout[103] - sout[111] + sout[119]
   + xorout[70] - xorout[71] >= -2
 R2898: - sout[79] + sout[87] + sout[95] + sout[103] - sout[111]
   + sout[119] + xorout[70] - xorout[71] >= -2
 R2899: sout[79] + sout[87] - sout[95] - sout[103] + sout[111] + sout[119]
   + xorout[70] - xorout[71] >= -2
 R2900: sout[79] - sout[87] + sout[95] - sout[103] + sout[111] + sout[119]
   + xorout[70] - xorout[71] >= -2
 R2901: - sout[79] + sout[87] + sout[95] - sout[103] + sout[111]
   + sout[119] + xorout[70] - xorout[71] >= -2
 R2902: sout[79] - sout[87] - sout[95] + sout[103] + sout[111] + sout[119]
   + xorout[70] - xorout[71] >= -2
 R2903: - sout[79] + sout[87] - sout[95] + sout[103] + sout[111]
   + sout[119] + xorout[70] - xorout[71] >= -2
 R2904: - sout[79] - sout[87] + sout[95] + sout[103] + sout[111]
   + sout[119] + xorout[70] - xorout[71] >= -2
 R2905: sout[65] + sout[81] + sout[89] + sout[105] + sout[113] - sout[121]
   + xorout[73] >= 0
 R2906: sout[65] + sout[81] + sout[89] + sout[105] - sout[113] + sout[121]
   + xorout[73] >= 0
 R2907: sout[65] + sout[81] + sout[89] - sout[105] + sout[113] + sout[121]
   + xorout[73] >= 0
 R2908: sout[65] + sout[81] - sout[89] + sout[105] + sout[113] + sout[121]
   + xorout[73] >= 0
 R2909: sout[65] - sout[81] + sout[89] + sout[105] + sout[113] + sout[121]
   + xorout[73] >= 0
 R2910: - sout[65] + sout[81] + sout[89] + sout[105] + sout[113]
   + sout[121] + xorout[73] >= 0
 R2911: sout[65] + sout[81] + sout[89] + sout[105] + sout[113] + sout[121]
   - xorout[73] >= 0
 R2912: sout[65] + sout[81] + sout[89] - sout[105] - sout[113] - sout[121]
   + xorout[73] >= -2
 R2913: sout[65] + sout[81] - sout[89] + sout[105] - sout[113] - sout[121]
   + xorout[73] >= -2
 R2914: sout[65] - sout[81] + sout[89] + sout[105] - sout[113] - sout[121]
   + xorout[73] >= -2
 R2915: - sout[65] + sout[81] + sout[89] + sout[105] - sout[113]
   - sout[121] + xorout[73] >= -2
 R2916: sout[65] + sout[81] - sout[89] - sout[105] + sout[113] - sout[121]
   + xorout[73] >= -2
 R2917: sout[65] - sout[81] + sout[89] - sout[105] + sout[113] - sout[121]
   + xorout[73] >= -2
 R2918: - sout[65] + sout[81] + sout[89] - sout[105] + sout[113]
   - sout[121] + xorout[73] >= -2
 R2919: sout[65] - sout[81] - sout[89] + sout[105] + sout[113] - sout[121]
   + xorout[73] >= -2
 R2920: - sout[65] + sout[81] - sout[89] + sout[105] + sout[113]
   - sout[121] + xorout[73] >= -2
 R2921: - sout[65] - sout[81] + sout[89] + sout[105] + sout[113]
   - sout[121] + xorout[73] >= -2
 R2922: sout[65] + sout[81] - sout[89] - sout[105] - sout[113] + sout[121]
   + xorout[73] >= -2
 R2923: sout[65] - sout[81] + sout[89] - sout[105] - sout[113] + sout[121]
   + xorout[73] >= -2
 R2924: - sout[65] + sout[81] + sout[89] - sout[105] - sout[113]
   + sout[121] + xorout[73] >= -2
 R2925: sout[65] - sout[81] - sout[89] + sout[105] - sout[113] + sout[121]
   + xorout[73] >= -2
 R2926: - sout[65] + sout[81] - sout[89] + sout[105] - sout[113]
   + sout[121] + xorout[73] >= -2
 R2927: - sout[65] - sout[81] + sout[89] + sout[105] - sout[113]
   + sout[121] + xorout[73] >= -2
 R2928: sout[65] - sout[81] - sout[89] - sout[105] + sout[113] + sout[121]
   + xorout[73] >= -2
 R2929: - sout[65] + sout[81] - sout[89] - sout[105] + sout[113]
   + sout[121] + xorout[73] >= -2
 R2930: - sout[65] - sout[81] + sout[89] - sout[105] + sout[113]
   + sout[121] + xorout[73] >= -2
 R2931: - sout[65] - sout[81] - sout[89] + sout[105] + sout[113]
   + sout[121] + xorout[73] >= -2
 R2932: sout[65] - sout[81] - sout[89] - sout[105] - sout[113] - sout[121]
   + xorout[73] >= -4
 R2933: - sout[65] + sout[81] - sout[89] - sout[105] - sout[113]
   - sout[121] + xorout[73] >= -4
 R2934: - sout[65] - sout[81] + sout[89] - sout[105] - sout[113]
   - sout[121] + xorout[73] >= -4
 R2935: - sout[65] - sout[81] - sout[89] + sout[105] - sout[113]
   - sout[121] + xorout[73] >= -4
 R2936: - sout[65] - sout[81] - sout[89] - sout[105] + sout[113]
   - sout[121] + xorout[73] >= -4
 R2937: - sout[65] - sout[81] - sout[89] - sout[105] - sout[113]
   + sout[121] + xorout[73] >= -4
 R2938: - sout[120] + xorout[72] >= 0
 R2939: - sout[112] + xorout[72] >= 0
 R2940: - sout[104] + xorout[72] >= 0
 R2941: - sout[88] + xorout[72] >= 0
 R2942: - sout[80] + xorout[72] >= 0
 R2943: - sout[64] + xorout[72] >= 0
 R2944: - sout[64] + sout[65] >= 0
 R2945: - sout[80] + sout[81] >= 0
 R2946: - sout[88] + sout[89] >= 0
 R2947: - sout[104] + sout[105] >= 0
 R2948: - sout[112] + sout[113] >= 0
 R2949: - sout[120] + sout[121] >= 0
 R2950: sout[64] + sout[80] + sout[88] + sout[104] + sout[112] + sout[120]
   - xorout[72] >= 0
 R2951: - xorout[72] + xorout[73] >= 0
 R2952: - sout[65] - sout[81] - sout[89] - sout[105] - sout[113]
   - sout[121] + xorout[72] - xorout[73] >= -6
 R2953: sout[65] + sout[81] - sout[89] - sout[105] - sout[113] - sout[121]
   + xorout[72] - xorout[73] >= -4
 R2954: sout[65] - sout[81] + sout[89] - sout[105] - sout[113] - sout[121]
   + xorout[72] - xorout[73] >= -4
 R2955: - sout[65] + sout[81] + sout[89] - sout[105] - sout[113]
   - sout[121] + xorout[72] - xorout[73] >= -4
 R2956: sout[65] - sout[81] - sout[89] + sout[105] - sout[113] - sout[121]
   + xorout[72] - xorout[73] >= -4
 R2957: - sout[65] + sout[81] - sout[89] + sout[105] - sout[113]
   - sout[121] + xorout[72] - xorout[73] >= -4
 R2958: - sout[65] - sout[81] + sout[89] + sout[105] - sout[113]
   - sout[121] + xorout[72] - xorout[73] >= -4
 R2959: sout[65] - sout[81] - sout[89] - sout[105] + sout[113] - sout[121]
   + xorout[72] - xorout[73] >= -4
 R2960: - sout[65] + sout[81] - sout[89] - sout[105] + sout[113]
   - sout[121] + xorout[72] - xorout[73] >= -4
 R2961: - sout[65] - sout[81] + sout[89] - sout[105] + sout[113]
   - sout[121] + xorout[72] - xorout[73] >= -4
 R2962: - sout[65] - sout[81] - sout[89] + sout[105] + sout[113]
   - sout[121] + xorout[72] - xorout[73] >= -4
 R2963: sout[65] - sout[81] - sout[89] - sout[105] - sout[113] + sout[121]
   + xorout[72] - xorout[73] >= -4
 R2964: - sout[65] + sout[81] - sout[89] - sout[105] - sout[113]
   + sout[121] + xorout[72] - xorout[73] >= -4
 R2965: - sout[65] - sout[81] + sout[89] - sout[105] - sout[113]
   + sout[121] + xorout[72] - xorout[73] >= -4
 R2966: - sout[65] - sout[81] - sout[89] + sout[105] - sout[113]
   + sout[121] + xorout[72] - xorout[73] >= -4
 R2967: - sout[65] - sout[81] - sout[89] - sout[105] + sout[113]
   + sout[121] + xorout[72] - xorout[73] >= -4
 R2968: sout[65] + sout[81] + sout[89] + sout[105] - sout[113] - sout[121]
   + xorout[72] - xorout[73] >= -2
 R2969: sout[65] + sout[81] + sout[89] - sout[105] + sout[113] - sout[121]
   + xorout[72] - xorout[73] >= -2
 R2970: sout[65] + sout[81] - sout[89] + sout[105] + sout[113] - sout[121]
   + xorout[72] - xorout[73] >= -2
 R2971: sout[65] - sout[81] + sout[89] + sout[105] + sout[113] - sout[121]
   + xorout[72] - xorout[73] >= -2
 R2972: - sout[65] + sout[81] + sout[89] + sout[105] + sout[113]
   - sout[121] + xorout[72] - xorout[73] >= -2
 R2973: sout[65] + sout[81] + sout[89] - sout[105] - sout[113] + sout[121]
   + xorout[72] - xorout[73] >= -2
 R2974: sout[65] + sout[81] - sout[89] + sout[105] - sout[113] + sout[121]
   + xorout[72] - xorout[73] >= -2
 R2975: sout[65] - sout[81] + sout[89] + sout[105] - sout[113] + sout[121]
   + xorout[72] - xorout[73] >= -2
 R2976: - sout[65] + sout[81] + sout[89] + sout[105] - sout[113]
   + sout[121] + xorout[72] - xorout[73] >= -2
 R2977: sout[65] + sout[81] - sout[89] - sout[105] + sout[113] + sout[121]
   + xorout[72] - xorout[73] >= -2
 R2978: sout[65] - sout[81] + sout[89] - sout[105] + sout[113] + sout[121]
   + xorout[72] - xorout[73] >= -2
 R2979: - sout[65] + sout[81] + sout[89] - sout[105] + sout[113]
   + sout[121] + xorout[72] - xorout[73] >= -2
 R2980: sout[65] - sout[81] - sout[89] + sout[105] + sout[113] + sout[121]
   + xorout[72] - xorout[73] >= -2
 R2981: - sout[65] + sout[81] - sout[89] + sout[105] + sout[113]
   + sout[121] + xorout[72] - xorout[73] >= -2
 R2982: - sout[65] - sout[81] + sout[89] + sout[105] + sout[113]
   + sout[121] + xorout[72] - xorout[73] >= -2
 R2983: sout[67] + sout[83] + sout[91] + sout[107] + sout[115] - sout[123]
   + xorout[75] >= 0
 R2984: sout[67] + sout[83] + sout[91] + sout[107] - sout[115] + sout[123]
   + xorout[75] >= 0
 R2985: sout[67] + sout[83] + sout[91] - sout[107] + sout[115] + sout[123]
   + xorout[75] >= 0
 R2986: sout[67] + sout[83] - sout[91] + sout[107] + sout[115] + sout[123]
   + xorout[75] >= 0
 R2987: sout[67] - sout[83] + sout[91] + sout[107] + sout[115] + sout[123]
   + xorout[75] >= 0
 R2988: - sout[67] + sout[83] + sout[91] + sout[107] + sout[115]
   + sout[123] + xorout[75] >= 0
 R2989: sout[67] + sout[83] + sout[91] + sout[107] + sout[115] + sout[123]
   - xorout[75] >= 0
 R2990: sout[67] + sout[83] + sout[91] - sout[107] - sout[115] - sout[123]
   + xorout[75] >= -2
 R2991: sout[67] + sout[83] - sout[91] + sout[107] - sout[115] - sout[123]
   + xorout[75] >= -2
 R2992: sout[67] - sout[83] + sout[91] + sout[107] - sout[115] - sout[123]
   + xorout[75] >= -2
 R2993: - sout[67] + sout[83] + sout[91] + sout[107] - sout[115]
   - sout[123] + xorout[75] >= -2
 R2994: sout[67] + sout[83] - sout[91] - sout[107] + sout[115] - sout[123]
   + xorout[75] >= -2
 R2995: sout[67] - sout[83] + sout[91] - sout[107] + sout[115] - sout[123]
   + xorout[75] >= -2
 R2996: - sout[67] + sout[83] + sout[91] - sout[107] + sout[115]
   - sout[123] + xorout[75] >= -2
 R2997: sout[67] - sout[83] - sout[91] + sout[107] + sout[115] - sout[123]
   + xorout[75] >= -2
 R2998: - sout[67] + sout[83] - sout[91] + sout[107] + sout[115]
   - sout[123] + xorout[75] >= -2
 R2999: - sout[67] - sout[83] + sout[91] + sout[107] + sout[115]
   - sout[123] + xorout[75] >= -2
 R3000: sout[67] + sout[83] - sout[91] - sout[107] - sout[115] + sout[123]
   + xorout[75] >= -2
 R3001: sout[67] - sout[83] + sout[91] - sout[107] - sout[115] + sout[123]
   + xorout[75] >= -2
 R3002: - sout[67] + sout[83] + sout[91] - sout[107] - sout[115]
   + sout[123] + xorout[75] >= -2
 R3003: sout[67] - sout[83] - sout[91] + sout[107] - sout[115] + sout[123]
   + xorout[75] >= -2
 R3004: - sout[67] + sout[83] - sout[91] + sout[107] - sout[115]
   + sout[123] + xorout[75] >= -2
 R3005: - sout[67] - sout[83] + sout[91] + sout[107] - sout[115]
   + sout[123] + xorout[75] >= -2
 R3006: sout[67] - sout[83] - sout[91] - sout[107] + sout[115] + sout[123]
   + xorout[75] >= -2
 R3007: - sout[67] + sout[83] - sout[91] - sout[107] + sout[115]
   + sout[123] + xorout[75] >= -2
 R3008: - sout[67] - sout[83] + sout[91] - sout[107] + sout[115]
   + sout[123] + xorout[75] >= -2
 R3009: - sout[67] - sout[83] - sout[91] + sout[107] + sout[115]
   + sout[123] + xorout[75] >= -2
 R3010: sout[67] - sout[83] - sout[91] - sout[107] - sout[115] - sout[123]
   + xorout[75] >= -4
 R3011: - sout[67] + sout[83] - sout[91] - sout[107] - sout[115]
   - sout[123] + xorout[75] >= -4
 R3012: - sout[67] - sout[83] + sout[91] - sout[107] - sout[115]
   - sout[123] + xorout[75] >= -4
 R3013: - sout[67] - sout[83] - sout[91] + sout[107] - sout[115]
   - sout[123] + xorout[75] >= -4
 R3014: - sout[67] - sout[83] - sout[91] - sout[107] + sout[115]
   - sout[123] + xorout[75] >= -4
 R3015: - sout[67] - sout[83] - sout[91] - sout[107] - sout[115]
   + sout[123] + xorout[75] >= -4
 R3016: - sout[122] + xorout[74] >= 0
 R3017: - sout[114] + xorout[74] >= 0
 R3018: - sout[106] + xorout[74] >= 0
 R3019: - sout[90] + xorout[74] >= 0
 R3020: - sout[82] + xorout[74] >= 0
 R3021: - sout[66] + xorout[74] >= 0
 R3022: - sout[66] + sout[67] >= 0
 R3023: - sout[82] + sout[83] >= 0
 R3024: - sout[90] + sout[91] >= 0
 R3025: - sout[106] + sout[107] >= 0
 R3026: - sout[114] + sout[115] >= 0
 R3027: - sout[122] + sout[123] >= 0
 R3028: sout[66] + sout[82] + sout[90] + sout[106] + sout[114] + sout[122]
   - xorout[74] >= 0
 R3029: - xorout[74] + xorout[75] >= 0
 R3030: - sout[67] - sout[83] - sout[91] - sout[107] - sout[115]
   - sout[123] + xorout[74] - xorout[75] >= -6
 R3031: sout[67] + sout[83] - sout[91] - sout[107] - sout[115] - sout[123]
   + xorout[74] - xorout[75] >= -4
 R3032: sout[67] - sout[83] + sout[91] - sout[107] - sout[115] - sout[123]
   + xorout[74] - xorout[75] >= -4
 R3033: - sout[67] + sout[83] + sout[91] - sout[107] - sout[115]
   - sout[123] + xorout[74] - xorout[75] >= -4
 R3034: sout[67] - sout[83] - sout[91] + sout[107] - sout[115] - sout[123]
   + xorout[74] - xorout[75] >= -4
 R3035: - sout[67] + sout[83] - sout[91] + sout[107] - sout[115]
   - sout[123] + xorout[74] - xorout[75] >= -4
 R3036: - sout[67] - sout[83] + sout[91] + sout[107] - sout[115]
   - sout[123] + xorout[74] - xorout[75] >= -4
 R3037: sout[67] - sout[83] - sout[91] - sout[107] + sout[115] - sout[123]
   + xorout[74] - xorout[75] >= -4
 R3038: - sout[67] + sout[83] - sout[91] - sout[107] + sout[115]
   - sout[123] + xorout[74] - xorout[75] >= -4
 R3039: - sout[67] - sout[83] + sout[91] - sout[107] + sout[115]
   - sout[123] + xorout[74] - xorout[75] >= -4
 R3040: - sout[67] - sout[83] - sout[91] + sout[107] + sout[115]
   - sout[123] + xorout[74] - xorout[75] >= -4
 R3041: sout[67] - sout[83] - sout[91] - sout[107] - sout[115] + sout[123]
   + xorout[74] - xorout[75] >= -4
 R3042: - sout[67] + sout[83] - sout[91] - sout[107] - sout[115]
   + sout[123] + xorout[74] - xorout[75] >= -4
 R3043: - sout[67] - sout[83] + sout[91] - sout[107] - sout[115]
   + sout[123] + xorout[74] - xorout[75] >= -4
 R3044: - sout[67] - sout[83] - sout[91] + sout[107] - sout[115]
   + sout[123] + xorout[74] - xorout[75] >= -4
 R3045: - sout[67] - sout[83] - sout[91] - sout[107] + sout[115]
   + sout[123] + xorout[74] - xorout[75] >= -4
 R3046: sout[67] + sout[83] + sout[91] + sout[107] - sout[115] - sout[123]
   + xorout[74] - xorout[75] >= -2
 R3047: sout[67] + sout[83] + sout[91] - sout[107] + sout[115] - sout[123]
   + xorout[74] - xorout[75] >= -2
 R3048: sout[67] + sout[83] - sout[91] + sout[107] + sout[115] - sout[123]
   + xorout[74] - xorout[75] >= -2
 R3049: sout[67] - sout[83] + sout[91] + sout[107] + sout[115] - sout[123]
   + xorout[74] - xorout[75] >= -2
 R3050: - sout[67] + sout[83] + sout[91] + sout[107] + sout[115]
   - sout[123] + xorout[74] - xorout[75] >= -2
 R3051: sout[67] + sout[83] + sout[91] - sout[107] - sout[115] + sout[123]
   + xorout[74] - xorout[75] >= -2
 R3052: sout[67] + sout[83] - sout[91] + sout[107] - sout[115] + sout[123]
   + xorout[74] - xorout[75] >= -2
 R3053: sout[67] - sout[83] + sout[91] + sout[107] - sout[115] + sout[123]
   + xorout[74] - xorout[75] >= -2
 R3054: - sout[67] + sout[83] + sout[91] + sout[107] - sout[115]
   + sout[123] + xorout[74] - xorout[75] >= -2
 R3055: sout[67] + sout[83] - sout[91] - sout[107] + sout[115] + sout[123]
   + xorout[74] - xorout[75] >= -2
 R3056: sout[67] - sout[83] + sout[91] - sout[107] + sout[115] + sout[123]
   + xorout[74] - xorout[75] >= -2
 R3057: - sout[67] + sout[83] + sout[91] - sout[107] + sout[115]
   + sout[123] + xorout[74] - xorout[75] >= -2
 R3058: sout[67] - sout[83] - sout[91] + sout[107] + sout[115] + sout[123]
   + xorout[74] - xorout[75] >= -2
 R3059: - sout[67] + sout[83] - sout[91] + sout[107] + sout[115]
   + sout[123] + xorout[74] - xorout[75] >= -2
 R3060: - sout[67] - sout[83] + sout[91] + sout[107] + sout[115]
   + sout[123] + xorout[74] - xorout[75] >= -2
 R3061: sout[69] + sout[85] + sout[93] + sout[109] + sout[117] - sout[125]
   + xorout[77] >= 0
 R3062: sout[69] + sout[85] + sout[93] + sout[109] - sout[117] + sout[125]
   + xorout[77] >= 0
 R3063: sout[69] + sout[85] + sout[93] - sout[109] + sout[117] + sout[125]
   + xorout[77] >= 0
 R3064: sout[69] + sout[85] - sout[93] + sout[109] + sout[117] + sout[125]
   + xorout[77] >= 0
 R3065: sout[69] - sout[85] + sout[93] + sout[109] + sout[117] + sout[125]
   + xorout[77] >= 0
 R3066: - sout[69] + sout[85] + sout[93] + sout[109] + sout[117]
   + sout[125] + xorout[77] >= 0
 R3067: sout[69] + sout[85] + sout[93] + sout[109] + sout[117] + sout[125]
   - xorout[77] >= 0
 R3068: sout[69] + sout[85] + sout[93] - sout[109] - sout[117] - sout[125]
   + xorout[77] >= -2
 R3069: sout[69] + sout[85] - sout[93] + sout[109] - sout[117] - sout[125]
   + xorout[77] >= -2
 R3070: sout[69] - sout[85] + sout[93] + sout[109] - sout[117] - sout[125]
   + xorout[77] >= -2
 R3071: - sout[69] + sout[85] + sout[93] + sout[109] - sout[117]
   - sout[125] + xorout[77] >= -2
 R3072: sout[69] + sout[85] - sout[93] - sout[109] + sout[117] - sout[125]
   + xorout[77] >= -2
 R3073: sout[69] - sout[85] + sout[93] - sout[109] + sout[117] - sout[125]
   + xorout[77] >= -2
 R3074: - sout[69] + sout[85] + sout[93] - sout[109] + sout[117]
   - sout[125] + xorout[77] >= -2
 R3075: sout[69] - sout[85] - sout[93] + sout[109] + sout[117] - sout[125]
   + xorout[77] >= -2
 R3076: - sout[69] + sout[85] - sout[93] + sout[109] + sout[117]
   - sout[125] + xorout[77] >= -2
 R3077: - sout[69] - sout[85] + sout[93] + sout[109] + sout[117]
   - sout[125] + xorout[77] >= -2
 R3078: sout[69] + sout[85] - sout[93] - sout[109] - sout[117] + sout[125]
   + xorout[77] >= -2
 R3079: sout[69] - sout[85] + sout[93] - sout[109] - sout[117] + sout[125]
   + xorout[77] >= -2
 R3080: - sout[69] + sout[85] + sout[93] - sout[109] - sout[117]
   + sout[125] + xorout[77] >= -2
 R3081: sout[69] - sout[85] - sout[93] + sout[109] - sout[117] + sout[125]
   + xorout[77] >= -2
 R3082: - sout[69] + sout[85] - sout[93] + sout[109] - sout[117]
   + sout[125] + xorout[77] >= -2
 R3083: - sout[69] - sout[85] + sout[93] + sout[109] - sout[117]
   + sout[125] + xorout[77] >= -2
 R3084: sout[69] - sout[85] - sout[93] - sout[109] + sout[117] + sout[125]
   + xorout[77] >= -2
 R3085: - sout[69] + sout[85] - sout[93] - sout[109] + sout[117]
   + sout[125] + xorout[77] >= -2
 R3086: - sout[69] - sout[85] + sout[93] - sout[109] + sout[117]
   + sout[125] + xorout[77] >= -2
 R3087: - sout[69] - sout[85] - sout[93] + sout[109] + sout[117]
   + sout[125] + xorout[77] >= -2
 R3088: sout[69] - sout[85] - sout[93] - sout[109] - sout[117] - sout[125]
   + xorout[77] >= -4
 R3089: - sout[69] + sout[85] - sout[93] - sout[109] - sout[117]
   - sout[125] + xorout[77] >= -4
 R3090: - sout[69] - sout[85] + sout[93] - sout[109] - sout[117]
   - sout[125] + xorout[77] >= -4
 R3091: - sout[69] - sout[85] - sout[93] + sout[109] - sout[117]
   - sout[125] + xorout[77] >= -4
 R3092: - sout[69] - sout[85] - sout[93] - sout[109] + sout[117]
   - sout[125] + xorout[77] >= -4
 R3093: - sout[69] - sout[85] - sout[93] - sout[109] - sout[117]
   + sout[125] + xorout[77] >= -4
 R3094: - sout[124] + xorout[76] >= 0
 R3095: - sout[116] + xorout[76] >= 0
 R3096: - sout[108] + xorout[76] >= 0
 R3097: - sout[92] + xorout[76] >= 0
 R3098: - sout[84] + xorout[76] >= 0
 R3099: - sout[68] + xorout[76] >= 0
 R3100: - sout[68] + sout[69] >= 0
 R3101: - sout[84] + sout[85] >= 0
 R3102: - sout[92] + sout[93] >= 0
 R3103: - sout[108] + sout[109] >= 0
 R3104: - sout[116] + sout[117] >= 0
 R3105: - sout[124] + sout[125] >= 0
 R3106: sout[68] + sout[84] + sout[92] + sout[108] + sout[116] + sout[124]
   - xorout[76] >= 0
 R3107: - xorout[76] + xorout[77] >= 0
 R3108: - sout[69] - sout[85] - sout[93] - sout[109] - sout[117]
   - sout[125] + xorout[76] - xorout[77] >= -6
 R3109: sout[69] + sout[85] - sout[93] - sout[109] - sout[117] - sout[125]
   + xorout[76] - xorout[77] >= -4
 R3110: sout[69] - sout[85] + sout[93] - sout[109] - sout[117] - sout[125]
   + xorout[76] - xorout[77] >= -4
 R3111: - sout[69] + sout[85] + sout[93] - sout[109] - sout[117]
   - sout[125] + xorout[76] - xorout[77] >= -4
 R3112: sout[69] - sout[85] - sout[93] + sout[109] - sout[117] - sout[125]
   + xorout[76] - xorout[77] >= -4
 R3113: - sout[69] + sout[85] - sout[93] + sout[109] - sout[117]
   - sout[125] + xorout[76] - xorout[77] >= -4
 R3114: - sout[69] - sout[85] + sout[93] + sout[109] - sout[117]
   - sout[125] + xorout[76] - xorout[77] >= -4
 R3115: sout[69] - sout[85] - sout[93] - sout[109] + sout[117] - sout[125]
   + xorout[76] - xorout[77] >= -4
 R3116: - sout[69] + sout[85] - sout[93] - sout[109] + sout[117]
   - sout[125] + xorout[76] - xorout[77] >= -4
 R3117: - sout[69] - sout[85] + sout[93] - sout[109] + sout[117]
   - sout[125] + xorout[76] - xorout[77] >= -4
 R3118: - sout[69] - sout[85] - sout[93] + sout[109] + sout[117]
   - sout[125] + xorout[76] - xorout[77] >= -4
 R3119: sout[69] - sout[85] - sout[93] - sout[109] - sout[117] + sout[125]
   + xorout[76] - xorout[77] >= -4
 R3120: - sout[69] + sout[85] - sout[93] - sout[109] - sout[117]
   + sout[125] + xorout[76] - xorout[77] >= -4
 R3121: - sout[69] - sout[85] + sout[93] - sout[109] - sout[117]
   + sout[125] + xorout[76] - xorout[77] >= -4
 R3122: - sout[69] - sout[85] - sout[93] + sout[109] - sout[117]
   + sout[125] + xorout[76] - xorout[77] >= -4
 R3123: - sout[69] - sout[85] - sout[93] - sout[109] + sout[117]
   + sout[125] + xorout[76] - xorout[77] >= -4
 R3124: sout[69] + sout[85] + sout[93] + sout[109] - sout[117] - sout[125]
   + xorout[76] - xorout[77] >= -2
 R3125: sout[69] + sout[85] + sout[93] - sout[109] + sout[117] - sout[125]
   + xorout[76] - xorout[77] >= -2
 R3126: sout[69] + sout[85] - sout[93] + sout[109] + sout[117] - sout[125]
   + xorout[76] - xorout[77] >= -2
 R3127: sout[69] - sout[85] + sout[93] + sout[109] + sout[117] - sout[125]
   + xorout[76] - xorout[77] >= -2
 R3128: - sout[69] + sout[85] + sout[93] + sout[109] + sout[117]
   - sout[125] + xorout[76] - xorout[77] >= -2
 R3129: sout[69] + sout[85] + sout[93] - sout[109] - sout[117] + sout[125]
   + xorout[76] - xorout[77] >= -2
 R3130: sout[69] + sout[85] - sout[93] + sout[109] - sout[117] + sout[125]
   + xorout[76] - xorout[77] >= -2
 R3131: sout[69] - sout[85] + sout[93] + sout[109] - sout[117] + sout[125]
   + xorout[76] - xorout[77] >= -2
 R3132: - sout[69] + sout[85] + sout[93] + sout[109] - sout[117]
   + sout[125] + xorout[76] - xorout[77] >= -2
 R3133: sout[69] + sout[85] - sout[93] - sout[109] + sout[117] + sout[125]
   + xorout[76] - xorout[77] >= -2
 R3134: sout[69] - sout[85] + sout[93] - sout[109] + sout[117] + sout[125]
   + xorout[76] - xorout[77] >= -2
 R3135: - sout[69] + sout[85] + sout[93] - sout[109] + sout[117]
   + sout[125] + xorout[76] - xorout[77] >= -2
 R3136: sout[69] - sout[85] - sout[93] + sout[109] + sout[117] + sout[125]
   + xorout[76] - xorout[77] >= -2
 R3137: - sout[69] + sout[85] - sout[93] + sout[109] + sout[117]
   + sout[125] + xorout[76] - xorout[77] >= -2
 R3138: - sout[69] - sout[85] + sout[93] + sout[109] + sout[117]
   + sout[125] + xorout[76] - xorout[77] >= -2
 R3139: sout[71] + sout[87] + sout[95] + sout[111] + sout[119] - sout[127]
   + xorout[79] >= 0
 R3140: sout[71] + sout[87] + sout[95] + sout[111] - sout[119] + sout[127]
   + xorout[79] >= 0
 R3141: sout[71] + sout[87] + sout[95] - sout[111] + sout[119] + sout[127]
   + xorout[79] >= 0
 R3142: sout[71] + sout[87] - sout[95] + sout[111] + sout[119] + sout[127]
   + xorout[79] >= 0
 R3143: sout[71] - sout[87] + sout[95] + sout[111] + sout[119] + sout[127]
   + xorout[79] >= 0
 R3144: - sout[71] + sout[87] + sout[95] + sout[111] + sout[119]
   + sout[127] + xorout[79] >= 0
 R3145: sout[71] + sout[87] + sout[95] + sout[111] + sout[119] + sout[127]
   - xorout[79] >= 0
 R3146: sout[71] + sout[87] + sout[95] - sout[111] - sout[119] - sout[127]
   + xorout[79] >= -2
 R3147: sout[71] + sout[87] - sout[95] + sout[111] - sout[119] - sout[127]
   + xorout[79] >= -2
 R3148: sout[71] - sout[87] + sout[95] + sout[111] - sout[119] - sout[127]
   + xorout[79] >= -2
 R3149: - sout[71] + sout[87] + sout[95] + sout[111] - sout[119]
   - sout[127] + xorout[79] >= -2
 R3150: sout[71] + sout[87] - sout[95] - sout[111] + sout[119] - sout[127]
   + xorout[79] >= -2
 R3151: sout[71] - sout[87] + sout[95] - sout[111] + sout[119] - sout[127]
   + xorout[79] >= -2
 R3152: - sout[71] + sout[87] + sout[95] - sout[111] + sout[119]
   - sout[127] + xorout[79] >= -2
 R3153: sout[71] - sout[87] - sout[95] + sout[111] + sout[119] - sout[127]
   + xorout[79] >= -2
 R3154: - sout[71] + sout[87] - sout[95] + sout[111] + sout[119]
   - sout[127] + xorout[79] >= -2
 R3155: - sout[71] - sout[87] + sout[95] + sout[111] + sout[119]
   - sout[127] + xorout[79] >= -2
 R3156: sout[71] + sout[87] - sout[95] - sout[111] - sout[119] + sout[127]
   + xorout[79] >= -2
 R3157: sout[71] - sout[87] + sout[95] - sout[111] - sout[119] + sout[127]
   + xorout[79] >= -2
 R3158: - sout[71] + sout[87] + sout[95] - sout[111] - sout[119]
   + sout[127] + xorout[79] >= -2
 R3159: sout[71] - sout[87] - sout[95] + sout[111] - sout[119] + sout[127]
   + xorout[79] >= -2
 R3160: - sout[71] + sout[87] - sout[95] + sout[111] - sout[119]
   + sout[127] + xorout[79] >= -2
 R3161: - sout[71] - sout[87] + sout[95] + sout[111] - sout[119]
   + sout[127] + xorout[79] >= -2
 R3162: sout[71] - sout[87] - sout[95] - sout[111] + sout[119] + sout[127]
   + xorout[79] >= -2
 R3163: - sout[71] + sout[87] - sout[95] - sout[111] + sout[119]
   + sout[127] + xorout[79] >= -2
 R3164: - sout[71] - sout[87] + sout[95] - sout[111] + sout[119]
   + sout[127] + xorout[79] >= -2
 R3165: - sout[71] - sout[87] - sout[95] + sout[111] + sout[119]
   + sout[127] + xorout[79] >= -2
 R3166: sout[71] - sout[87] - sout[95] - sout[111] - sout[119] - sout[127]
   + xorout[79] >= -4
 R3167: - sout[71] + sout[87] - sout[95] - sout[111] - sout[119]
   - sout[127] + xorout[79] >= -4
 R3168: - sout[71] - sout[87] + sout[95] - sout[111] - sout[119]
   - sout[127] + xorout[79] >= -4
 R3169: - sout[71] - sout[87] - sout[95] + sout[111] - sout[119]
   - sout[127] + xorout[79] >= -4
 R3170: - sout[71] - sout[87] - sout[95] - sout[111] + sout[119]
   - sout[127] + xorout[79] >= -4
 R3171: - sout[71] - sout[87] - sout[95] - sout[111] - sout[119]
   + sout[127] + xorout[79] >= -4
 R3172: - sout[126] + xorout[78] >= 0
 R3173: - sout[118] + xorout[78] >= 0
 R3174: - sout[110] + xorout[78] >= 0
 R3175: - sout[94] + xorout[78] >= 0
 R3176: - sout[86] + xorout[78] >= 0
 R3177: - sout[70] + xorout[78] >= 0
 R3178: - sout[70] + sout[71] >= 0
 R3179: - sout[86] + sout[87] >= 0
 R3180: - sout[94] + sout[95] >= 0
 R3181: - sout[110] + sout[111] >= 0
 R3182: - sout[118] + sout[119] >= 0
 R3183: - sout[126] + sout[127] >= 0
 R3184: sout[70] + sout[86] + sout[94] + sout[110] + sout[118] + sout[126]
   - xorout[78] >= 0
 R3185: - xorout[78] + xorout[79] >= 0
 R3186: - sout[71] - sout[87] - sout[95] - sout[111] - sout[119]
   - sout[127] + xorout[78] - xorout[79] >= -6
 R3187: sout[71] + sout[87] - sout[95] - sout[111] - sout[119] - sout[127]
   + xorout[78] - xorout[79] >= -4
 R3188: sout[71] - sout[87] + sout[95] - sout[111] - sout[119] - sout[127]
   + xorout[78] - xorout[79] >= -4
 R3189: - sout[71] + sout[87] + sout[95] - sout[111] - sout[119]
   - sout[127] + xorout[78] - xorout[79] >= -4
 R3190: sout[71] - sout[87] - sout[95] + sout[111] - sout[119] - sout[127]
   + xorout[78] - xorout[79] >= -4
 R3191: - sout[71] + sout[87] - sout[95] + sout[111] - sout[119]
   - sout[127] + xorout[78] - xorout[79] >= -4
 R3192: - sout[71] - sout[87] + sout[95] + sout[111] - sout[119]
   - sout[127] + xorout[78] - xorout[79] >= -4
 R3193: sout[71] - sout[87] - sout[95] - sout[111] + sout[119] - sout[127]
   + xorout[78] - xorout[79] >= -4
 R3194: - sout[71] + sout[87] - sout[95] - sout[111] + sout[119]
   - sout[127] + xorout[78] - xorout[79] >= -4
 R3195: - sout[71] - sout[87] + sout[95] - sout[111] + sout[119]
   - sout[127] + xorout[78] - xorout[79] >= -4
 R3196: - sout[71] - sout[87] - sout[95] + sout[111] + sout[119]
   - sout[127] + xorout[78] - xorout[79] >= -4
 R3197: sout[71] - sout[87] - sout[95] - sout[111] - sout[119] + sout[127]
   + xorout[78] - xorout[79] >= -4
 R3198: - sout[71] + sout[87] - sout[95] - sout[111] - sout[119]
   + sout[127] + xorout[78] - xorout[79] >= -4
 R3199: - sout[71] - sout[87] + sout[95] - sout[111] - sout[119]
   + sout[127] + xorout[78] - xorout[79] >= -4
 R3200: - sout[71] - sout[87] - sout[95] + sout[111] - sout[119]
   + sout[127] + xorout[78] - xorout[79] >= -4
 R3201: - sout[71] - sout[87] - sout[95] - sout[111] + sout[119]
   + sout[127] + xorout[78] - xorout[79] >= -4
 R3202: sout[71] + sout[87] + sout[95] + sout[111] - sout[119] - sout[127]
   + xorout[78] - xorout[79] >= -2
 R3203: sout[71] + sout[87] + sout[95] - sout[111] + sout[119] - sout[127]
   + xorout[78] - xorout[79] >= -2
 R3204: sout[71] + sout[87] - sout[95] + sout[111] + sout[119] - sout[127]
   + xorout[78] - xorout[79] >= -2
 R3205: sout[71] - sout[87] + sout[95] + sout[111] + sout[119] - sout[127]
   + xorout[78] - xorout[79] >= -2
 R3206: - sout[71] + sout[87] + sout[95] + sout[111] + sout[119]
   - sout[127] + xorout[78] - xorout[79] >= -2
 R3207: sout[71] + sout[87] + sout[95] - sout[111] - sout[119] + sout[127]
   + xorout[78] - xorout[79] >= -2
 R3208: sout[71] + sout[87] - sout[95] + sout[111] - sout[119] + sout[127]
   + xorout[78] - xorout[79] >= -2
 R3209: sout[71] - sout[87] + sout[95] + sout[111] - sout[119] + sout[127]
   + xorout[78] - xorout[79] >= -2
 R3210: - sout[71] + sout[87] + sout[95] + sout[111] - sout[119]
   + sout[127] + xorout[78] - xorout[79] >= -2
 R3211: sout[71] + sout[87] - sout[95] - sout[111] + sout[119] + sout[127]
   + xorout[78] - xorout[79] >= -2
 R3212: sout[71] - sout[87] + sout[95] - sout[111] + sout[119] + sout[127]
   + xorout[78] - xorout[79] >= -2
 R3213: - sout[71] + sout[87] + sout[95] - sout[111] + sout[119]
   + sout[127] + xorout[78] - xorout[79] >= -2
 R3214: sout[71] - sout[87] - sout[95] + sout[111] + sout[119] + sout[127]
   + xorout[78] - xorout[79] >= -2
 R3215: - sout[71] + sout[87] - sout[95] + sout[111] + sout[119]
   + sout[127] + xorout[78] - xorout[79] >= -2
 R3216: - sout[71] - sout[87] + sout[95] + sout[111] + sout[119]
   + sout[127] + xorout[78] - xorout[79] >= -2
 R3217: sout[65] + sout[73] + sout[89] + sout[97] + sout[113] - sout[121]
   + xorout[81] >= 0
 R3218: sout[65] + sout[73] + sout[89] + sout[97] - sout[113] + sout[121]
   + xorout[81] >= 0
 R3219: sout[65] + sout[73] + sout[89] - sout[97] + sout[113] + sout[121]
   + xorout[81] >= 0
 R3220: sout[65] + sout[73] - sout[89] + sout[97] + sout[113] + sout[121]
   + xorout[81] >= 0
 R3221: sout[65] - sout[73] + sout[89] + sout[97] + sout[113] + sout[121]
   + xorout[81] >= 0
 R3222: - sout[65] + sout[73] + sout[89] + sout[97] + sout[113] + sout[121]
   + xorout[81] >= 0
 R3223: sout[65] + sout[73] + sout[89] + sout[97] + sout[113] + sout[121]
   - xorout[81] >= 0
 R3224: sout[65] + sout[73] + sout[89] - sout[97] - sout[113] - sout[121]
   + xorout[81] >= -2
 R3225: sout[65] + sout[73] - sout[89] + sout[97] - sout[113] - sout[121]
   + xorout[81] >= -2
 R3226: sout[65] - sout[73] + sout[89] + sout[97] - sout[113] - sout[121]
   + xorout[81] >= -2
 R3227: - sout[65] + sout[73] + sout[89] + sout[97] - sout[113] - sout[121]
   + xorout[81] >= -2
 R3228: sout[65] + sout[73] - sout[89] - sout[97] + sout[113] - sout[121]
   + xorout[81] >= -2
 R3229: sout[65] - sout[73] + sout[89] - sout[97] + sout[113] - sout[121]
   + xorout[81] >= -2
 R3230: - sout[65] + sout[73] + sout[89] - sout[97] + sout[113] - sout[121]
   + xorout[81] >= -2
 R3231: sout[65] - sout[73] - sout[89] + sout[97] + sout[113] - sout[121]
   + xorout[81] >= -2
 R3232: - sout[65] + sout[73] - sout[89] + sout[97] + sout[113] - sout[121]
   + xorout[81] >= -2
 R3233: - sout[65] - sout[73] + sout[89] + sout[97] + sout[113] - sout[121]
   + xorout[81] >= -2
 R3234: sout[65] + sout[73] - sout[89] - sout[97] - sout[113] + sout[121]
   + xorout[81] >= -2
 R3235: sout[65] - sout[73] + sout[89] - sout[97] - sout[113] + sout[121]
   + xorout[81] >= -2
 R3236: - sout[65] + sout[73] + sout[89] - sout[97] - sout[113] + sout[121]
   + xorout[81] >= -2
 R3237: sout[65] - sout[73] - sout[89] + sout[97] - sout[113] + sout[121]
   + xorout[81] >= -2
 R3238: - sout[65] + sout[73] - sout[89] + sout[97] - sout[113] + sout[121]
   + xorout[81] >= -2
 R3239: - sout[65] - sout[73] + sout[89] + sout[97] - sout[113] + sout[121]
   + xorout[81] >= -2
 R3240: sout[65] - sout[73] - sout[89] - sout[97] + sout[113] + sout[121]
   + xorout[81] >= -2
 R3241: - sout[65] + sout[73] - sout[89] - sout[97] + sout[113] + sout[121]
   + xorout[81] >= -2
 R3242: - sout[65] - sout[73] + sout[89] - sout[97] + sout[113] + sout[121]
   + xorout[81] >= -2
 R3243: - sout[65] - sout[73] - sout[89] + sout[97] + sout[113] + sout[121]
   + xorout[81] >= -2
 R3244: sout[65] - sout[73] - sout[89] - sout[97] - sout[113] - sout[121]
   + xorout[81] >= -4
 R3245: - sout[65] + sout[73] - sout[89] - sout[97] - sout[113] - sout[121]
   + xorout[81] >= -4
 R3246: - sout[65] - sout[73] + sout[89] - sout[97] - sout[113] - sout[121]
   + xorout[81] >= -4
 R3247: - sout[65] - sout[73] - sout[89] + sout[97] - sout[113] - sout[121]
   + xorout[81] >= -4
 R3248: - sout[65] - sout[73] - sout[89] - sout[97] + sout[113] - sout[121]
   + xorout[81] >= -4
 R3249: - sout[65] - sout[73] - sout[89] - sout[97] - sout[113] + sout[121]
   + xorout[81] >= -4
 R3250: - sout[120] + xorout[80] >= 0
 R3251: - sout[112] + xorout[80] >= 0
 R3252: - sout[96] + xorout[80] >= 0
 R3253: - sout[88] + xorout[80] >= 0
 R3254: - sout[72] + xorout[80] >= 0
 R3255: - sout[64] + xorout[80] >= 0
 R3256: - sout[64] + sout[65] >= 0
 R3257: - sout[72] + sout[73] >= 0
 R3258: - sout[88] + sout[89] >= 0
 R3259: - sout[96] + sout[97] >= 0
 R3260: - sout[112] + sout[113] >= 0
 R3261: - sout[120] + sout[121] >= 0
 R3262: sout[64] + sout[72] + sout[88] + sout[96] + sout[112] + sout[120]
   - xorout[80] >= 0
 R3263: - xorout[80] + xorout[81] >= 0
 R3264: - sout[65] - sout[73] - sout[89] - sout[97] - sout[113] - sout[121]
   + xorout[80] - xorout[81] >= -6
 R3265: sout[65] + sout[73] - sout[89] - sout[97] - sout[113] - sout[121]
   + xorout[80] - xorout[81] >= -4
 R3266: sout[65] - sout[73] + sout[89] - sout[97] - sout[113] - sout[121]
   + xorout[80] - xorout[81] >= -4
 R3267: - sout[65] + sout[73] + sout[89] - sout[97] - sout[113] - sout[121]
   + xorout[80] - xorout[81] >= -4
 R3268: sout[65] - sout[73] - sout[89] + sout[97] - sout[113] - sout[121]
   + xorout[80] - xorout[81] >= -4
 R3269: - sout[65] + sout[73] - sout[89] + sout[97] - sout[113] - sout[121]
   + xorout[80] - xorout[81] >= -4
 R3270: - sout[65] - sout[73] + sout[89] + sout[97] - sout[113] - sout[121]
   + xorout[80] - xorout[81] >= -4
 R3271: sout[65] - sout[73] - sout[89] - sout[97] + sout[113] - sout[121]
   + xorout[80] - xorout[81] >= -4
 R3272: - sout[65] + sout[73] - sout[89] - sout[97] + sout[113] - sout[121]
   + xorout[80] - xorout[81] >= -4
 R3273: - sout[65] - sout[73] + sout[89] - sout[97] + sout[113] - sout[121]
   + xorout[80] - xorout[81] >= -4
 R3274: - sout[65] - sout[73] - sout[89] + sout[97] + sout[113] - sout[121]
   + xorout[80] - xorout[81] >= -4
 R3275: sout[65] - sout[73] - sout[89] - sout[97] - sout[113] + sout[121]
   + xorout[80] - xorout[81] >= -4
 R3276: - sout[65] + sout[73] - sout[89] - sout[97] - sout[113] + sout[121]
   + xorout[80] - xorout[81] >= -4
 R3277: - sout[65] - sout[73] + sout[89] - sout[97] - sout[113] + sout[121]
   + xorout[80] - xorout[81] >= -4
 R3278: - sout[65] - sout[73] - sout[89] + sout[97] - sout[113] + sout[121]
   + xorout[80] - xorout[81] >= -4
 R3279: - sout[65] - sout[73] - sout[89] - sout[97] + sout[113] + sout[121]
   + xorout[80] - xorout[81] >= -4
 R3280: sout[65] + sout[73] + sout[89] + sout[97] - sout[113] - sout[121]
   + xorout[80] - xorout[81] >= -2
 R3281: sout[65] + sout[73] + sout[89] - sout[97] + sout[113] - sout[121]
   + xorout[80] - xorout[81] >= -2
 R3282: sout[65] + sout[73] - sout[89] + sout[97] + sout[113] - sout[121]
   + xorout[80] - xorout[81] >= -2
 R3283: sout[65] - sout[73] + sout[89] + sout[97] + sout[113] - sout[121]
   + xorout[80] - xorout[81] >= -2
 R3284: - sout[65] + sout[73] + sout[89] + sout[97] + sout[113] - sout[121]
   + xorout[80] - xorout[81] >= -2
 R3285: sout[65] + sout[73] + sout[89] - sout[97] - sout[113] + sout[121]
   + xorout[80] - xorout[81] >= -2
 R3286: sout[65] + sout[73] - sout[89] + sout[97] - sout[113] + sout[121]
   + xorout[80] - xorout[81] >= -2
 R3287: sout[65] - sout[73] + sout[89] + sout[97] - sout[113] + sout[121]
   + xorout[80] - xorout[81] >= -2
 R3288: - sout[65] + sout[73] + sout[89] + sout[97] - sout[113] + sout[121]
   + xorout[80] - xorout[81] >= -2
 R3289: sout[65] + sout[73] - sout[89] - sout[97] + sout[113] + sout[121]
   + xorout[80] - xorout[81] >= -2
 R3290: sout[65] - sout[73] + sout[89] - sout[97] + sout[113] + sout[121]
   + xorout[80] - xorout[81] >= -2
 R3291: - sout[65] + sout[73] + sout[89] - sout[97] + sout[113] + sout[121]
   + xorout[80] - xorout[81] >= -2
 R3292: sout[65] - sout[73] - sout[89] + sout[97] + sout[113] + sout[121]
   + xorout[80] - xorout[81] >= -2
 R3293: - sout[65] + sout[73] - sout[89] + sout[97] + sout[113] + sout[121]
   + xorout[80] - xorout[81] >= -2
 R3294: - sout[65] - sout[73] + sout[89] + sout[97] + sout[113] + sout[121]
   + xorout[80] - xorout[81] >= -2
 R3295: sout[67] + sout[75] + sout[91] + sout[99] + sout[115] - sout[123]
   + xorout[83] >= 0
 R3296: sout[67] + sout[75] + sout[91] + sout[99] - sout[115] + sout[123]
   + xorout[83] >= 0
 R3297: sout[67] + sout[75] + sout[91] - sout[99] + sout[115] + sout[123]
   + xorout[83] >= 0
 R3298: sout[67] + sout[75] - sout[91] + sout[99] + sout[115] + sout[123]
   + xorout[83] >= 0
 R3299: sout[67] - sout[75] + sout[91] + sout[99] + sout[115] + sout[123]
   + xorout[83] >= 0
 R3300: - sout[67] + sout[75] + sout[91] + sout[99] + sout[115] + sout[123]
   + xorout[83] >= 0
 R3301: sout[67] + sout[75] + sout[91] + sout[99] + sout[115] + sout[123]
   - xorout[83] >= 0
 R3302: sout[67] + sout[75] + sout[91] - sout[99] - sout[115] - sout[123]
   + xorout[83] >= -2
 R3303: sout[67] + sout[75] - sout[91] + sout[99] - sout[115] - sout[123]
   + xorout[83] >= -2
 R3304: sout[67] - sout[75] + sout[91] + sout[99] - sout[115] - sout[123]
   + xorout[83] >= -2
 R3305: - sout[67] + sout[75] + sout[91] + sout[99] - sout[115] - sout[123]
   + xorout[83] >= -2
 R3306: sout[67] + sout[75] - sout[91] - sout[99] + sout[115] - sout[123]
   + xorout[83] >= -2
 R3307: sout[67] - sout[75] + sout[91] - sout[99] + sout[115] - sout[123]
   + xorout[83] >= -2
 R3308: - sout[67] + sout[75] + sout[91] - sout[99] + sout[115] - sout[123]
   + xorout[83] >= -2
 R3309: sout[67] - sout[75] - sout[91] + sout[99] + sout[115] - sout[123]
   + xorout[83] >= -2
 R3310: - sout[67] + sout[75] - sout[91] + sout[99] + sout[115] - sout[123]
   + xorout[83] >= -2
 R3311: - sout[67] - sout[75] + sout[91] + sout[99] + sout[115] - sout[123]
   + xorout[83] >= -2
 R3312: sout[67] + sout[75] - sout[91] - sout[99] - sout[115] + sout[123]
   + xorout[83] >= -2
 R3313: sout[67] - sout[75] + sout[91] - sout[99] - sout[115] + sout[123]
   + xorout[83] >= -2
 R3314: - sout[67] + sout[75] + sout[91] - sout[99] - sout[115] + sout[123]
   + xorout[83] >= -2
 R3315: sout[67] - sout[75] - sout[91] + sout[99] - sout[115] + sout[123]
   + xorout[83] >= -2
 R3316: - sout[67] + sout[75] - sout[91] + sout[99] - sout[115] + sout[123]
   + xorout[83] >= -2
 R3317: - sout[67] - sout[75] + sout[91] + sout[99] - sout[115] + sout[123]
   + xorout[83] >= -2
 R3318: sout[67] - sout[75] - sout[91] - sout[99] + sout[115] + sout[123]
   + xorout[83] >= -2
 R3319: - sout[67] + sout[75] - sout[91] - sout[99] + sout[115] + sout[123]
   + xorout[83] >= -2
 R3320: - sout[67] - sout[75] + sout[91] - sout[99] + sout[115] + sout[123]
   + xorout[83] >= -2
 R3321: - sout[67] - sout[75] - sout[91] + sout[99] + sout[115] + sout[123]
   + xorout[83] >= -2
 R3322: sout[67] - sout[75] - sout[91] - sout[99] - sout[115] - sout[123]
   + xorout[83] >= -4
 R3323: - sout[67] + sout[75] - sout[91] - sout[99] - sout[115] - sout[123]
   + xorout[83] >= -4
 R3324: - sout[67] - sout[75] + sout[91] - sout[99] - sout[115] - sout[123]
   + xorout[83] >= -4
 R3325: - sout[67] - sout[75] - sout[91] + sout[99] - sout[115] - sout[123]
   + xorout[83] >= -4
 R3326: - sout[67] - sout[75] - sout[91] - sout[99] + sout[115] - sout[123]
   + xorout[83] >= -4
 R3327: - sout[67] - sout[75] - sout[91] - sout[99] - sout[115] + sout[123]
   + xorout[83] >= -4
 R3328: - sout[122] + xorout[82] >= 0
 R3329: - sout[114] + xorout[82] >= 0
 R3330: - sout[98] + xorout[82] >= 0
 R3331: - sout[90] + xorout[82] >= 0
 R3332: - sout[74] + xorout[82] >= 0
 R3333: - sout[66] + xorout[82] >= 0
 R3334: - sout[66] + sout[67] >= 0
 R3335: - sout[74] + sout[75] >= 0
 R3336: - sout[90] + sout[91] >= 0
 R3337: - sout[98] + sout[99] >= 0
 R3338: - sout[114] + sout[115] >= 0
 R3339: - sout[122] + sout[123] >= 0
 R3340: sout[66] + sout[74] + sout[90] + sout[98] + sout[114] + sout[122]
   - xorout[82] >= 0
 R3341: - xorout[82] + xorout[83] >= 0
 R3342: - sout[67] - sout[75] - sout[91] - sout[99] - sout[115] - sout[123]
   + xorout[82] - xorout[83] >= -6
 R3343: sout[67] + sout[75] - sout[91] - sout[99] - sout[115] - sout[123]
   + xorout[82] - xorout[83] >= -4
 R3344: sout[67] - sout[75] + sout[91] - sout[99] - sout[115] - sout[123]
   + xorout[82] - xorout[83] >= -4
 R3345: - sout[67] + sout[75] + sout[91] - sout[99] - sout[115] - sout[123]
   + xorout[82] - xorout[83] >= -4
 R3346: sout[67] - sout[75] - sout[91] + sout[99] - sout[115] - sout[123]
   + xorout[82] - xorout[83] >= -4
 R3347: - sout[67] + sout[75] - sout[91] + sout[99] - sout[115] - sout[123]
   + xorout[82] - xorout[83] >= -4
 R3348: - sout[67] - sout[75] + sout[91] + sout[99] - sout[115] - sout[123]
   + xorout[82] - xorout[83] >= -4
 R3349: sout[67] - sout[75] - sout[91] - sout[99] + sout[115] - sout[123]
   + xorout[82] - xorout[83] >= -4
 R3350: - sout[67] + sout[75] - sout[91] - sout[99] + sout[115] - sout[123]
   + xorout[82] - xorout[83] >= -4
 R3351: - sout[67] - sout[75] + sout[91] - sout[99] + sout[115] - sout[123]
   + xorout[82] - xorout[83] >= -4
 R3352: - sout[67] - sout[75] - sout[91] + sout[99] + sout[115] - sout[123]
   + xorout[82] - xorout[83] >= -4
 R3353: sout[67] - sout[75] - sout[91] - sout[99] - sout[115] + sout[123]
   + xorout[82] - xorout[83] >= -4
 R3354: - sout[67] + sout[75] - sout[91] - sout[99] - sout[115] + sout[123]
   + xorout[82] - xorout[83] >= -4
 R3355: - sout[67] - sout[75] + sout[91] - sout[99] - sout[115] + sout[123]
   + xorout[82] - xorout[83] >= -4
 R3356: - sout[67] - sout[75] - sout[91] + sout[99] - sout[115] + sout[123]
   + xorout[82] - xorout[83] >= -4
 R3357: - sout[67] - sout[75] - sout[91] - sout[99] + sout[115] + sout[123]
   + xorout[82] - xorout[83] >= -4
 R3358: sout[67] + sout[75] + sout[91] + sout[99] - sout[115] - sout[123]
   + xorout[82] - xorout[83] >= -2
 R3359: sout[67] + sout[75] + sout[91] - sout[99] + sout[115] - sout[123]
   + xorout[82] - xorout[83] >= -2
 R3360: sout[67] + sout[75] - sout[91] + sout[99] + sout[115] - sout[123]
   + xorout[82] - xorout[83] >= -2
 R3361: sout[67] - sout[75] + sout[91] + sout[99] + sout[115] - sout[123]
   + xorout[82] - xorout[83] >= -2
 R3362: - sout[67] + sout[75] + sout[91] + sout[99] + sout[115] - sout[123]
   + xorout[82] - xorout[83] >= -2
 R3363: sout[67] + sout[75] + sout[91] - sout[99] - sout[115] + sout[123]
   + xorout[82] - xorout[83] >= -2
 R3364: sout[67] + sout[75] - sout[91] + sout[99] - sout[115] + sout[123]
   + xorout[82] - xorout[83] >= -2
 R3365: sout[67] - sout[75] + sout[91] + sout[99] - sout[115] + sout[123]
   + xorout[82] - xorout[83] >= -2
 R3366: - sout[67] + sout[75] + sout[91] + sout[99] - sout[115] + sout[123]
   + xorout[82] - xorout[83] >= -2
 R3367: sout[67] + sout[75] - sout[91] - sout[99] + sout[115] + sout[123]
   + xorout[82] - xorout[83] >= -2
 R3368: sout[67] - sout[75] + sout[91] - sout[99] + sout[115] + sout[123]
   + xorout[82] - xorout[83] >= -2
 R3369: - sout[67] + sout[75] + sout[91] - sout[99] + sout[115] + sout[123]
   + xorout[82] - xorout[83] >= -2
 R3370: sout[67] - sout[75] - sout[91] + sout[99] + sout[115] + sout[123]
   + xorout[82] - xorout[83] >= -2
 R3371: - sout[67] + sout[75] - sout[91] + sout[99] + sout[115] + sout[123]
   + xorout[82] - xorout[83] >= -2
 R3372: - sout[67] - sout[75] + sout[91] + sout[99] + sout[115] + sout[123]
   + xorout[82] - xorout[83] >= -2
 R3373: sout[69] + sout[77] + sout[93] + sout[101] + sout[117] - sout[125]
   + xorout[85] >= 0
 R3374: sout[69] + sout[77] + sout[93] + sout[101] - sout[117] + sout[125]
   + xorout[85] >= 0
 R3375: sout[69] + sout[77] + sout[93] - sout[101] + sout[117] + sout[125]
   + xorout[85] >= 0
 R3376: sout[69] + sout[77] - sout[93] + sout[101] + sout[117] + sout[125]
   + xorout[85] >= 0
 R3377: sout[69] - sout[77] + sout[93] + sout[101] + sout[117] + sout[125]
   + xorout[85] >= 0
 R3378: - sout[69] + sout[77] + sout[93] + sout[101] + sout[117]
   + sout[125] + xorout[85] >= 0
 R3379: sout[69] + sout[77] + sout[93] + sout[101] + sout[117] + sout[125]
   - xorout[85] >= 0
 R3380: sout[69] + sout[77] + sout[93] - sout[101] - sout[117] - sout[125]
   + xorout[85] >= -2
 R3381: sout[69] + sout[77] - sout[93] + sout[101] - sout[117] - sout[125]
   + xorout[85] >= -2
 R3382: sout[69] - sout[77] + sout[93] + sout[101] - sout[117] - sout[125]
   + xorout[85] >= -2
 R3383: - sout[69] + sout[77] + sout[93] + sout[101] - sout[117]
   - sout[125] + xorout[85] >= -2
 R3384: sout[69] + sout[77] - sout[93] - sout[101] + sout[117] - sout[125]
   + xorout[85] >= -2
 R3385: sout[69] - sout[77] + sout[93] - sout[101] + sout[117] - sout[125]
   + xorout[85] >= -2
 R3386: - sout[69] + sout[77] + sout[93] - sout[101] + sout[117]
   - sout[125] + xorout[85] >= -2
 R3387: sout[69] - sout[77] - sout[93] + sout[101] + sout[117] - sout[125]
   + xorout[85] >= -2
 R3388: - sout[69] + sout[77] - sout[93] + sout[101] + sout[117]
   - sout[125] + xorout[85] >= -2
 R3389: - sout[69] - sout[77] + sout[93] + sout[101] + sout[117]
   - sout[125] + xorout[85] >= -2
 R3390: sout[69] + sout[77] - sout[93] - sout[101] - sout[117] + sout[125]
   + xorout[85] >= -2
 R3391: sout[69] - sout[77] + sout[93] - sout[101] - sout[117] + sout[125]
   + xorout[85] >= -2
 R3392: - sout[69] + sout[77] + sout[93] - sout[101] - sout[117]
   + sout[125] + xorout[85] >= -2
 R3393: sout[69] - sout[77] - sout[93] + sout[101] - sout[117] + sout[125]
   + xorout[85] >= -2
 R3394: - sout[69] + sout[77] - sout[93] + sout[101] - sout[117]
   + sout[125] + xorout[85] >= -2
 R3395: - sout[69] - sout[77] + sout[93] + sout[101] - sout[117]
   + sout[125] + xorout[85] >= -2
 R3396: sout[69] - sout[77] - sout[93] - sout[101] + sout[117] + sout[125]
   + xorout[85] >= -2
 R3397: - sout[69] + sout[77] - sout[93] - sout[101] + sout[117]
   + sout[125] + xorout[85] >= -2
 R3398: - sout[69] - sout[77] + sout[93] - sout[101] + sout[117]
   + sout[125] + xorout[85] >= -2
 R3399: - sout[69] - sout[77] - sout[93] + sout[101] + sout[117]
   + sout[125] + xorout[85] >= -2
 R3400: sout[69] - sout[77] - sout[93] - sout[101] - sout[117] - sout[125]
   + xorout[85] >= -4
 R3401: - sout[69] + sout[77] - sout[93] - sout[101] - sout[117]
   - sout[125] + xorout[85] >= -4
 R3402: - sout[69] - sout[77] + sout[93] - sout[101] - sout[117]
   - sout[125] + xorout[85] >= -4
 R3403: - sout[69] - sout[77] - sout[93] + sout[101] - sout[117]
   - sout[125] + xorout[85] >= -4
 R3404: - sout[69] - sout[77] - sout[93] - sout[101] + sout[117]
   - sout[125] + xorout[85] >= -4
 R3405: - sout[69] - sout[77] - sout[93] - sout[101] - sout[117]
   + sout[125] + xorout[85] >= -4
 R3406: - sout[124] + xorout[84] >= 0
 R3407: - sout[116] + xorout[84] >= 0
 R3408: - sout[100] + xorout[84] >= 0
 R3409: - sout[92] + xorout[84] >= 0
 R3410: - sout[76] + xorout[84] >= 0
 R3411: - sout[68] + xorout[84] >= 0
 R3412: - sout[68] + sout[69] >= 0
 R3413: - sout[76] + sout[77] >= 0
 R3414: - sout[92] + sout[93] >= 0
 R3415: - sout[100] + sout[101] >= 0
 R3416: - sout[116] + sout[117] >= 0
 R3417: - sout[124] + sout[125] >= 0
 R3418: sout[68] + sout[76] + sout[92] + sout[100] + sout[116] + sout[124]
   - xorout[84] >= 0
 R3419: - xorout[84] + xorout[85] >= 0
 R3420: - sout[69] - sout[77] - sout[93] - sout[101] - sout[117]
   - sout[125] + xorout[84] - xorout[85] >= -6
 R3421: sout[69] + sout[77] - sout[93] - sout[101] - sout[117] - sout[125]
   + xorout[84] - xorout[85] >= -4
 R3422: sout[69] - sout[77] + sout[93] - sout[101] - sout[117] - sout[125]
   + xorout[84] - xorout[85] >= -4
 R3423: - sout[69] + sout[77] + sout[93] - sout[101] - sout[117]
   - sout[125] + xorout[84] - xorout[85] >= -4
 R3424: sout[69] - sout[77] - sout[93] + sout[101] - sout[117] - sout[125]
   + xorout[84] - xorout[85] >= -4
 R3425: - sout[69] + sout[77] - sout[93] + sout[101] - sout[117]
   - sout[125] + xorout[84] - xorout[85] >= -4
 R3426: - sout[69] - sout[77] + sout[93] + sout[101] - sout[117]
   - sout[125] + xorout[84] - xorout[85] >= -4
 R3427: sout[69] - sout[77] - sout[93] - sout[101] + sout[117] - sout[125]
   + xorout[84] - xorout[85] >= -4
 R3428: - sout[69] + sout[77] - sout[93] - sout[101] + sout[117]
   - sout[125] + xorout[84] - xorout[85] >= -4
 R3429: - sout[69] - sout[77] + sout[93] - sout[101] + sout[117]
   - sout[125] + xorout[84] - xorout[85] >= -4
 R3430: - sout[69] - sout[77] - sout[93] + sout[101] + sout[117]
   - sout[125] + xorout[84] - xorout[85] >= -4
 R3431: sout[69] - sout[77] - sout[93] - sout[101] - sout[117] + sout[125]
   + xorout[84] - xorout[85] >= -4
 R3432: - sout[69] + sout[77] - sout[93] - sout[101] - sout[117]
   + sout[125] + xorout[84] - xorout[85] >= -4
 R3433: - sout[69] - sout[77] + sout[93] - sout[101] - sout[117]
   + sout[125] + xorout[84] - xorout[85] >= -4
 R3434: - sout[69] - sout[77] - sout[93] + sout[101] - sout[117]
   + sout[125] + xorout[84] - xorout[85] >= -4
 R3435: - sout[69] - sout[77] - sout[93] - sout[101] + sout[117]
   + sout[125] + xorout[84] - xorout[85] >= -4
 R3436: sout[69] + sout[77] + sout[93] + sout[101] - sout[117] - sout[125]
   + xorout[84] - xorout[85] >= -2
 R3437: sout[69] + sout[77] + sout[93] - sout[101] + sout[117] - sout[125]
   + xorout[84] - xorout[85] >= -2
 R3438: sout[69] + sout[77] - sout[93] + sout[101] + sout[117] - sout[125]
   + xorout[84] - xorout[85] >= -2
 R3439: sout[69] - sout[77] + sout[93] + sout[101] + sout[117] - sout[125]
   + xorout[84] - xorout[85] >= -2
 R3440: - sout[69] + sout[77] + sout[93] + sout[101] + sout[117]
   - sout[125] + xorout[84] - xorout[85] >= -2
 R3441: sout[69] + sout[77] + sout[93] - sout[101] - sout[117] + sout[125]
   + xorout[84] - xorout[85] >= -2
 R3442: sout[69] + sout[77] - sout[93] + sout[101] - sout[117] + sout[125]
   + xorout[84] - xorout[85] >= -2
 R3443: sout[69] - sout[77] + sout[93] + sout[101] - sout[117] + sout[125]
   + xorout[84] - xorout[85] >= -2
 R3444: - sout[69] + sout[77] + sout[93] + sout[101] - sout[117]
   + sout[125] + xorout[84] - xorout[85] >= -2
 R3445: sout[69] + sout[77] - sout[93] - sout[101] + sout[117] + sout[125]
   + xorout[84] - xorout[85] >= -2
 R3446: sout[69] - sout[77] + sout[93] - sout[101] + sout[117] + sout[125]
   + xorout[84] - xorout[85] >= -2
 R3447: - sout[69] + sout[77] + sout[93] - sout[101] + sout[117]
   + sout[125] + xorout[84] - xorout[85] >= -2
 R3448: sout[69] - sout[77] - sout[93] + sout[101] + sout[117] + sout[125]
   + xorout[84] - xorout[85] >= -2
 R3449: - sout[69] + sout[77] - sout[93] + sout[101] + sout[117]
   + sout[125] + xorout[84] - xorout[85] >= -2
 R3450: - sout[69] - sout[77] + sout[93] + sout[101] + sout[117]
   + sout[125] + xorout[84] - xorout[85] >= -2
 R3451: sout[71] + sout[79] + sout[95] + sout[103] + sout[119] - sout[127]
   + xorout[87] >= 0
 R3452: sout[71] + sout[79] + sout[95] + sout[103] - sout[119] + sout[127]
   + xorout[87] >= 0
 R3453: sout[71] + sout[79] + sout[95] - sout[103] + sout[119] + sout[127]
   + xorout[87] >= 0
 R3454: sout[71] + sout[79] - sout[95] + sout[103] + sout[119] + sout[127]
   + xorout[87] >= 0
 R3455: sout[71] - sout[79] + sout[95] + sout[103] + sout[119] + sout[127]
   + xorout[87] >= 0
 R3456: - sout[71] + sout[79] + sout[95] + sout[103] + sout[119]
   + sout[127] + xorout[87] >= 0
 R3457: sout[71] + sout[79] + sout[95] + sout[103] + sout[119] + sout[127]
   - xorout[87] >= 0
 R3458: sout[71] + sout[79] + sout[95] - sout[103] - sout[119] - sout[127]
   + xorout[87] >= -2
 R3459: sout[71] + sout[79] - sout[95] + sout[103] - sout[119] - sout[127]
   + xorout[87] >= -2
 R3460: sout[71] - sout[79] + sout[95] + sout[103] - sout[119] - sout[127]
   + xorout[87] >= -2
 R3461: - sout[71] + sout[79] + sout[95] + sout[103] - sout[119]
   - sout[127] + xorout[87] >= -2
 R3462: sout[71] + sout[79] - sout[95] - sout[103] + sout[119] - sout[127]
   + xorout[87] >= -2
 R3463: sout[71] - sout[79] + sout[95] - sout[103] + sout[119] - sout[127]
   + xorout[87] >= -2
 R3464: - sout[71] + sout[79] + sout[95] - sout[103] + sout[119]
   - sout[127] + xorout[87] >= -2
 R3465: sout[71] - sout[79] - sout[95] + sout[103] + sout[119] - sout[127]
   + xorout[87] >= -2
 R3466: - sout[71] + sout[79] - sout[95] + sout[103] + sout[119]
   - sout[127] + xorout[87] >= -2
 R3467: - sout[71] - sout[79] + sout[95] + sout[103] + sout[119]
   - sout[127] + xorout[87] >= -2
 R3468: sout[71] + sout[79] - sout[95] - sout[103] - sout[119] + sout[127]
   + xorout[87] >= -2
 R3469: sout[71] - sout[79] + sout[95] - sout[103] - sout[119] + sout[127]
   + xorout[87] >= -2
 R3470: - sout[71] + sout[79] + sout[95] - sout[103] - sout[119]
   + sout[127] + xorout[87] >= -2
 R3471: sout[71] - sout[79] - sout[95] + sout[103] - sout[119] + sout[127]
   + xorout[87] >= -2
 R3472: - sout[71] + sout[79] - sout[95] + sout[103] - sout[119]
   + sout[127] + xorout[87] >= -2
 R3473: - sout[71] - sout[79] + sout[95] + sout[103] - sout[119]
   + sout[127] + xorout[87] >= -2
 R3474: sout[71] - sout[79] - sout[95] - sout[103] + sout[119] + sout[127]
   + xorout[87] >= -2
 R3475: - sout[71] + sout[79] - sout[95] - sout[103] + sout[119]
   + sout[127] + xorout[87] >= -2
 R3476: - sout[71] - sout[79] + sout[95] - sout[103] + sout[119]
   + sout[127] + xorout[87] >= -2
 R3477: - sout[71] - sout[79] - sout[95] + sout[103] + sout[119]
   + sout[127] + xorout[87] >= -2
 R3478: sout[71] - sout[79] - sout[95] - sout[103] - sout[119] - sout[127]
   + xorout[87] >= -4
 R3479: - sout[71] + sout[79] - sout[95] - sout[103] - sout[119]
   - sout[127] + xorout[87] >= -4
 R3480: - sout[71] - sout[79] + sout[95] - sout[103] - sout[119]
   - sout[127] + xorout[87] >= -4
 R3481: - sout[71] - sout[79] - sout[95] + sout[103] - sout[119]
   - sout[127] + xorout[87] >= -4
 R3482: - sout[71] - sout[79] - sout[95] - sout[103] + sout[119]
   - sout[127] + xorout[87] >= -4
 R3483: - sout[71] - sout[79] - sout[95] - sout[103] - sout[119]
   + sout[127] + xorout[87] >= -4
 R3484: - sout[126] + xorout[86] >= 0
 R3485: - sout[118] + xorout[86] >= 0
 R3486: - sout[102] + xorout[86] >= 0
 R3487: - sout[94] + xorout[86] >= 0
 R3488: - sout[78] + xorout[86] >= 0
 R3489: - sout[70] + xorout[86] >= 0
 R3490: - sout[70] + sout[71] >= 0
 R3491: - sout[78] + sout[79] >= 0
 R3492: - sout[94] + sout[95] >= 0
 R3493: - sout[102] + sout[103] >= 0
 R3494: - sout[118] + sout[119] >= 0
 R3495: - sout[126] + sout[127] >= 0
 R3496: sout[70] + sout[78] + sout[94] + sout[102] + sout[118] + sout[126]
   - xorout[86] >= 0
 R3497: - xorout[86] + xorout[87] >= 0
 R3498: - sout[71] - sout[79] - sout[95] - sout[103] - sout[119]
   - sout[127] + xorout[86] - xorout[87] >= -6
 R3499: sout[71] + sout[79] - sout[95] - sout[103] - sout[119] - sout[127]
   + xorout[86] - xorout[87] >= -4
 R3500: sout[71] - sout[79] + sout[95] - sout[103] - sout[119] - sout[127]
   + xorout[86] - xorout[87] >= -4
 R3501: - sout[71] + sout[79] + sout[95] - sout[103] - sout[119]
   - sout[127] + xorout[86] - xorout[87] >= -4
 R3502: sout[71] - sout[79] - sout[95] + sout[103] - sout[119] - sout[127]
   + xorout[86] - xorout[87] >= -4
 R3503: - sout[71] + sout[79] - sout[95] + sout[103] - sout[119]
   - sout[127] + xorout[86] - xorout[87] >= -4
 R3504: - sout[71] - sout[79] + sout[95] + sout[103] - sout[119]
   - sout[127] + xorout[86] - xorout[87] >= -4
 R3505: sout[71] - sout[79] - sout[95] - sout[103] + sout[119] - sout[127]
   + xorout[86] - xorout[87] >= -4
 R3506: - sout[71] + sout[79] - sout[95] - sout[103] + sout[119]
   - sout[127] + xorout[86] - xorout[87] >= -4
 R3507: - sout[71] - sout[79] + sout[95] - sout[103] + sout[119]
   - sout[127] + xorout[86] - xorout[87] >= -4
 R3508: - sout[71] - sout[79] - sout[95] + sout[103] + sout[119]
   - sout[127] + xorout[86] - xorout[87] >= -4
 R3509: sout[71] - sout[79] - sout[95] - sout[103] - sout[119] + sout[127]
   + xorout[86] - xorout[87] >= -4
 R3510: - sout[71] + sout[79] - sout[95] - sout[103] - sout[119]
   + sout[127] + xorout[86] - xorout[87] >= -4
 R3511: - sout[71] - sout[79] + sout[95] - sout[103] - sout[119]
   + sout[127] + xorout[86] - xorout[87] >= -4
 R3512: - sout[71] - sout[79] - sout[95] + sout[103] - sout[119]
   + sout[127] + xorout[86] - xorout[87] >= -4
 R3513: - sout[71] - sout[79] - sout[95] - sout[103] + sout[119]
   + sout[127] + xorout[86] - xorout[87] >= -4
 R3514: sout[71] + sout[79] + sout[95] + sout[103] - sout[119] - sout[127]
   + xorout[86] - xorout[87] >= -2
 R3515: sout[71] + sout[79] + sout[95] - sout[103] + sout[119] - sout[127]
   + xorout[86] - xorout[87] >= -2
 R3516: sout[71] + sout[79] - sout[95] + sout[103] + sout[119] - sout[127]
   + xorout[86] - xorout[87] >= -2
 R3517: sout[71] - sout[79] + sout[95] + sout[103] + sout[119] - sout[127]
   + xorout[86] - xorout[87] >= -2
 R3518: - sout[71] + sout[79] + sout[95] + sout[103] + sout[119]
   - sout[127] + xorout[86] - xorout[87] >= -2
 R3519: sout[71] + sout[79] + sout[95] - sout[103] - sout[119] + sout[127]
   + xorout[86] - xorout[87] >= -2
 R3520: sout[71] + sout[79] - sout[95] + sout[103] - sout[119] + sout[127]
   + xorout[86] - xorout[87] >= -2
 R3521: sout[71] - sout[79] + sout[95] + sout[103] - sout[119] + sout[127]
   + xorout[86] - xorout[87] >= -2
 R3522: - sout[71] + sout[79] + sout[95] + sout[103] - sout[119]
   + sout[127] + xorout[86] - xorout[87] >= -2
 R3523: sout[71] + sout[79] - sout[95] - sout[103] + sout[119] + sout[127]
   + xorout[86] - xorout[87] >= -2
 R3524: sout[71] - sout[79] + sout[95] - sout[103] + sout[119] + sout[127]
   + xorout[86] - xorout[87] >= -2
 R3525: - sout[71] + sout[79] + sout[95] - sout[103] + sout[119]
   + sout[127] + xorout[86] - xorout[87] >= -2
 R3526: sout[71] - sout[79] - sout[95] + sout[103] + sout[119] + sout[127]
   + xorout[86] - xorout[87] >= -2
 R3527: - sout[71] + sout[79] - sout[95] + sout[103] + sout[119]
   + sout[127] + xorout[86] - xorout[87] >= -2
 R3528: - sout[71] - sout[79] + sout[95] + sout[103] + sout[119]
   + sout[127] + xorout[86] - xorout[87] >= -2
 R3529: sout[65] + sout[73] + sout[81] + sout[97] + sout[105] - sout[121]
   + xorout[89] >= 0
 R3530: sout[65] + sout[73] + sout[81] + sout[97] - sout[105] + sout[121]
   + xorout[89] >= 0
 R3531: sout[65] + sout[73] + sout[81] - sout[97] + sout[105] + sout[121]
   + xorout[89] >= 0
 R3532: sout[65] + sout[73] - sout[81] + sout[97] + sout[105] + sout[121]
   + xorout[89] >= 0
 R3533: sout[65] - sout[73] + sout[81] + sout[97] + sout[105] + sout[121]
   + xorout[89] >= 0
 R3534: - sout[65] + sout[73] + sout[81] + sout[97] + sout[105] + sout[121]
   + xorout[89] >= 0
 R3535: sout[65] + sout[73] + sout[81] + sout[97] + sout[105] + sout[121]
   - xorout[89] >= 0
 R3536: sout[65] + sout[73] + sout[81] - sout[97] - sout[105] - sout[121]
   + xorout[89] >= -2
 R3537: sout[65] + sout[73] - sout[81] + sout[97] - sout[105] - sout[121]
   + xorout[89] >= -2
 R3538: sout[65] - sout[73] + sout[81] + sout[97] - sout[105] - sout[121]
   + xorout[89] >= -2
 R3539: - sout[65] + sout[73] + sout[81] + sout[97] - sout[105] - sout[121]
   + xorout[89] >= -2
 R3540: sout[65] + sout[73] - sout[81] - sout[97] + sout[105] - sout[121]
   + xorout[89] >= -2
 R3541: sout[65] - sout[73] + sout[81] - sout[97] + sout[105] - sout[121]
   + xorout[89] >= -2
 R3542: - sout[65] + sout[73] + sout[81] - sout[97] + sout[105] - sout[121]
   + xorout[89] >= -2
 R3543: sout[65] - sout[73] - sout[81] + sout[97] + sout[105] - sout[121]
   + xorout[89] >= -2
 R3544: - sout[65] + sout[73] - sout[81] + sout[97] + sout[105] - sout[121]
   + xorout[89] >= -2
 R3545: - sout[65] - sout[73] + sout[81] + sout[97] + sout[105] - sout[121]
   + xorout[89] >= -2
 R3546: sout[65] + sout[73] - sout[81] - sout[97] - sout[105] + sout[121]
   + xorout[89] >= -2
 R3547: sout[65] - sout[73] + sout[81] - sout[97] - sout[105] + sout[121]
   + xorout[89] >= -2
 R3548: - sout[65] + sout[73] + sout[81] - sout[97] - sout[105] + sout[121]
   + xorout[89] >= -2
 R3549: sout[65] - sout[73] - sout[81] + sout[97] - sout[105] + sout[121]
   + xorout[89] >= -2
 R3550: - sout[65] + sout[73] - sout[81] + sout[97] - sout[105] + sout[121]
   + xorout[89] >= -2
 R3551: - sout[65] - sout[73] + sout[81] + sout[97] - sout[105] + sout[121]
   + xorout[89] >= -2
 R3552: sout[65] - sout[73] - sout[81] - sout[97] + sout[105] + sout[121]
   + xorout[89] >= -2
 R3553: - sout[65] + sout[73] - sout[81] - sout[97] + sout[105] + sout[121]
   + xorout[89] >= -2
 R3554: - sout[65] - sout[73] + sout[81] - sout[97] + sout[105] + sout[121]
   + xorout[89] >= -2
 R3555: - sout[65] - sout[73] - sout[81] + sout[97] + sout[105] + sout[121]
   + xorout[89] >= -2
 R3556: sout[65] - sout[73] - sout[81] - sout[97] - sout[105] - sout[121]
   + xorout[89] >= -4
 R3557: - sout[65] + sout[73] - sout[81] - sout[97] - sout[105] - sout[121]
   + xorout[89] >= -4
 R3558: - sout[65] - sout[73] + sout[81] - sout[97] - sout[105] - sout[121]
   + xorout[89] >= -4
 R3559: - sout[65] - sout[73] - sout[81] + sout[97] - sout[105] - sout[121]
   + xorout[89] >= -4
 R3560: - sout[65] - sout[73] - sout[81] - sout[97] + sout[105] - sout[121]
   + xorout[89] >= -4
 R3561: - sout[65] - sout[73] - sout[81] - sout[97] - sout[105] + sout[121]
   + xorout[89] >= -4
 R3562: - sout[120] + xorout[88] >= 0
 R3563: - sout[104] + xorout[88] >= 0
 R3564: - sout[96] + xorout[88] >= 0
 R3565: - sout[80] + xorout[88] >= 0
 R3566: - sout[72] + xorout[88] >= 0
 R3567: - sout[64] + xorout[88] >= 0
 R3568: - sout[64] + sout[65] >= 0
 R3569: - sout[72] + sout[73] >= 0
 R3570: - sout[80] + sout[81] >= 0
 R3571: - sout[96] + sout[97] >= 0
 R3572: - sout[104] + sout[105] >= 0
 R3573: - sout[120] + sout[121] >= 0
 R3574: sout[64] + sout[72] + sout[80] + sout[96] + sout[104] + sout[120]
   - xorout[88] >= 0
 R3575: - xorout[88] + xorout[89] >= 0
 R3576: - sout[65] - sout[73] - sout[81] - sout[97] - sout[105] - sout[121]
   + xorout[88] - xorout[89] >= -6
 R3577: sout[65] + sout[73] - sout[81] - sout[97] - sout[105] - sout[121]
   + xorout[88] - xorout[89] >= -4
 R3578: sout[65] - sout[73] + sout[81] - sout[97] - sout[105] - sout[121]
   + xorout[88] - xorout[89] >= -4
 R3579: - sout[65] + sout[73] + sout[81] - sout[97] - sout[105] - sout[121]
   + xorout[88] - xorout[89] >= -4
 R3580: sout[65] - sout[73] - sout[81] + sout[97] - sout[105] - sout[121]
   + xorout[88] - xorout[89] >= -4
 R3581: - sout[65] + sout[73] - sout[81] + sout[97] - sout[105] - sout[121]
   + xorout[88] - xorout[89] >= -4
 R3582: - sout[65] - sout[73] + sout[81] + sout[97] - sout[105] - sout[121]
   + xorout[88] - xorout[89] >= -4
 R3583: sout[65] - sout[73] - sout[81] - sout[97] + sout[105] - sout[121]
   + xorout[88] - xorout[89] >= -4
 R3584: - sout[65] + sout[73] - sout[81] - sout[97] + sout[105] - sout[121]
   + xorout[88] - xorout[89] >= -4
 R3585: - sout[65] - sout[73] + sout[81] - sout[97] + sout[105] - sout[121]
   + xorout[88] - xorout[89] >= -4
 R3586: - sout[65] - sout[73] - sout[81] + sout[97] + sout[105] - sout[121]
   + xorout[88] - xorout[89] >= -4
 R3587: sout[65] - sout[73] - sout[81] - sout[97] - sout[105] + sout[121]
   + xorout[88] - xorout[89] >= -4
 R3588: - sout[65] + sout[73] - sout[81] - sout[97] - sout[105] + sout[121]
   + xorout[88] - xorout[89] >= -4
 R3589: - sout[65] - sout[73] + sout[81] - sout[97] - sout[105] + sout[121]
   + xorout[88] - xorout[89] >= -4
 R3590: - sout[65] - sout[73] - sout[81] + sout[97] - sout[105] + sout[121]
   + xorout[88] - xorout[89] >= -4
 R3591: - sout[65] - sout[73] - sout[81] - sout[97] + sout[105] + sout[121]
   + xorout[88] - xorout[89] >= -4
 R3592: sout[65] + sout[73] + sout[81] + sout[97] - sout[105] - sout[121]
   + xorout[88] - xorout[89] >= -2
 R3593: sout[65] + sout[73] + sout[81] - sout[97] + sout[105] - sout[121]
   + xorout[88] - xorout[89] >= -2
 R3594: sout[65] + sout[73] - sout[81] + sout[97] + sout[105] - sout[121]
   + xorout[88] - xorout[89] >= -2
 R3595: sout[65] - sout[73] + sout[81] + sout[97] + sout[105] - sout[121]
   + xorout[88] - xorout[89] >= -2
 R3596: - sout[65] + sout[73] + sout[81] + sout[97] + sout[105] - sout[121]
   + xorout[88] - xorout[89] >= -2
 R3597: sout[65] + sout[73] + sout[81] - sout[97] - sout[105] + sout[121]
   + xorout[88] - xorout[89] >= -2
 R3598: sout[65] + sout[73] - sout[81] + sout[97] - sout[105] + sout[121]
   + xorout[88] - xorout[89] >= -2
 R3599: sout[65] - sout[73] + sout[81] + sout[97] - sout[105] + sout[121]
   + xorout[88] - xorout[89] >= -2
 R3600: - sout[65] + sout[73] + sout[81] + sout[97] - sout[105] + sout[121]
   + xorout[88] - xorout[89] >= -2
 R3601: sout[65] + sout[73] - sout[81] - sout[97] + sout[105] + sout[121]
   + xorout[88] - xorout[89] >= -2
 R3602: sout[65] - sout[73] + sout[81] - sout[97] + sout[105] + sout[121]
   + xorout[88] - xorout[89] >= -2
 R3603: - sout[65] + sout[73] + sout[81] - sout[97] + sout[105] + sout[121]
   + xorout[88] - xorout[89] >= -2
 R3604: sout[65] - sout[73] - sout[81] + sout[97] + sout[105] + sout[121]
   + xorout[88] - xorout[89] >= -2
 R3605: - sout[65] + sout[73] - sout[81] + sout[97] + sout[105] + sout[121]
   + xorout[88] - xorout[89] >= -2
 R3606: - sout[65] - sout[73] + sout[81] + sout[97] + sout[105] + sout[121]
   + xorout[88] - xorout[89] >= -2
 R3607: sout[67] + sout[75] + sout[83] + sout[99] + sout[107] - sout[123]
   + xorout[91] >= 0
 R3608: sout[67] + sout[75] + sout[83] + sout[99] - sout[107] + sout[123]
   + xorout[91] >= 0
 R3609: sout[67] + sout[75] + sout[83] - sout[99] + sout[107] + sout[123]
   + xorout[91] >= 0
 R3610: sout[67] + sout[75] - sout[83] + sout[99] + sout[107] + sout[123]
   + xorout[91] >= 0
 R3611: sout[67] - sout[75] + sout[83] + sout[99] + sout[107] + sout[123]
   + xorout[91] >= 0
 R3612: - sout[67] + sout[75] + sout[83] + sout[99] + sout[107] + sout[123]
   + xorout[91] >= 0
 R3613: sout[67] + sout[75] + sout[83] + sout[99] + sout[107] + sout[123]
   - xorout[91] >= 0
 R3614: sout[67] + sout[75] + sout[83] - sout[99] - sout[107] - sout[123]
   + xorout[91] >= -2
 R3615: sout[67] + sout[75] - sout[83] + sout[99] - sout[107] - sout[123]
   + xorout[91] >= -2
 R3616: sout[67] - sout[75] + sout[83] + sout[99] - sout[107] - sout[123]
   + xorout[91] >= -2
 R3617: - sout[67] + sout[75] + sout[83] + sout[99] - sout[107] - sout[123]
   + xorout[91] >= -2
 R3618: sout[67] + sout[75] - sout[83] - sout[99] + sout[107] - sout[123]
   + xorout[91] >= -2
 R3619: sout[67] - sout[75] + sout[83] - sout[99] + sout[107] - sout[123]
   + xorout[91] >= -2
 R3620: - sout[67] + sout[75] + sout[83] - sout[99] + sout[107] - sout[123]
   + xorout[91] >= -2
 R3621: sout[67] - sout[75] - sout[83] + sout[99] + sout[107] - sout[123]
   + xorout[91] >= -2
 R3622: - sout[67] + sout[75] - sout[83] + sout[99] + sout[107] - sout[123]
   + xorout[91] >= -2
 R3623: - sout[67] - sout[75] + sout[83] + sout[99] + sout[107] - sout[123]
   + xorout[91] >= -2
 R3624: sout[67] + sout[75] - sout[83] - sout[99] - sout[107] + sout[123]
   + xorout[91] >= -2
 R3625: sout[67] - sout[75] + sout[83] - sout[99] - sout[107] + sout[123]
   + xorout[91] >= -2
 R3626: - sout[67] + sout[75] + sout[83] - sout[99] - sout[107] + sout[123]
   + xorout[91] >= -2
 R3627: sout[67] - sout[75] - sout[83] + sout[99] - sout[107] + sout[123]
   + xorout[91] >= -2
 R3628: - sout[67] + sout[75] - sout[83] + sout[99] - sout[107] + sout[123]
   + xorout[91] >= -2
 R3629: - sout[67] - sout[75] + sout[83] + sout[99] - sout[107] + sout[123]
   + xorout[91] >= -2
 R3630: sout[67] - sout[75] - sout[83] - sout[99] + sout[107] + sout[123]
   + xorout[91] >= -2
 R3631: - sout[67] + sout[75] - sout[83] - sout[99] + sout[107] + sout[123]
   + xorout[91] >= -2
 R3632: - sout[67] - sout[75] + sout[83] - sout[99] + sout[107] + sout[123]
   + xorout[91] >= -2
 R3633: - sout[67] - sout[75] - sout[83] + sout[99] + sout[107] + sout[123]
   + xorout[91] >= -2
 R3634: sout[67] - sout[75] - sout[83] - sout[99] - sout[107] - sout[123]
   + xorout[91] >= -4
 R3635: - sout[67] + sout[75] - sout[83] - sout[99] - sout[107] - sout[123]
   + xorout[91] >= -4
 R3636: - sout[67] - sout[75] + sout[83] - sout[99] - sout[107] - sout[123]
   + xorout[91] >= -4
 R3637: - sout[67] - sout[75] - sout[83] + sout[99] - sout[107] - sout[123]
   + xorout[91] >= -4
 R3638: - sout[67] - sout[75] - sout[83] - sout[99] + sout[107] - sout[123]
   + xorout[91] >= -4
 R3639: - sout[67] - sout[75] - sout[83] - sout[99] - sout[107] + sout[123]
   + xorout[91] >= -4
 R3640: - sout[122] + xorout[90] >= 0
 R3641: - sout[106] + xorout[90] >= 0
 R3642: - sout[98] + xorout[90] >= 0
 R3643: - sout[82] + xorout[90] >= 0
 R3644: - sout[74] + xorout[90] >= 0
 R3645: - sout[66] + xorout[90] >= 0
 R3646: - sout[66] + sout[67] >= 0
 R3647: - sout[74] + sout[75] >= 0
 R3648: - sout[82] + sout[83] >= 0
 R3649: - sout[98] + sout[99] >= 0
 R3650: - sout[106] + sout[107] >= 0
 R3651: - sout[122] + sout[123] >= 0
 R3652: sout[66] + sout[74] + sout[82] + sout[98] + sout[106] + sout[122]
   - xorout[90] >= 0
 R3653: - xorout[90] + xorout[91] >= 0
 R3654: - sout[67] - sout[75] - sout[83] - sout[99] - sout[107] - sout[123]
   + xorout[90] - xorout[91] >= -6
 R3655: sout[67] + sout[75] - sout[83] - sout[99] - sout[107] - sout[123]
   + xorout[90] - xorout[91] >= -4
 R3656: sout[67] - sout[75] + sout[83] - sout[99] - sout[107] - sout[123]
   + xorout[90] - xorout[91] >= -4
 R3657: - sout[67] + sout[75] + sout[83] - sout[99] - sout[107] - sout[123]
   + xorout[90] - xorout[91] >= -4
 R3658: sout[67] - sout[75] - sout[83] + sout[99] - sout[107] - sout[123]
   + xorout[90] - xorout[91] >= -4
 R3659: - sout[67] + sout[75] - sout[83] + sout[99] - sout[107] - sout[123]
   + xorout[90] - xorout[91] >= -4
 R3660: - sout[67] - sout[75] + sout[83] + sout[99] - sout[107] - sout[123]
   + xorout[90] - xorout[91] >= -4
 R3661: sout[67] - sout[75] - sout[83] - sout[99] + sout[107] - sout[123]
   + xorout[90] - xorout[91] >= -4
 R3662: - sout[67] + sout[75] - sout[83] - sout[99] + sout[107] - sout[123]
   + xorout[90] - xorout[91] >= -4
 R3663: - sout[67] - sout[75] + sout[83] - sout[99] + sout[107] - sout[123]
   + xorout[90] - xorout[91] >= -4
 R3664: - sout[67] - sout[75] - sout[83] + sout[99] + sout[107] - sout[123]
   + xorout[90] - xorout[91] >= -4
 R3665: sout[67] - sout[75] - sout[83] - sout[99] - sout[107] + sout[123]
   + xorout[90] - xorout[91] >= -4
 R3666: - sout[67] + sout[75] - sout[83] - sout[99] - sout[107] + sout[123]
   + xorout[90] - xorout[91] >= -4
 R3667: - sout[67] - sout[75] + sout[83] - sout[99] - sout[107] + sout[123]
   + xorout[90] - xorout[91] >= -4
 R3668: - sout[67] - sout[75] - sout[83] + sout[99] - sout[107] + sout[123]
   + xorout[90] - xorout[91] >= -4
 R3669: - sout[67] - sout[75] - sout[83] - sout[99] + sout[107] + sout[123]
   + xorout[90] - xorout[91] >= -4
 R3670: sout[67] + sout[75] + sout[83] + sout[99] - sout[107] - sout[123]
   + xorout[90] - xorout[91] >= -2
 R3671: sout[67] + sout[75] + sout[83] - sout[99] + sout[107] - sout[123]
   + xorout[90] - xorout[91] >= -2
 R3672: sout[67] + sout[75] - sout[83] + sout[99] + sout[107] - sout[123]
   + xorout[90] - xorout[91] >= -2
 R3673: sout[67] - sout[75] + sout[83] + sout[99] + sout[107] - sout[123]
   + xorout[90] - xorout[91] >= -2
 R3674: - sout[67] + sout[75] + sout[83] + sout[99] + sout[107] - sout[123]
   + xorout[90] - xorout[91] >= -2
 R3675: sout[67] + sout[75] + sout[83] - sout[99] - sout[107] + sout[123]
   + xorout[90] - xorout[91] >= -2
 R3676: sout[67] + sout[75] - sout[83] + sout[99] - sout[107] + sout[123]
   + xorout[90] - xorout[91] >= -2
 R3677: sout[67] - sout[75] + sout[83] + sout[99] - sout[107] + sout[123]
   + xorout[90] - xorout[91] >= -2
 R3678: - sout[67] + sout[75] + sout[83] + sout[99] - sout[107] + sout[123]
   + xorout[90] - xorout[91] >= -2
 R3679: sout[67] + sout[75] - sout[83] - sout[99] + sout[107] + sout[123]
   + xorout[90] - xorout[91] >= -2
 R3680: sout[67] - sout[75] + sout[83] - sout[99] + sout[107] + sout[123]
   + xorout[90] - xorout[91] >= -2
 R3681: - sout[67] + sout[75] + sout[83] - sout[99] + sout[107] + sout[123]
   + xorout[90] - xorout[91] >= -2
 R3682: sout[67] - sout[75] - sout[83] + sout[99] + sout[107] + sout[123]
   + xorout[90] - xorout[91] >= -2
 R3683: - sout[67] + sout[75] - sout[83] + sout[99] + sout[107] + sout[123]
   + xorout[90] - xorout[91] >= -2
 R3684: - sout[67] - sout[75] + sout[83] + sout[99] + sout[107] + sout[123]
   + xorout[90] - xorout[91] >= -2
 R3685: sout[69] + sout[77] + sout[85] + sout[101] + sout[109] - sout[125]
   + xorout[93] >= 0
 R3686: sout[69] + sout[77] + sout[85] + sout[101] - sout[109] + sout[125]
   + xorout[93] >= 0
 R3687: sout[69] + sout[77] + sout[85] - sout[101] + sout[109] + sout[125]
   + xorout[93] >= 0
 R3688: sout[69] + sout[77] - sout[85] + sout[101] + sout[109] + sout[125]
   + xorout[93] >= 0
 R3689: sout[69] - sout[77] + sout[85] + sout[101] + sout[109] + sout[125]
   + xorout[93] >= 0
 R3690: - sout[69] + sout[77] + sout[85] + sout[101] + sout[109]
   + sout[125] + xorout[93] >= 0
 R3691: sout[69] + sout[77] + sout[85] + sout[101] + sout[109] + sout[125]
   - xorout[93] >= 0
 R3692: sout[69] + sout[77] + sout[85] - sout[101] - sout[109] - sout[125]
   + xorout[93] >= -2
 R3693: sout[69] + sout[77] - sout[85] + sout[101] - sout[109] - sout[125]
   + xorout[93] >= -2
 R3694: sout[69] - sout[77] + sout[85] + sout[101] - sout[109] - sout[125]
   + xorout[93] >= -2
 R3695: - sout[69] + sout[77] + sout[85] + sout[101] - sout[109]
   - sout[125] + xorout[93] >= -2
 R3696: sout[69] + sout[77] - sout[85] - sout[101] + sout[109] - sout[125]
   + xorout[93] >= -2
 R3697: sout[69] - sout[77] + sout[85] - sout[101] + sout[109] - sout[125]
   + xorout[93] >= -2
 R3698: - sout[69] + sout[77] + sout[85] - sout[101] + sout[109]
   - sout[125] + xorout[93] >= -2
 R3699: sout[69] - sout[77] - sout[85] + sout[101] + sout[109] - sout[125]
   + xorout[93] >= -2
 R3700: - sout[69] + sout[77] - sout[85] + sout[101] + sout[109]
   - sout[125] + xorout[93] >= -2
 R3701: - sout[69] - sout[77] + sout[85] + sout[101] + sout[109]
   - sout[125] + xorout[93] >= -2
 R3702: sout[69] + sout[77] - sout[85] - sout[101] - sout[109] + sout[125]
   + xorout[93] >= -2
 R3703: sout[69] - sout[77] + sout[85] - sout[101] - sout[109] + sout[125]
   + xorout[93] >= -2
 R3704: - sout[69] + sout[77] + sout[85] - sout[101] - sout[109]
   + sout[125] + xorout[93] >= -2
 R3705: sout[69] - sout[77] - sout[85] + sout[101] - sout[109] + sout[125]
   + xorout[93] >= -2
 R3706: - sout[69] + sout[77] - sout[85] + sout[101] - sout[109]
   + sout[125] + xorout[93] >= -2
 R3707: - sout[69] - sout[77] + sout[85] + sout[101] - sout[109]
   + sout[125] + xorout[93] >= -2
 R3708: sout[69] - sout[77] - sout[85] - sout[101] + sout[109] + sout[125]
   + xorout[93] >= -2
 R3709: - sout[69] + sout[77] - sout[85] - sout[101] + sout[109]
   + sout[125] + xorout[93] >= -2
 R3710: - sout[69] - sout[77] + sout[85] - sout[101] + sout[109]
   + sout[125] + xorout[93] >= -2
 R3711: - sout[69] - sout[77] - sout[85] + sout[101] + sout[109]
   + sout[125] + xorout[93] >= -2
 R3712: sout[69] - sout[77] - sout[85] - sout[101] - sout[109] - sout[125]
   + xorout[93] >= -4
 R3713: - sout[69] + sout[77] - sout[85] - sout[101] - sout[109]
   - sout[125] + xorout[93] >= -4
 R3714: - sout[69] - sout[77] + sout[85] - sout[101] - sout[109]
   - sout[125] + xorout[93] >= -4
 R3715: - sout[69] - sout[77] - sout[85] + sout[101] - sout[109]
   - sout[125] + xorout[93] >= -4
 R3716: - sout[69] - sout[77] - sout[85] - sout[101] + sout[109]
   - sout[125] + xorout[93] >= -4
 R3717: - sout[69] - sout[77] - sout[85] - sout[101] - sout[109]
   + sout[125] + xorout[93] >= -4
 R3718: - sout[124] + xorout[92] >= 0
 R3719: - sout[108] + xorout[92] >= 0
 R3720: - sout[100] + xorout[92] >= 0
 R3721: - sout[84] + xorout[92] >= 0
 R3722: - sout[76] + xorout[92] >= 0
 R3723: - sout[68] + xorout[92] >= 0
 R3724: - sout[68] + sout[69] >= 0
 R3725: - sout[76] + sout[77] >= 0
 R3726: - sout[84] + sout[85] >= 0
 R3727: - sout[100] + sout[101] >= 0
 R3728: - sout[108] + sout[109] >= 0
 R3729: - sout[124] + sout[125] >= 0
 R3730: sout[68] + sout[76] + sout[84] + sout[100] + sout[108] + sout[124]
   - xorout[92] >= 0
 R3731: - xorout[92] + xorout[93] >= 0
 R3732: - sout[69] - sout[77] - sout[85] - sout[101] - sout[109]
   - sout[125] + xorout[92] - xorout[93] >= -6
 R3733: sout[69] + sout[77] - sout[85] - sout[101] - sout[109] - sout[125]
   + xorout[92] - xorout[93] >= -4
 R3734: sout[69] - sout[77] + sout[85] - sout[101] - sout[109] - sout[125]
   + xorout[92] - xorout[93] >= -4
 R3735: - sout[69] + sout[77] + sout[85] - sout[101] - sout[109]
   - sout[125] + xorout[92] - xorout[93] >= -4
 R3736: sout[69] - sout[77] - sout[85] + sout[101] - sout[109] - sout[125]
   + xorout[92] - xorout[93] >= -4
 R3737: - sout[69] + sout[77] - sout[85] + sout[101] - sout[109]
   - sout[125] + xorout[92] - xorout[93] >= -4
 R3738: - sout[69] - sout[77] + sout[85] + sout[101] - sout[109]
   - sout[125] + xorout[92] - xorout[93] >= -4
 R3739: sout[69] - sout[77] - sout[85] - sout[101] + sout[109] - sout[125]
   + xorout[92] - xorout[93] >= -4
 R3740: - sout[69] + sout[77] - sout[85] - sout[101] + sout[109]
   - sout[125] + xorout[92] - xorout[93] >= -4
 R3741: - sout[69] - sout[77] + sout[85] - sout[101] + sout[109]
   - sout[125] + xorout[92] - xorout[93] >= -4
 R3742: - sout[69] - sout[77] - sout[85] + sout[101] + sout[109]
   - sout[125] + xorout[92] - xorout[93] >= -4
 R3743: sout[69] - sout[77] - sout[85] - sout[101] - sout[109] + sout[125]
   + xorout[92] - xorout[93] >= -4
 R3744: - sout[69] + sout[77] - sout[85] - sout[101] - sout[109]
   + sout[125] + xorout[92] - xorout[93] >= -4
 R3745: - sout[69] - sout[77] + sout[85] - sout[101] - sout[109]
   + sout[125] + xorout[92] - xorout[93] >= -4
 R3746: - sout[69] - sout[77] - sout[85] + sout[101] - sout[109]
   + sout[125] + xorout[92] - xorout[93] >= -4
 R3747: - sout[69] - sout[77] - sout[85] - sout[101] + sout[109]
   + sout[125] + xorout[92] - xorout[93] >= -4
 R3748: sout[69] + sout[77] + sout[85] + sout[101] - sout[109] - sout[125]
   + xorout[92] - xorout[93] >= -2
 R3749: sout[69] + sout[77] + sout[85] - sout[101] + sout[109] - sout[125]
   + xorout[92] - xorout[93] >= -2
 R3750: sout[69] + sout[77] - sout[85] + sout[101] + sout[109] - sout[125]
   + xorout[92] - xorout[93] >= -2
 R3751: sout[69] - sout[77] + sout[85] + sout[101] + sout[109] - sout[125]
   + xorout[92] - xorout[93] >= -2
 R3752: - sout[69] + sout[77] + sout[85] + sout[101] + sout[109]
   - sout[125] + xorout[92] - xorout[93] >= -2
 R3753: sout[69] + sout[77] + sout[85] - sout[101] - sout[109] + sout[125]
   + xorout[92] - xorout[93] >= -2
 R3754: sout[69] + sout[77] - sout[85] + sout[101] - sout[109] + sout[125]
   + xorout[92] - xorout[93] >= -2
 R3755: sout[69] - sout[77] + sout[85] + sout[101] - sout[109] + sout[125]
   + xorout[92] - xorout[93] >= -2
 R3756: - sout[69] + sout[77] + sout[85] + sout[101] - sout[109]
   + sout[125] + xorout[92] - xorout[93] >= -2
 R3757: sout[69] + sout[77] - sout[85] - sout[101] + sout[109] + sout[125]
   + xorout[92] - xorout[93] >= -2
 R3758: sout[69] - sout[77] + sout[85] - sout[101] + sout[109] + sout[125]
   + xorout[92] - xorout[93] >= -2
 R3759: - sout[69] + sout[77] + sout[85] - sout[101] + sout[109]
   + sout[125] + xorout[92] - xorout[93] >= -2
 R3760: sout[69] - sout[77] - sout[85] + sout[101] + sout[109] + sout[125]
   + xorout[92] - xorout[93] >= -2
 R3761: - sout[69] + sout[77] - sout[85] + sout[101] + sout[109]
   + sout[125] + xorout[92] - xorout[93] >= -2
 R3762: - sout[69] - sout[77] + sout[85] + sout[101] + sout[109]
   + sout[125] + xorout[92] - xorout[93] >= -2
 R3763: sout[71] + sout[79] + sout[87] + sout[103] + sout[111] - sout[127]
   + xorout[95] >= 0
 R3764: sout[71] + sout[79] + sout[87] + sout[103] - sout[111] + sout[127]
   + xorout[95] >= 0
 R3765: sout[71] + sout[79] + sout[87] - sout[103] + sout[111] + sout[127]
   + xorout[95] >= 0
 R3766: sout[71] + sout[79] - sout[87] + sout[103] + sout[111] + sout[127]
   + xorout[95] >= 0
 R3767: sout[71] - sout[79] + sout[87] + sout[103] + sout[111] + sout[127]
   + xorout[95] >= 0
 R3768: - sout[71] + sout[79] + sout[87] + sout[103] + sout[111]
   + sout[127] + xorout[95] >= 0
 R3769: sout[71] + sout[79] + sout[87] + sout[103] + sout[111] + sout[127]
   - xorout[95] >= 0
 R3770: sout[71] + sout[79] + sout[87] - sout[103] - sout[111] - sout[127]
   + xorout[95] >= -2
 R3771: sout[71] + sout[79] - sout[87] + sout[103] - sout[111] - sout[127]
   + xorout[95] >= -2
 R3772: sout[71] - sout[79] + sout[87] + sout[103] - sout[111] - sout[127]
   + xorout[95] >= -2
 R3773: - sout[71] + sout[79] + sout[87] + sout[103] - sout[111]
   - sout[127] + xorout[95] >= -2
 R3774: sout[71] + sout[79] - sout[87] - sout[103] + sout[111] - sout[127]
   + xorout[95] >= -2
 R3775: sout[71] - sout[79] + sout[87] - sout[103] + sout[111] - sout[127]
   + xorout[95] >= -2
 R3776: - sout[71] + sout[79] + sout[87] - sout[103] + sout[111]
   - sout[127] + xorout[95] >= -2
 R3777: sout[71] - sout[79] - sout[87] + sout[103] + sout[111] - sout[127]
   + xorout[95] >= -2
 R3778: - sout[71] + sout[79] - sout[87] + sout[103] + sout[111]
   - sout[127] + xorout[95] >= -2
 R3779: - sout[71] - sout[79] + sout[87] + sout[103] + sout[111]
   - sout[127] + xorout[95] >= -2
 R3780: sout[71] + sout[79] - sout[87] - sout[103] - sout[111] + sout[127]
   + xorout[95] >= -2
 R3781: sout[71] - sout[79] + sout[87] - sout[103] - sout[111] + sout[127]
   + xorout[95] >= -2
 R3782: - sout[71] + sout[79] + sout[87] - sout[103] - sout[111]
   + sout[127] + xorout[95] >= -2
 R3783: sout[71] - sout[79] - sout[87] + sout[103] - sout[111] + sout[127]
   + xorout[95] >= -2
 R3784: - sout[71] + sout[79] - sout[87] + sout[103] - sout[111]
   + sout[127] + xorout[95] >= -2
 R3785: - sout[71] - sout[79] + sout[87] + sout[103] - sout[111]
   + sout[127] + xorout[95] >= -2
 R3786: sout[71] - sout[79] - sout[87] - sout[103] + sout[111] + sout[127]
   + xorout[95] >= -2
 R3787: - sout[71] + sout[79] - sout[87] - sout[103] + sout[111]
   + sout[127] + xorout[95] >= -2
 R3788: - sout[71] - sout[79] + sout[87] - sout[103] + sout[111]
   + sout[127] + xorout[95] >= -2
 R3789: - sout[71] - sout[79] - sout[87] + sout[103] + sout[111]
   + sout[127] + xorout[95] >= -2
 R3790: sout[71] - sout[79] - sout[87] - sout[103] - sout[111] - sout[127]
   + xorout[95] >= -4
 R3791: - sout[71] + sout[79] - sout[87] - sout[103] - sout[111]
   - sout[127] + xorout[95] >= -4
 R3792: - sout[71] - sout[79] + sout[87] - sout[103] - sout[111]
   - sout[127] + xorout[95] >= -4
 R3793: - sout[71] - sout[79] - sout[87] + sout[103] - sout[111]
   - sout[127] + xorout[95] >= -4
 R3794: - sout[71] - sout[79] - sout[87] - sout[103] + sout[111]
   - sout[127] + xorout[95] >= -4
 R3795: - sout[71] - sout[79] - sout[87] - sout[103] - sout[111]
   + sout[127] + xorout[95] >= -4
 R3796: - sout[126] + xorout[94] >= 0
 R3797: - sout[110] + xorout[94] >= 0
 R3798: - sout[102] + xorout[94] >= 0
 R3799: - sout[86] + xorout[94] >= 0
 R3800: - sout[78] + xorout[94] >= 0
 R3801: - sout[70] + xorout[94] >= 0
 R3802: - sout[70] + sout[71] >= 0
 R3803: - sout[78] + sout[79] >= 0
 R3804: - sout[86] + sout[87] >= 0
 R3805: - sout[102] + sout[103] >= 0
 R3806: - sout[110] + sout[111] >= 0
 R3807: - sout[126] + sout[127] >= 0
 R3808: sout[70] + sout[78] + sout[86] + sout[102] + sout[110] + sout[126]
   - xorout[94] >= 0
 R3809: - xorout[94] + xorout[95] >= 0
 R3810: - sout[71] - sout[79] - sout[87] - sout[103] - sout[111]
   - sout[127] + xorout[94] - xorout[95] >= -6
 R3811: sout[71] + sout[79] - sout[87] - sout[103] - sout[111] - sout[127]
   + xorout[94] - xorout[95] >= -4
 R3812: sout[71] - sout[79] + sout[87] - sout[103] - sout[111] - sout[127]
   + xorout[94] - xorout[95] >= -4
 R3813: - sout[71] + sout[79] + sout[87] - sout[103] - sout[111]
   - sout[127] + xorout[94] - xorout[95] >= -4
 R3814: sout[71] - sout[79] - sout[87] + sout[103] - sout[111] - sout[127]
   + xorout[94] - xorout[95] >= -4
 R3815: - sout[71] + sout[79] - sout[87] + sout[103] - sout[111]
   - sout[127] + xorout[94] - xorout[95] >= -4
 R3816: - sout[71] - sout[79] + sout[87] + sout[103] - sout[111]
   - sout[127] + xorout[94] - xorout[95] >= -4
 R3817: sout[71] - sout[79] - sout[87] - sout[103] + sout[111] - sout[127]
   + xorout[94] - xorout[95] >= -4
 R3818: - sout[71] + sout[79] - sout[87] - sout[103] + sout[111]
   - sout[127] + xorout[94] - xorout[95] >= -4
 R3819: - sout[71] - sout[79] + sout[87] - sout[103] + sout[111]
   - sout[127] + xorout[94] - xorout[95] >= -4
 R3820: - sout[71] - sout[79] - sout[87] + sout[103] + sout[111]
   - sout[127] + xorout[94] - xorout[95] >= -4
 R3821: sout[71] - sout[79] - sout[87] - sout[103] - sout[111] + sout[127]
   + xorout[94] - xorout[95] >= -4
 R3822: - sout[71] + sout[79] - sout[87] - sout[103] - sout[111]
   + sout[127] + xorout[94] - xorout[95] >= -4
 R3823: - sout[71] - sout[79] + sout[87] - sout[103] - sout[111]
   + sout[127] + xorout[94] - xorout[95] >= -4
 R3824: - sout[71] - sout[79] - sout[87] + sout[103] - sout[111]
   + sout[127] + xorout[94] - xorout[95] >= -4
 R3825: - sout[71] - sout[79] - sout[87] - sout[103] + sout[111]
   + sout[127] + xorout[94] - xorout[95] >= -4
 R3826: sout[71] + sout[79] + sout[87] + sout[103] - sout[111] - sout[127]
   + xorout[94] - xorout[95] >= -2
 R3827: sout[71] + sout[79] + sout[87] - sout[103] + sout[111] - sout[127]
   + xorout[94] - xorout[95] >= -2
 R3828: sout[71] + sout[79] - sout[87] + sout[103] + sout[111] - sout[127]
   + xorout[94] - xorout[95] >= -2
 R3829: sout[71] - sout[79] + sout[87] + sout[103] + sout[111] - sout[127]
   + xorout[94] - xorout[95] >= -2
 R3830: - sout[71] + sout[79] + sout[87] + sout[103] + sout[111]
   - sout[127] + xorout[94] - xorout[95] >= -2
 R3831: sout[71] + sout[79] + sout[87] - sout[103] - sout[111] + sout[127]
   + xorout[94] - xorout[95] >= -2
 R3832: sout[71] + sout[79] - sout[87] + sout[103] - sout[111] + sout[127]
   + xorout[94] - xorout[95] >= -2
 R3833: sout[71] - sout[79] + sout[87] + sout[103] - sout[111] + sout[127]
   + xorout[94] - xorout[95] >= -2
 R3834: - sout[71] + sout[79] + sout[87] + sout[103] - sout[111]
   + sout[127] + xorout[94] - xorout[95] >= -2
 R3835: sout[71] + sout[79] - sout[87] - sout[103] + sout[111] + sout[127]
   + xorout[94] - xorout[95] >= -2
 R3836: sout[71] - sout[79] + sout[87] - sout[103] + sout[111] + sout[127]
   + xorout[94] - xorout[95] >= -2
 R3837: - sout[71] + sout[79] + sout[87] - sout[103] + sout[111]
   + sout[127] + xorout[94] - xorout[95] >= -2
 R3838: sout[71] - sout[79] - sout[87] + sout[103] + sout[111] + sout[127]
   + xorout[94] - xorout[95] >= -2
 R3839: - sout[71] + sout[79] - sout[87] + sout[103] + sout[111]
   + sout[127] + xorout[94] - xorout[95] >= -2
 R3840: - sout[71] - sout[79] + sout[87] + sout[103] + sout[111]
   + sout[127] + xorout[94] - xorout[95] >= -2
 R3841: sout[65] + sout[73] + sout[89] + sout[97] - sout[105] + xorout[97]
   >= 0
 R3842: sout[65] + sout[73] + sout[89] - sout[97] + sout[105] + xorout[97]
   >= 0
 R3843: sout[65] + sout[73] - sout[89] + sout[97] + sout[105] + xorout[97]
   >= 0
 R3844: sout[65] - sout[73] + sout[89] + sout[97] + sout[105] + xorout[97]
   >= 0
 R3845: - sout[65] + sout[73] + sout[89] + sout[97] + sout[105]
   + xorout[97] >= 0
 R3846: sout[65] + sout[73] + sout[89] + sout[97] + sout[105] - xorout[97]
   >= 0
 R3847: sout[65] + sout[73] - sout[89] - sout[97] - sout[105] + xorout[97]
   >= -2
 R3848: sout[65] - sout[73] + sout[89] - sout[97] - sout[105] + xorout[97]
   >= -2
 R3849: - sout[65] + sout[73] + sout[89] - sout[97] - sout[105]
   + xorout[97] >= -2
 R3850: sout[65] - sout[73] - sout[89] + sout[97] - sout[105] + xorout[97]
   >= -2
 R3851: - sout[65] + sout[73] - sout[89] + sout[97] - sout[105]
   + xorout[97] >= -2
 R3852: - sout[65] - sout[73] + sout[89] + sout[97] - sout[105]
   + xorout[97] >= -2
 R3853: sout[65] - sout[73] - sout[89] - sout[97] + sout[105] + xorout[97]
   >= -2
 R3854: - sout[65] + sout[73] - sout[89] - sout[97] + sout[105]
   + xorout[97] >= -2
 R3855: - sout[65] - sout[73] + sout[89] - sout[97] + sout[105]
   + xorout[97] >= -2
 R3856: - sout[65] - sout[73] - sout[89] + sout[97] + sout[105]
   + xorout[97] >= -2
 R3857: - sout[65] - sout[73] - sout[89] - sout[97] - sout[105]
   + xorout[97] >= -4
 R3858: - sout[104] + xorout[96] >= 0
 R3859: - sout[96] + xorout[96] >= 0
 R3860: - sout[88] + xorout[96] >= 0
 R3861: - sout[72] + xorout[96] >= 0
 R3862: - sout[64] + xorout[96] >= 0
 R3863: - sout[64] + sout[65] >= 0
 R3864: - sout[72] + sout[73] >= 0
 R3865: - sout[88] + sout[89] >= 0
 R3866: - sout[96] + sout[97] >= 0
 R3867: - sout[104] + sout[105] >= 0
 R3868: sout[64] + sout[72] + sout[88] + sout[96] + sout[104] - xorout[96]
   >= 0
 R3869: - xorout[96] + xorout[97] >= 0
 R3870: sout[65] - sout[73] - sout[89] - sout[97] - sout[105] + xorout[96]
   - xorout[97] >= -4
 R3871: - sout[65] + sout[73] - sout[89] - sout[97] - sout[105]
   + xorout[96] - xorout[97] >= -4
 R3872: - sout[65] - sout[73] + sout[89] - sout[97] - sout[105]
   + xorout[96] - xorout[97] >= -4
 R3873: - sout[65] - sout[73] - sout[89] + sout[97] - sout[105]
   + xorout[96] - xorout[97] >= -4
 R3874: - sout[65] - sout[73] - sout[89] - sout[97] + sout[105]
   + xorout[96] - xorout[97] >= -4
 R3875: sout[65] + sout[73] + sout[89] - sout[97] - sout[105] + xorout[96]
   - xorout[97] >= -2
 R3876: sout[65] + sout[73] - sout[89] + sout[97] - sout[105] + xorout[96]
   - xorout[97] >= -2
 R3877: sout[65] - sout[73] + sout[89] + sout[97] - sout[105] + xorout[96]
   - xorout[97] >= -2
 R3878: - sout[65] + sout[73] + sout[89] + sout[97] - sout[105]
   + xorout[96] - xorout[97] >= -2
 R3879: sout[65] + sout[73] - sout[89] - sout[97] + sout[105] + xorout[96]
   - xorout[97] >= -2
 R3880: sout[65] - sout[73] + sout[89] - sout[97] + sout[105] + xorout[96]
   - xorout[97] >= -2
 R3881: - sout[65] + sout[73] + sout[89] - sout[97] + sout[105]
   + xorout[96] - xorout[97] >= -2
 R3882: sout[65] - sout[73] - sout[89] + sout[97] + sout[105] + xorout[96]
   - xorout[97] >= -2
 R3883: - sout[65] + sout[73] - sout[89] + sout[97] + sout[105]
   + xorout[96] - xorout[97] >= -2
 R3884: - sout[65] - sout[73] + sout[89] + sout[97] + sout[105]
   + xorout[96] - xorout[97] >= -2
 R3885: sout[67] + sout[75] + sout[91] + sout[99] - sout[107] + xorout[99]
   >= 0
 R3886: sout[67] + sout[75] + sout[91] - sout[99] + sout[107] + xorout[99]
   >= 0
 R3887: sout[67] + sout[75] - sout[91] + sout[99] + sout[107] + xorout[99]
   >= 0
 R3888: sout[67] - sout[75] + sout[91] + sout[99] + sout[107] + xorout[99]
   >= 0
 R3889: - sout[67] + sout[75] + sout[91] + sout[99] + sout[107]
   + xorout[99] >= 0
 R3890: sout[67] + sout[75] + sout[91] + sout[99] + sout[107] - xorout[99]
   >= 0
 R3891: sout[67] + sout[75] - sout[91] - sout[99] - sout[107] + xorout[99]
   >= -2
 R3892: sout[67] - sout[75] + sout[91] - sout[99] - sout[107] + xorout[99]
   >= -2
 R3893: - sout[67] + sout[75] + sout[91] - sout[99] - sout[107]
   + xorout[99] >= -2
 R3894: sout[67] - sout[75] - sout[91] + sout[99] - sout[107] + xorout[99]
   >= -2
 R3895: - sout[67] + sout[75] - sout[91] + sout[99] - sout[107]
   + xorout[99] >= -2
 R3896: - sout[67] - sout[75] + sout[91] + sout[99] - sout[107]
   + xorout[99] >= -2
 R3897: sout[67] - sout[75] - sout[91] - sout[99] + sout[107] + xorout[99]
   >= -2
 R3898: - sout[67] + sout[75] - sout[91] - sout[99] + sout[107]
   + xorout[99] >= -2
 R3899: - sout[67] - sout[75] + sout[91] - sout[99] + sout[107]
   + xorout[99] >= -2
 R3900: - sout[67] - sout[75] - sout[91] + sout[99] + sout[107]
   + xorout[99] >= -2
 R3901: - sout[67] - sout[75] - sout[91] - sout[99] - sout[107]
   + xorout[99] >= -4
 R3902: - sout[106] + xorout[98] >= 0
 R3903: - sout[98] + xorout[98] >= 0
 R3904: - sout[90] + xorout[98] >= 0
 R3905: - sout[74] + xorout[98] >= 0
 R3906: - sout[66] + xorout[98] >= 0
 R3907: - sout[66] + sout[67] >= 0
 R3908: - sout[74] + sout[75] >= 0
 R3909: - sout[90] + sout[91] >= 0
 R3910: - sout[98] + sout[99] >= 0
 R3911: - sout[106] + sout[107] >= 0
 R3912: sout[66] + sout[74] + sout[90] + sout[98] + sout[106] - xorout[98]
   >= 0
 R3913: - xorout[98] + xorout[99] >= 0
 R3914: sout[67] - sout[75] - sout[91] - sout[99] - sout[107] + xorout[98]
   - xorout[99] >= -4
 R3915: - sout[67] + sout[75] - sout[91] - sout[99] - sout[107]
   + xorout[98] - xorout[99] >= -4
 R3916: - sout[67] - sout[75] + sout[91] - sout[99] - sout[107]
   + xorout[98] - xorout[99] >= -4
 R3917: - sout[67] - sout[75] - sout[91] + sout[99] - sout[107]
   + xorout[98] - xorout[99] >= -4
 R3918: - sout[67] - sout[75] - sout[91] - sout[99] + sout[107]
   + xorout[98] - xorout[99] >= -4
 R3919: sout[67] + sout[75] + sout[91] - sout[99] - sout[107] + xorout[98]
   - xorout[99] >= -2
 R3920: sout[67] + sout[75] - sout[91] + sout[99] - sout[107] + xorout[98]
   - xorout[99] >= -2
 R3921: sout[67] - sout[75] + sout[91] + sout[99] - sout[107] + xorout[98]
   - xorout[99] >= -2
 R3922: - sout[67] + sout[75] + sout[91] + sout[99] - sout[107]
   + xorout[98] - xorout[99] >= -2
 R3923: sout[67] + sout[75] - sout[91] - sout[99] + sout[107] + xorout[98]
   - xorout[99] >= -2
 R3924: sout[67] - sout[75] + sout[91] - sout[99] + sout[107] + xorout[98]
   - xorout[99] >= -2
 R3925: - sout[67] + sout[75] + sout[91] - sout[99] + sout[107]
   + xorout[98] - xorout[99] >= -2
 R3926: sout[67] - sout[75] - sout[91] + sout[99] + sout[107] + xorout[98]
   - xorout[99] >= -2
 R3927: - sout[67] + sout[75] - sout[91] + sout[99] + sout[107]
   + xorout[98] - xorout[99] >= -2
 R3928: - sout[67] - sout[75] + sout[91] + sout[99] + sout[107]
   + xorout[98] - xorout[99] >= -2
 R3929: sout[69] + sout[77] + sout[93] + sout[101] - sout[109]
   + xorout[101] >= 0
 R3930: sout[69] + sout[77] + sout[93] - sout[101] + sout[109]
   + xorout[101] >= 0
 R3931: sout[69] + sout[77] - sout[93] + sout[101] + sout[109]
   + xorout[101] >= 0
 R3932: sout[69] - sout[77] + sout[93] + sout[101] + sout[109]
   + xorout[101] >= 0
 R3933: - sout[69] + sout[77] + sout[93] + sout[101] + sout[109]
   + xorout[101] >= 0
 R3934: sout[69] + sout[77] + sout[93] + sout[101] + sout[109]
   - xorout[101] >= 0
 R3935: sout[69] + sout[77] - sout[93] - sout[101] - sout[109]
   + xorout[101] >= -2
 R3936: sout[69] - sout[77] + sout[93] - sout[101] - sout[109]
   + xorout[101] >= -2
 R3937: - sout[69] + sout[77] + sout[93] - sout[101] - sout[109]
   + xorout[101] >= -2
 R3938: sout[69] - sout[77] - sout[93] + sout[101] - sout[109]
   + xorout[101] >= -2
 R3939: - sout[69] + sout[77] - sout[93] + sout[101] - sout[109]
   + xorout[101] >= -2
 R3940: - sout[69] - sout[77] + sout[93] + sout[101] - sout[109]
   + xorout[101] >= -2
 R3941: sout[69] - sout[77] - sout[93] - sout[101] + sout[109]
   + xorout[101] >= -2
 R3942: - sout[69] + sout[77] - sout[93] - sout[101] + sout[109]
   + xorout[101] >= -2
 R3943: - sout[69] - sout[77] + sout[93] - sout[101] + sout[109]
   + xorout[101] >= -2
 R3944: - sout[69] - sout[77] - sout[93] + sout[101] + sout[109]
   + xorout[101] >= -2
 R3945: - sout[69] - sout[77] - sout[93] - sout[101] - sout[109]
   + xorout[101] >= -4
 R3946: - sout[108] + xorout[100] >= 0
 R3947: - sout[100] + xorout[100] >= 0
 R3948: - sout[92] + xorout[100] >= 0
 R3949: - sout[76] + xorout[100] >= 0
 R3950: - sout[68] + xorout[100] >= 0
 R3951: - sout[68] + sout[69] >= 0
 R3952: - sout[76] + sout[77] >= 0
 R3953: - sout[92] + sout[93] >= 0
 R3954: - sout[100] + sout[101] >= 0
 R3955: - sout[108] + sout[109] >= 0
 R3956: sout[68] + sout[76] + sout[92] + sout[100] + sout[108]
   - xorout[100] >= 0
 R3957: - xorout[100] + xorout[101] >= 0
 R3958: sout[69] - sout[77] - sout[93] - sout[101] - sout[109]
   + xorout[100] - xorout[101] >= -4
 R3959: - sout[69] + sout[77] - sout[93] - sout[101] - sout[109]
   + xorout[100] - xorout[101] >= -4
 R3960: - sout[69] - sout[77] + sout[93] - sout[101] - sout[109]
   + xorout[100] - xorout[101] >= -4
 R3961: - sout[69] - sout[77] - sout[93] + sout[101] - sout[109]
   + xorout[100] - xorout[101] >= -4
 R3962: - sout[69] - sout[77] - sout[93] - sout[101] + sout[109]
   + xorout[100] - xorout[101] >= -4
 R3963: sout[69] + sout[77] + sout[93] - sout[101] - sout[109]
   + xorout[100] - xorout[101] >= -2
 R3964: sout[69] + sout[77] - sout[93] + sout[101] - sout[109]
   + xorout[100] - xorout[101] >= -2
 R3965: sout[69] - sout[77] + sout[93] + sout[101] - sout[109]
   + xorout[100] - xorout[101] >= -2
 R3966: - sout[69] + sout[77] + sout[93] + sout[101] - sout[109]
   + xorout[100] - xorout[101] >= -2
 R3967: sout[69] + sout[77] - sout[93] - sout[101] + sout[109]
   + xorout[100] - xorout[101] >= -2
 R3968: sout[69] - sout[77] + sout[93] - sout[101] + sout[109]
   + xorout[100] - xorout[101] >= -2
 R3969: - sout[69] + sout[77] + sout[93] - sout[101] + sout[109]
   + xorout[100] - xorout[101] >= -2
 R3970: sout[69] - sout[77] - sout[93] + sout[101] + sout[109]
   + xorout[100] - xorout[101] >= -2
 R3971: - sout[69] + sout[77] - sout[93] + sout[101] + sout[109]
   + xorout[100] - xorout[101] >= -2
 R3972: - sout[69] - sout[77] + sout[93] + sout[101] + sout[109]
   + xorout[100] - xorout[101] >= -2
 R3973: sout[71] + sout[79] + sout[95] + sout[103] - sout[111]
   + xorout[103] >= 0
 R3974: sout[71] + sout[79] + sout[95] - sout[103] + sout[111]
   + xorout[103] >= 0
 R3975: sout[71] + sout[79] - sout[95] + sout[103] + sout[111]
   + xorout[103] >= 0
 R3976: sout[71] - sout[79] + sout[95] + sout[103] + sout[111]
   + xorout[103] >= 0
 R3977: - sout[71] + sout[79] + sout[95] + sout[103] + sout[111]
   + xorout[103] >= 0
 R3978: sout[71] + sout[79] + sout[95] + sout[103] + sout[111]
   - xorout[103] >= 0
 R3979: sout[71] + sout[79] - sout[95] - sout[103] - sout[111]
   + xorout[103] >= -2
 R3980: sout[71] - sout[79] + sout[95] - sout[103] - sout[111]
   + xorout[103] >= -2
 R3981: - sout[71] + sout[79] + sout[95] - sout[103] - sout[111]
   + xorout[103] >= -2
 R3982: sout[71] - sout[79] - sout[95] + sout[103] - sout[111]
   + xorout[103] >= -2
 R3983: - sout[71] + sout[79] - sout[95] + sout[103] - sout[111]
   + xorout[103] >= -2
 R3984: - sout[71] - sout[79] + sout[95] + sout[103] - sout[111]
   + xorout[103] >= -2
 R3985: sout[71] - sout[79] - sout[95] - sout[103] + sout[111]
   + xorout[103] >= -2
 R3986: - sout[71] + sout[79] - sout[95] - sout[103] + sout[111]
   + xorout[103] >= -2
 R3987: - sout[71] - sout[79] + sout[95] - sout[103] + sout[111]
   + xorout[103] >= -2
 R3988: - sout[71] - sout[79] - sout[95] + sout[103] + sout[111]
   + xorout[103] >= -2
 R3989: - sout[71] - sout[79] - sout[95] - sout[103] - sout[111]
   + xorout[103] >= -4
 R3990: - sout[110] + xorout[102] >= 0
 R3991: - sout[102] + xorout[102] >= 0
 R3992: - sout[94] + xorout[102] >= 0
 R3993: - sout[78] + xorout[102] >= 0
 R3994: - sout[70] + xorout[102] >= 0
 R3995: - sout[70] + sout[71] >= 0
 R3996: - sout[78] + sout[79] >= 0
 R3997: - sout[94] + sout[95] >= 0
 R3998: - sout[102] + sout[103] >= 0
 R3999: - sout[110] + sout[111] >= 0
 R4000: sout[70] + sout[78] + sout[94] + sout[102] + sout[110]
   - xorout[102] >= 0
 R4001: - xorout[102] + xorout[103] >= 0
 R4002: sout[71] - sout[79] - sout[95] - sout[103] - sout[111]
   + xorout[102] - xorout[103] >= -4
 R4003: - sout[71] + sout[79] - sout[95] - sout[103] - sout[111]
   + xorout[102] - xorout[103] >= -4
 R4004: - sout[71] - sout[79] + sout[95] - sout[103] - sout[111]
   + xorout[102] - xorout[103] >= -4
 R4005: - sout[71] - sout[79] - sout[95] + sout[103] - sout[111]
   + xorout[102] - xorout[103] >= -4
 R4006: - sout[71] - sout[79] - sout[95] - sout[103] + sout[111]
   + xorout[102] - xorout[103] >= -4
 R4007: sout[71] + sout[79] + sout[95] - sout[103] - sout[111]
   + xorout[102] - xorout[103] >= -2
 R4008: sout[71] + sout[79] - sout[95] + sout[103] - sout[111]
   + xorout[102] - xorout[103] >= -2
 R4009: sout[71] - sout[79] + sout[95] + sout[103] - sout[111]
   + xorout[102] - xorout[103] >= -2
 R4010: - sout[71] + sout[79] + sout[95] + sout[103] - sout[111]
   + xorout[102] - xorout[103] >= -2
 R4011: sout[71] + sout[79] - sout[95] - sout[103] + sout[111]
   + xorout[102] - xorout[103] >= -2
 R4012: sout[71] - sout[79] + sout[95] - sout[103] + sout[111]
   + xorout[102] - xorout[103] >= -2
 R4013: - sout[71] + sout[79] + sout[95] - sout[103] + sout[111]
   + xorout[102] - xorout[103] >= -2
 R4014: sout[71] - sout[79] - sout[95] + sout[103] + sout[111]
   + xorout[102] - xorout[103] >= -2
 R4015: - sout[71] + sout[79] - sout[95] + sout[103] + sout[111]
   + xorout[102] - xorout[103] >= -2
 R4016: - sout[71] - sout[79] + sout[95] + sout[103] + sout[111]
   + xorout[102] - xorout[103] >= -2
 R4017: sout[65] + sout[73] + sout[81] + sout[105] - sout[113]
   + xorout[105] >= 0
 R4018: sout[65] + sout[73] + sout[81] - sout[105] + sout[113]
   + xorout[105] >= 0
 R4019: sout[65] + sout[73] - sout[81] + sout[105] + sout[113]
   + xorout[105] >= 0
 R4020: sout[65] - sout[73] + sout[81] + sout[105] + sout[113]
   + xorout[105] >= 0
 R4021: - sout[65] + sout[73] + sout[81] + sout[105] + sout[113]
   + xorout[105] >= 0
 R4022: sout[65] + sout[73] + sout[81] + sout[105] + sout[113]
   - xorout[105] >= 0
 R4023: sout[65] + sout[73] - sout[81] - sout[105] - sout[113]
   + xorout[105] >= -2
 R4024: sout[65] - sout[73] + sout[81] - sout[105] - sout[113]
   + xorout[105] >= -2
 R4025: - sout[65] + sout[73] + sout[81] - sout[105] - sout[113]
   + xorout[105] >= -2
 R4026: sout[65] - sout[73] - sout[81] + sout[105] - sout[113]
   + xorout[105] >= -2
 R4027: - sout[65] + sout[73] - sout[81] + sout[105] - sout[113]
   + xorout[105] >= -2
 R4028: - sout[65] - sout[73] + sout[81] + sout[105] - sout[113]
   + xorout[105] >= -2
 R4029: sout[65] - sout[73] - sout[81] - sout[105] + sout[113]
   + xorout[105] >= -2
 R4030: - sout[65] + sout[73] - sout[81] - sout[105] + sout[113]
   + xorout[105] >= -2
 R4031: - sout[65] - sout[73] + sout[81] - sout[105] + sout[113]
   + xorout[105] >= -2
 R4032: - sout[65] - sout[73] - sout[81] + sout[105] + sout[113]
   + xorout[105] >= -2
 R4033: - sout[65] - sout[73] - sout[81] - sout[105] - sout[113]
   + xorout[105] >= -4
 R4034: - sout[112] + xorout[104] >= 0
 R4035: - sout[104] + xorout[104] >= 0
 R4036: - sout[80] + xorout[104] >= 0
 R4037: - sout[72] + xorout[104] >= 0
 R4038: - sout[64] + xorout[104] >= 0
 R4039: - sout[64] + sout[65] >= 0
 R4040: - sout[72] + sout[73] >= 0
 R4041: - sout[80] + sout[81] >= 0
 R4042: - sout[104] + sout[105] >= 0
 R4043: - sout[112] + sout[113] >= 0
 R4044: sout[64] + sout[72] + sout[80] + sout[104] + sout[112]
   - xorout[104] >= 0
 R4045: - xorout[104] + xorout[105] >= 0
 R4046: sout[65] - sout[73] - sout[81] - sout[105] - sout[113]
   + xorout[104] - xorout[105] >= -4
 R4047: - sout[65] + sout[73] - sout[81] - sout[105] - sout[113]
   + xorout[104] - xorout[105] >= -4
 R4048: - sout[65] - sout[73] + sout[81] - sout[105] - sout[113]
   + xorout[104] - xorout[105] >= -4
 R4049: - sout[65] - sout[73] - sout[81] + sout[105] - sout[113]
   + xorout[104] - xorout[105] >= -4
 R4050: - sout[65] - sout[73] - sout[81] - sout[105] + sout[113]
   + xorout[104] - xorout[105] >= -4
 R4051: sout[65] + sout[73] + sout[81] - sout[105] - sout[113]
   + xorout[104] - xorout[105] >= -2
 R4052: sout[65] + sout[73] - sout[81] + sout[105] - sout[113]
   + xorout[104] - xorout[105] >= -2
 R4053: sout[65] - sout[73] + sout[81] + sout[105] - sout[113]
   + xorout[104] - xorout[105] >= -2
 R4054: - sout[65] + sout[73] + sout[81] + sout[105] - sout[113]
   + xorout[104] - xorout[105] >= -2
 R4055: sout[65] + sout[73] - sout[81] - sout[105] + sout[113]
   + xorout[104] - xorout[105] >= -2
 R4056: sout[65] - sout[73] + sout[81] - sout[105] + sout[113]
   + xorout[104] - xorout[105] >= -2
 R4057: - sout[65] + sout[73] + sout[81] - sout[105] + sout[113]
   + xorout[104] - xorout[105] >= -2
 R4058: sout[65] - sout[73] - sout[81] + sout[105] + sout[113]
   + xorout[104] - xorout[105] >= -2
 R4059: - sout[65] + sout[73] - sout[81] + sout[105] + sout[113]
   + xorout[104] - xorout[105] >= -2
 R4060: - sout[65] - sout[73] + sout[81] + sout[105] + sout[113]
   + xorout[104] - xorout[105] >= -2
 R4061: sout[67] + sout[75] + sout[83] + sout[107] - sout[115]
   + xorout[107] >= 0
 R4062: sout[67] + sout[75] + sout[83] - sout[107] + sout[115]
   + xorout[107] >= 0
 R4063: sout[67] + sout[75] - sout[83] + sout[107] + sout[115]
   + xorout[107] >= 0
 R4064: sout[67] - sout[75] + sout[83] + sout[107] + sout[115]
   + xorout[107] >= 0
 R4065: - sout[67] + sout[75] + sout[83] + sout[107] + sout[115]
   + xorout[107] >= 0
 R4066: sout[67] + sout[75] + sout[83] + sout[107] + sout[115]
   - xorout[107] >= 0
 R4067: sout[67] + sout[75] - sout[83] - sout[107] - sout[115]
   + xorout[107] >= -2
 R4068: sout[67] - sout[75] + sout[83] - sout[107] - sout[115]
   + xorout[107] >= -2
 R4069: - sout[67] + sout[75] + sout[83] - sout[107] - sout[115]
   + xorout[107] >= -2
 R4070: sout[67] - sout[75] - sout[83] + sout[107] - sout[115]
   + xorout[107] >= -2
 R4071: - sout[67] + sout[75] - sout[83] + sout[107] - sout[115]
   + xorout[107] >= -2
 R4072: - sout[67] - sout[75] + sout[83] + sout[107] - sout[115]
   + xorout[107] >= -2
 R4073: sout[67] - sout[75] - sout[83] - sout[107] + sout[115]
   + xorout[107] >= -2
 R4074: - sout[67] + sout[75] - sout[83] - sout[107] + sout[115]
   + xorout[107] >= -2
 R4075: - sout[67] - sout[75] + sout[83] - sout[107] + sout[115]
   + xorout[107] >= -2
 R4076: - sout[67] - sout[75] - sout[83] + sout[107] + sout[115]
   + xorout[107] >= -2
 R4077: - sout[67] - sout[75] - sout[83] - sout[107] - sout[115]
   + xorout[107] >= -4
 R4078: - sout[114] + xorout[106] >= 0
 R4079: - sout[106] + xorout[106] >= 0
 R4080: - sout[82] + xorout[106] >= 0
 R4081: - sout[74] + xorout[106] >= 0
 R4082: - sout[66] + xorout[106] >= 0
 R4083: - sout[66] + sout[67] >= 0
 R4084: - sout[74] + sout[75] >= 0
 R4085: - sout[82] + sout[83] >= 0
 R4086: - sout[106] + sout[107] >= 0
 R4087: - sout[114] + sout[115] >= 0
 R4088: sout[66] + sout[74] + sout[82] + sout[106] + sout[114]
   - xorout[106] >= 0
 R4089: - xorout[106] + xorout[107] >= 0
 R4090: sout[67] - sout[75] - sout[83] - sout[107] - sout[115]
   + xorout[106] - xorout[107] >= -4
 R4091: - sout[67] + sout[75] - sout[83] - sout[107] - sout[115]
   + xorout[106] - xorout[107] >= -4
 R4092: - sout[67] - sout[75] + sout[83] - sout[107] - sout[115]
   + xorout[106] - xorout[107] >= -4
 R4093: - sout[67] - sout[75] - sout[83] + sout[107] - sout[115]
   + xorout[106] - xorout[107] >= -4
 R4094: - sout[67] - sout[75] - sout[83] - sout[107] + sout[115]
   + xorout[106] - xorout[107] >= -4
 R4095: sout[67] + sout[75] + sout[83] - sout[107] - sout[115]
   + xorout[106] - xorout[107] >= -2
 R4096: sout[67] + sout[75] - sout[83] + sout[107] - sout[115]
   + xorout[106] - xorout[107] >= -2
 R4097: sout[67] - sout[75] + sout[83] + sout[107] - sout[115]
   + xorout[106] - xorout[107] >= -2
 R4098: - sout[67] + sout[75] + sout[83] + sout[107] - sout[115]
   + xorout[106] - xorout[107] >= -2
 R4099: sout[67] + sout[75] - sout[83] - sout[107] + sout[115]
   + xorout[106] - xorout[107] >= -2
 R4100: sout[67] - sout[75] + sout[83] - sout[107] + sout[115]
   + xorout[106] - xorout[107] >= -2
 R4101: - sout[67] + sout[75] + sout[83] - sout[107] + sout[115]
   + xorout[106] - xorout[107] >= -2
 R4102: sout[67] - sout[75] - sout[83] + sout[107] + sout[115]
   + xorout[106] - xorout[107] >= -2
 R4103: - sout[67] + sout[75] - sout[83] + sout[107] + sout[115]
   + xorout[106] - xorout[107] >= -2
 R4104: - sout[67] - sout[75] + sout[83] + sout[107] + sout[115]
   + xorout[106] - xorout[107] >= -2
 R4105: sout[69] + sout[77] + sout[85] + sout[109] - sout[117]
   + xorout[109] >= 0
 R4106: sout[69] + sout[77] + sout[85] - sout[109] + sout[117]
   + xorout[109] >= 0
 R4107: sout[69] + sout[77] - sout[85] + sout[109] + sout[117]
   + xorout[109] >= 0
 R4108: sout[69] - sout[77] + sout[85] + sout[109] + sout[117]
   + xorout[109] >= 0
 R4109: - sout[69] + sout[77] + sout[85] + sout[109] + sout[117]
   + xorout[109] >= 0
 R4110: sout[69] + sout[77] + sout[85] + sout[109] + sout[117]
   - xorout[109] >= 0
 R4111: sout[69] + sout[77] - sout[85] - sout[109] - sout[117]
   + xorout[109] >= -2
 R4112: sout[69] - sout[77] + sout[85] - sout[109] - sout[117]
   + xorout[109] >= -2
 R4113: - sout[69] + sout[77] + sout[85] - sout[109] - sout[117]
   + xorout[109] >= -2
 R4114: sout[69] - sout[77] - sout[85] + sout[109] - sout[117]
   + xorout[109] >= -2
 R4115: - sout[69] + sout[77] - sout[85] + sout[109] - sout[117]
   + xorout[109] >= -2
 R4116: - sout[69] - sout[77] + sout[85] + sout[109] - sout[117]
   + xorout[109] >= -2
 R4117: sout[69] - sout[77] - sout[85] - sout[109] + sout[117]
   + xorout[109] >= -2
 R4118: - sout[69] + sout[77] - sout[85] - sout[109] + sout[117]
   + xorout[109] >= -2
 R4119: - sout[69] - sout[77] + sout[85] - sout[109] + sout[117]
   + xorout[109] >= -2
 R4120: - sout[69] - sout[77] - sout[85] + sout[109] + sout[117]
   + xorout[109] >= -2
 R4121: - sout[69] - sout[77] - sout[85] - sout[109] - sout[117]
   + xorout[109] >= -4
 R4122: - sout[116] + xorout[108] >= 0
 R4123: - sout[108] + xorout[108] >= 0
 R4124: - sout[84] + xorout[108] >= 0
 R4125: - sout[76] + xorout[108] >= 0
 R4126: - sout[68] + xorout[108] >= 0
 R4127: - sout[68] + sout[69] >= 0
 R4128: - sout[76] + sout[77] >= 0
 R4129: - sout[84] + sout[85] >= 0
 R4130: - sout[108] + sout[109] >= 0
 R4131: - sout[116] + sout[117] >= 0
 R4132: sout[68] + sout[76] + sout[84] + sout[108] + sout[116]
   - xorout[108] >= 0
 R4133: - xorout[108] + xorout[109] >= 0
 R4134: sout[69] - sout[77] - sout[85] - sout[109] - sout[117]
   + xorout[108] - xorout[109] >= -4
 R4135: - sout[69] + sout[77] - sout[85] - sout[109] - sout[117]
   + xorout[108] - xorout[109] >= -4
 R4136: - sout[69] - sout[77] + sout[85] - sout[109] - sout[117]
   + xorout[108] - xorout[109] >= -4
 R4137: - sout[69] - sout[77] - sout[85] + sout[109] - sout[117]
   + xorout[108] - xorout[109] >= -4
 R4138: - sout[69] - sout[77] - sout[85] - sout[109] + sout[117]
   + xorout[108] - xorout[109] >= -4
 R4139: sout[69] + sout[77] + sout[85] - sout[109] - sout[117]
   + xorout[108] - xorout[109] >= -2
 R4140: sout[69] + sout[77] - sout[85] + sout[109] - sout[117]
   + xorout[108] - xorout[109] >= -2
 R4141: sout[69] - sout[77] + sout[85] + sout[109] - sout[117]
   + xorout[108] - xorout[109] >= -2
 R4142: - sout[69] + sout[77] + sout[85] + sout[109] - sout[117]
   + xorout[108] - xorout[109] >= -2
 R4143: sout[69] + sout[77] - sout[85] - sout[109] + sout[117]
   + xorout[108] - xorout[109] >= -2
 R4144: sout[69] - sout[77] + sout[85] - sout[109] + sout[117]
   + xorout[108] - xorout[109] >= -2
 R4145: - sout[69] + sout[77] + sout[85] - sout[109] + sout[117]
   + xorout[108] - xorout[109] >= -2
 R4146: sout[69] - sout[77] - sout[85] + sout[109] + sout[117]
   + xorout[108] - xorout[109] >= -2
 R4147: - sout[69] + sout[77] - sout[85] + sout[109] + sout[117]
   + xorout[108] - xorout[109] >= -2
 R4148: - sout[69] - sout[77] + sout[85] + sout[109] + sout[117]
   + xorout[108] - xorout[109] >= -2
 R4149: sout[71] + sout[79] + sout[87] + sout[111] - sout[119]
   + xorout[111] >= 0
 R4150: sout[71] + sout[79] + sout[87] - sout[111] + sout[119]
   + xorout[111] >= 0
 R4151: sout[71] + sout[79] - sout[87] + sout[111] + sout[119]
   + xorout[111] >= 0
 R4152: sout[71] - sout[79] + sout[87] + sout[111] + sout[119]
   + xorout[111] >= 0
 R4153: - sout[71] + sout[79] + sout[87] + sout[111] + sout[119]
   + xorout[111] >= 0
 R4154: sout[71] + sout[79] + sout[87] + sout[111] + sout[119]
   - xorout[111] >= 0
 R4155: sout[71] + sout[79] - sout[87] - sout[111] - sout[119]
   + xorout[111] >= -2
 R4156: sout[71] - sout[79] + sout[87] - sout[111] - sout[119]
   + xorout[111] >= -2
 R4157: - sout[71] + sout[79] + sout[87] - sout[111] - sout[119]
   + xorout[111] >= -2
 R4158: sout[71] - sout[79] - sout[87] + sout[111] - sout[119]
   + xorout[111] >= -2
 R4159: - sout[71] + sout[79] - sout[87] + sout[111] - sout[119]
   + xorout[111] >= -2
 R4160: - sout[71] - sout[79] + sout[87] + sout[111] - sout[119]
   + xorout[111] >= -2
 R4161: sout[71] - sout[79] - sout[87] - sout[111] + sout[119]
   + xorout[111] >= -2
 R4162: - sout[71] + sout[79] - sout[87] - sout[111] + sout[119]
   + xorout[111] >= -2
 R4163: - sout[71] - sout[79] + sout[87] - sout[111] + sout[119]
   + xorout[111] >= -2
 R4164: - sout[71] - sout[79] - sout[87] + sout[111] + sout[119]
   + xorout[111] >= -2
 R4165: - sout[71] - sout[79] - sout[87] - sout[111] - sout[119]
   + xorout[111] >= -4
 R4166: - sout[118] + xorout[110] >= 0
 R4167: - sout[110] + xorout[110] >= 0
 R4168: - sout[86] + xorout[110] >= 0
 R4169: - sout[78] + xorout[110] >= 0
 R4170: - sout[70] + xorout[110] >= 0
 R4171: - sout[70] + sout[71] >= 0
 R4172: - sout[78] + sout[79] >= 0
 R4173: - sout[86] + sout[87] >= 0
 R4174: - sout[110] + sout[111] >= 0
 R4175: - sout[118] + sout[119] >= 0
 R4176: sout[70] + sout[78] + sout[86] + sout[110] + sout[118]
   - xorout[110] >= 0
 R4177: - xorout[110] + xorout[111] >= 0
 R4178: sout[71] - sout[79] - sout[87] - sout[111] - sout[119]
   + xorout[110] - xorout[111] >= -4
 R4179: - sout[71] + sout[79] - sout[87] - sout[111] - sout[119]
   + xorout[110] - xorout[111] >= -4
 R4180: - sout[71] - sout[79] + sout[87] - sout[111] - sout[119]
   + xorout[110] - xorout[111] >= -4
 R4181: - sout[71] - sout[79] - sout[87] + sout[111] - sout[119]
   + xorout[110] - xorout[111] >= -4
 R4182: - sout[71] - sout[79] - sout[87] - sout[111] + sout[119]
   + xorout[110] - xorout[111] >= -4
 R4183: sout[71] + sout[79] + sout[87] - sout[111] - sout[119]
   + xorout[110] - xorout[111] >= -2
 R4184: sout[71] + sout[79] - sout[87] + sout[111] - sout[119]
   + xorout[110] - xorout[111] >= -2
 R4185: sout[71] - sout[79] + sout[87] + sout[111] - sout[119]
   + xorout[110] - xorout[111] >= -2
 R4186: - sout[71] + sout[79] + sout[87] + sout[111] - sout[119]
   + xorout[110] - xorout[111] >= -2
 R4187: sout[71] + sout[79] - sout[87] - sout[111] + sout[119]
   + xorout[110] - xorout[111] >= -2
 R4188: sout[71] - sout[79] + sout[87] - sout[111] + sout[119]
   + xorout[110] - xorout[111] >= -2
 R4189: - sout[71] + sout[79] + sout[87] - sout[111] + sout[119]
   + xorout[110] - xorout[111] >= -2
 R4190: sout[71] - sout[79] - sout[87] + sout[111] + sout[119]
   + xorout[110] - xorout[111] >= -2
 R4191: - sout[71] + sout[79] - sout[87] + sout[111] + sout[119]
   + xorout[110] - xorout[111] >= -2
 R4192: - sout[71] - sout[79] + sout[87] + sout[111] + sout[119]
   + xorout[110] - xorout[111] >= -2
 R4193: sout[73] + sout[81] + sout[89] + sout[113] - sout[121]
   + xorout[113] >= 0
 R4194: sout[73] + sout[81] + sout[89] - sout[113] + sout[121]
   + xorout[113] >= 0
 R4195: sout[73] + sout[81] - sout[89] + sout[113] + sout[121]
   + xorout[113] >= 0
 R4196: sout[73] - sout[81] + sout[89] + sout[113] + sout[121]
   + xorout[113] >= 0
 R4197: - sout[73] + sout[81] + sout[89] + sout[113] + sout[121]
   + xorout[113] >= 0
 R4198: sout[73] + sout[81] + sout[89] + sout[113] + sout[121]
   - xorout[113] >= 0
 R4199: sout[73] + sout[81] - sout[89] - sout[113] - sout[121]
   + xorout[113] >= -2
 R4200: sout[73] - sout[81] + sout[89] - sout[113] - sout[121]
   + xorout[113] >= -2
 R4201: - sout[73] + sout[81] + sout[89] - sout[113] - sout[121]
   + xorout[113] >= -2
 R4202: sout[73] - sout[81] - sout[89] + sout[113] - sout[121]
   + xorout[113] >= -2
 R4203: - sout[73] + sout[81] - sout[89] + sout[113] - sout[121]
   + xorout[113] >= -2
 R4204: - sout[73] - sout[81] + sout[89] + sout[113] - sout[121]
   + xorout[113] >= -2
 R4205: sout[73] - sout[81] - sout[89] - sout[113] + sout[121]
   + xorout[113] >= -2
 R4206: - sout[73] + sout[81] - sout[89] - sout[113] + sout[121]
   + xorout[113] >= -2
 R4207: - sout[73] - sout[81] + sout[89] - sout[113] + sout[121]
   + xorout[113] >= -2
 R4208: - sout[73] - sout[81] - sout[89] + sout[113] + sout[121]
   + xorout[113] >= -2
 R4209: - sout[73] - sout[81] - sout[89] - sout[113] - sout[121]
   + xorout[113] >= -4
 R4210: - sout[120] + xorout[112] >= 0
 R4211: - sout[112] + xorout[112] >= 0
 R4212: - sout[88] + xorout[112] >= 0
 R4213: - sout[80] + xorout[112] >= 0
 R4214: - sout[72] + xorout[112] >= 0
 R4215: - sout[72] + sout[73] >= 0
 R4216: - sout[80] + sout[81] >= 0
 R4217: - sout[88] + sout[89] >= 0
 R4218: - sout[112] + sout[113] >= 0
 R4219: - sout[120] + sout[121] >= 0
 R4220: sout[72] + sout[80] + sout[88] + sout[112] + sout[120]
   - xorout[112] >= 0
 R4221: - xorout[112] + xorout[113] >= 0
 R4222: sout[73] - sout[81] - sout[89] - sout[113] - sout[121]
   + xorout[112] - xorout[113] >= -4
 R4223: - sout[73] + sout[81] - sout[89] - sout[113] - sout[121]
   + xorout[112] - xorout[113] >= -4
 R4224: - sout[73] - sout[81] + sout[89] - sout[113] - sout[121]
   + xorout[112] - xorout[113] >= -4
 R4225: - sout[73] - sout[81] - sout[89] + sout[113] - sout[121]
   + xorout[112] - xorout[113] >= -4
 R4226: - sout[73] - sout[81] - sout[89] - sout[113] + sout[121]
   + xorout[112] - xorout[113] >= -4
 R4227: sout[73] + sout[81] + sout[89] - sout[113] - sout[121]
   + xorout[112] - xorout[113] >= -2
 R4228: sout[73] + sout[81] - sout[89] + sout[113] - sout[121]
   + xorout[112] - xorout[113] >= -2
 R4229: sout[73] - sout[81] + sout[89] + sout[113] - sout[121]
   + xorout[112] - xorout[113] >= -2
 R4230: - sout[73] + sout[81] + sout[89] + sout[113] - sout[121]
   + xorout[112] - xorout[113] >= -2
 R4231: sout[73] + sout[81] - sout[89] - sout[113] + sout[121]
   + xorout[112] - xorout[113] >= -2
 R4232: sout[73] - sout[81] + sout[89] - sout[113] + sout[121]
   + xorout[112] - xorout[113] >= -2
 R4233: - sout[73] + sout[81] + sout[89] - sout[113] + sout[121]
   + xorout[112] - xorout[113] >= -2
 R4234: sout[73] - sout[81] - sout[89] + sout[113] + sout[121]
   + xorout[112] - xorout[113] >= -2
 R4235: - sout[73] + sout[81] - sout[89] + sout[113] + sout[121]
   + xorout[112] - xorout[113] >= -2
 R4236: - sout[73] - sout[81] + sout[89] + sout[113] + sout[121]
   + xorout[112] - xorout[113] >= -2
 R4237: sout[75] + sout[83] + sout[91] + sout[115] - sout[123]
   + xorout[115] >= 0
 R4238: sout[75] + sout[83] + sout[91] - sout[115] + sout[123]
   + xorout[115] >= 0
 R4239: sout[75] + sout[83] - sout[91] + sout[115] + sout[123]
   + xorout[115] >= 0
 R4240: sout[75] - sout[83] + sout[91] + sout[115] + sout[123]
   + xorout[115] >= 0
 R4241: - sout[75] + sout[83] + sout[91] + sout[115] + sout[123]
   + xorout[115] >= 0
 R4242: sout[75] + sout[83] + sout[91] + sout[115] + sout[123]
   - xorout[115] >= 0
 R4243: sout[75] + sout[83] - sout[91] - sout[115] - sout[123]
   + xorout[115] >= -2
 R4244: sout[75] - sout[83] + sout[91] - sout[115] - sout[123]
   + xorout[115] >= -2
 R4245: - sout[75] + sout[83] + sout[91] - sout[115] - sout[123]
   + xorout[115] >= -2
 R4246: sout[75] - sout[83] - sout[91] + sout[115] - sout[123]
   + xorout[115] >= -2
 R4247: - sout[75] + sout[83] - sout[91] + sout[115] - sout[123]
   + xorout[115] >= -2
 R4248: - sout[75] - sout[83] + sout[91] + sout[115] - sout[123]
   + xorout[115] >= -2
 R4249: sout[75] - sout[83] - sout[91] - sout[115] + sout[123]
   + xorout[115] >= -2
 R4250: - sout[75] + sout[83] - sout[91] - sout[115] + sout[123]
   + xorout[115] >= -2
 R4251: - sout[75] - sout[83] + sout[91] - sout[115] + sout[123]
   + xorout[115] >= -2
 R4252: - sout[75] - sout[83] - sout[91] + sout[115] + sout[123]
   + xorout[115] >= -2
 R4253: - sout[75] - sout[83] - sout[91] - sout[115] - sout[123]
   + xorout[115] >= -4
 R4254: - sout[122] + xorout[114] >= 0
 R4255: - sout[114] + xorout[114] >= 0
 R4256: - sout[90] + xorout[114] >= 0
 R4257: - sout[82] + xorout[114] >= 0
 R4258: - sout[74] + xorout[114] >= 0
 R4259: - sout[74] + sout[75] >= 0
 R4260: - sout[82] + sout[83] >= 0
 R4261: - sout[90] + sout[91] >= 0
 R4262: - sout[114] + sout[115] >= 0
 R4263: - sout[122] + sout[123] >= 0
 R4264: sout[74] + sout[82] + sout[90] + sout[114] + sout[122]
   - xorout[114] >= 0
 R4265: - xorout[114] + xorout[115] >= 0
 R4266: sout[75] - sout[83] - sout[91] - sout[115] - sout[123]
   + xorout[114] - xorout[115] >= -4
 R4267: - sout[75] + sout[83] - sout[91] - sout[115] - sout[123]
   + xorout[114] - xorout[115] >= -4
 R4268: - sout[75] - sout[83] + sout[91] - sout[115] - sout[123]
   + xorout[114] - xorout[115] >= -4
 R4269: - sout[75] - sout[83] - sout[91] + sout[115] - sout[123]
   + xorout[114] - xorout[115] >= -4
 R4270: - sout[75] - sout[83] - sout[91] - sout[115] + sout[123]
   + xorout[114] - xorout[115] >= -4
 R4271: sout[75] + sout[83] + sout[91] - sout[115] - sout[123]
   + xorout[114] - xorout[115] >= -2
 R4272: sout[75] + sout[83] - sout[91] + sout[115] - sout[123]
   + xorout[114] - xorout[115] >= -2
 R4273: sout[75] - sout[83] + sout[91] + sout[115] - sout[123]
   + xorout[114] - xorout[115] >= -2
 R4274: - sout[75] + sout[83] + sout[91] + sout[115] - sout[123]
   + xorout[114] - xorout[115] >= -2
 R4275: sout[75] + sout[83] - sout[91] - sout[115] + sout[123]
   + xorout[114] - xorout[115] >= -2
 R4276: sout[75] - sout[83] + sout[91] - sout[115] + sout[123]
   + xorout[114] - xorout[115] >= -2
 R4277: - sout[75] + sout[83] + sout[91] - sout[115] + sout[123]
   + xorout[114] - xorout[115] >= -2
 R4278: sout[75] - sout[83] - sout[91] + sout[115] + sout[123]
   + xorout[114] - xorout[115] >= -2
 R4279: - sout[75] + sout[83] - sout[91] + sout[115] + sout[123]
   + xorout[114] - xorout[115] >= -2
 R4280: - sout[75] - sout[83] + sout[91] + sout[115] + sout[123]
   + xorout[114] - xorout[115] >= -2
 R4281: sout[77] + sout[85] + sout[93] + sout[117] - sout[125]
   + xorout[117] >= 0
 R4282: sout[77] + sout[85] + sout[93] - sout[117] + sout[125]
   + xorout[117] >= 0
 R4283: sout[77] + sout[85] - sout[93] + sout[117] + sout[125]
   + xorout[117] >= 0
 R4284: sout[77] - sout[85] + sout[93] + sout[117] + sout[125]
   + xorout[117] >= 0
 R4285: - sout[77] + sout[85] + sout[93] + sout[117] + sout[125]
   + xorout[117] >= 0
 R4286: sout[77] + sout[85] + sout[93] + sout[117] + sout[125]
   - xorout[117] >= 0
 R4287: sout[77] + sout[85] - sout[93] - sout[117] - sout[125]
   + xorout[117] >= -2
 R4288: sout[77] - sout[85] + sout[93] - sout[117] - sout[125]
   + xorout[117] >= -2
 R4289: - sout[77] + sout[85] + sout[93] - sout[117] - sout[125]
   + xorout[117] >= -2
 R4290: sout[77] - sout[85] - sout[93] + sout[117] - sout[125]
   + xorout[117] >= -2
 R4291: - sout[77] + sout[85] - sout[93] + sout[117] - sout[125]
   + xorout[117] >= -2
 R4292: - sout[77] - sout[85] + sout[93] + sout[117] - sout[125]
   + xorout[117] >= -2
 R4293: sout[77] - sout[85] - sout[93] - sout[117] + sout[125]
   + xorout[117] >= -2
 R4294: - sout[77] + sout[85] - sout[93] - sout[117] + sout[125]
   + xorout[117] >= -2
 R4295: - sout[77] - sout[85] + sout[93] - sout[117] + sout[125]
   + xorout[117] >= -2
 R4296: - sout[77] - sout[85] - sout[93] + sout[117] + sout[125]
   + xorout[117] >= -2
 R4297: - sout[77] - sout[85] - sout[93] - sout[117] - sout[125]
   + xorout[117] >= -4
 R4298: - sout[124] + xorout[116] >= 0
 R4299: - sout[116] + xorout[116] >= 0
 R4300: - sout[92] + xorout[116] >= 0
 R4301: - sout[84] + xorout[116] >= 0
 R4302: - sout[76] + xorout[116] >= 0
 R4303: - sout[76] + sout[77] >= 0
 R4304: - sout[84] + sout[85] >= 0
 R4305: - sout[92] + sout[93] >= 0
 R4306: - sout[116] + sout[117] >= 0
 R4307: - sout[124] + sout[125] >= 0
 R4308: sout[76] + sout[84] + sout[92] + sout[116] + sout[124]
   - xorout[116] >= 0
 R4309: - xorout[116] + xorout[117] >= 0
 R4310: sout[77] - sout[85] - sout[93] - sout[117] - sout[125]
   + xorout[116] - xorout[117] >= -4
 R4311: - sout[77] + sout[85] - sout[93] - sout[117] - sout[125]
   + xorout[116] - xorout[117] >= -4
 R4312: - sout[77] - sout[85] + sout[93] - sout[117] - sout[125]
   + xorout[116] - xorout[117] >= -4
 R4313: - sout[77] - sout[85] - sout[93] + sout[117] - sout[125]
   + xorout[116] - xorout[117] >= -4
 R4314: - sout[77] - sout[85] - sout[93] - sout[117] + sout[125]
   + xorout[116] - xorout[117] >= -4
 R4315: sout[77] + sout[85] + sout[93] - sout[117] - sout[125]
   + xorout[116] - xorout[117] >= -2
 R4316: sout[77] + sout[85] - sout[93] + sout[117] - sout[125]
   + xorout[116] - xorout[117] >= -2
 R4317: sout[77] - sout[85] + sout[93] + sout[117] - sout[125]
   + xorout[116] - xorout[117] >= -2
 R4318: - sout[77] + sout[85] + sout[93] + sout[117] - sout[125]
   + xorout[116] - xorout[117] >= -2
 R4319: sout[77] + sout[85] - sout[93] - sout[117] + sout[125]
   + xorout[116] - xorout[117] >= -2
 R4320: sout[77] - sout[85] + sout[93] - sout[117] + sout[125]
   + xorout[116] - xorout[117] >= -2
 R4321: - sout[77] + sout[85] + sout[93] - sout[117] + sout[125]
   + xorout[116] - xorout[117] >= -2
 R4322: sout[77] - sout[85] - sout[93] + sout[117] + sout[125]
   + xorout[116] - xorout[117] >= -2
 R4323: - sout[77] + sout[85] - sout[93] + sout[117] + sout[125]
   + xorout[116] - xorout[117] >= -2
 R4324: - sout[77] - sout[85] + sout[93] + sout[117] + sout[125]
   + xorout[116] - xorout[117] >= -2
 R4325: sout[79] + sout[87] + sout[95] + sout[119] - sout[127]
   + xorout[119] >= 0
 R4326: sout[79] + sout[87] + sout[95] - sout[119] + sout[127]
   + xorout[119] >= 0
 R4327: sout[79] + sout[87] - sout[95] + sout[119] + sout[127]
   + xorout[119] >= 0
 R4328: sout[79] - sout[87] + sout[95] + sout[119] + sout[127]
   + xorout[119] >= 0
 R4329: - sout[79] + sout[87] + sout[95] + sout[119] + sout[127]
   + xorout[119] >= 0
 R4330: sout[79] + sout[87] + sout[95] + sout[119] + sout[127]
   - xorout[119] >= 0
 R4331: sout[79] + sout[87] - sout[95] - sout[119] - sout[127]
   + xorout[119] >= -2
 R4332: sout[79] - sout[87] + sout[95] - sout[119] - sout[127]
   + xorout[119] >= -2
 R4333: - sout[79] + sout[87] + sout[95] - sout[119] - sout[127]
   + xorout[119] >= -2
 R4334: sout[79] - sout[87] - sout[95] + sout[119] - sout[127]
   + xorout[119] >= -2
 R4335: - sout[79] + sout[87] - sout[95] + sout[119] - sout[127]
   + xorout[119] >= -2
 R4336: - sout[79] - sout[87] + sout[95] + sout[119] - sout[127]
   + xorout[119] >= -2
 R4337: sout[79] - sout[87] - sout[95] - sout[119] + sout[127]
   + xorout[119] >= -2
 R4338: - sout[79] + sout[87] - sout[95] - sout[119] + sout[127]
   + xorout[119] >= -2
 R4339: - sout[79] - sout[87] + sout[95] - sout[119] + sout[127]
   + xorout[119] >= -2
 R4340: - sout[79] - sout[87] - sout[95] + sout[119] + sout[127]
   + xorout[119] >= -2
 R4341: - sout[79] - sout[87] - sout[95] - sout[119] - sout[127]
   + xorout[119] >= -4
 R4342: - sout[126] + xorout[118] >= 0
 R4343: - sout[118] + xorout[118] >= 0
 R4344: - sout[94] + xorout[118] >= 0
 R4345: - sout[86] + xorout[118] >= 0
 R4346: - sout[78] + xorout[118] >= 0
 R4347: - sout[78] + sout[79] >= 0
 R4348: - sout[86] + sout[87] >= 0
 R4349: - sout[94] + sout[95] >= 0
 R4350: - sout[118] + sout[119] >= 0
 R4351: - sout[126] + sout[127] >= 0
 R4352: sout[78] + sout[86] + sout[94] + sout[118] + sout[126]
   - xorout[118] >= 0
 R4353: - xorout[118] + xorout[119] >= 0
 R4354: sout[79] - sout[87] - sout[95] - sout[119] - sout[127]
   + xorout[118] - xorout[119] >= -4
 R4355: - sout[79] + sout[87] - sout[95] - sout[119] - sout[127]
   + xorout[118] - xorout[119] >= -4
 R4356: - sout[79] - sout[87] + sout[95] - sout[119] - sout[127]
   + xorout[118] - xorout[119] >= -4
 R4357: - sout[79] - sout[87] - sout[95] + sout[119] - sout[127]
   + xorout[118] - xorout[119] >= -4
 R4358: - sout[79] - sout[87] - sout[95] - sout[119] + sout[127]
   + xorout[118] - xorout[119] >= -4
 R4359: sout[79] + sout[87] + sout[95] - sout[119] - sout[127]
   + xorout[118] - xorout[119] >= -2
 R4360: sout[79] + sout[87] - sout[95] + sout[119] - sout[127]
   + xorout[118] - xorout[119] >= -2
 R4361: sout[79] - sout[87] + sout[95] + sout[119] - sout[127]
   + xorout[118] - xorout[119] >= -2
 R4362: - sout[79] + sout[87] + sout[95] + sout[119] - sout[127]
   + xorout[118] - xorout[119] >= -2
 R4363: sout[79] + sout[87] - sout[95] - sout[119] + sout[127]
   + xorout[118] - xorout[119] >= -2
 R4364: sout[79] - sout[87] + sout[95] - sout[119] + sout[127]
   + xorout[118] - xorout[119] >= -2
 R4365: - sout[79] + sout[87] + sout[95] - sout[119] + sout[127]
   + xorout[118] - xorout[119] >= -2
 R4366: sout[79] - sout[87] - sout[95] + sout[119] + sout[127]
   + xorout[118] - xorout[119] >= -2
 R4367: - sout[79] + sout[87] - sout[95] + sout[119] + sout[127]
   + xorout[118] - xorout[119] >= -2
 R4368: - sout[79] - sout[87] + sout[95] + sout[119] + sout[127]
   + xorout[118] - xorout[119] >= -2
 R4369: sout[65] + sout[81] + sout[89] + sout[97] - sout[121] + xorout[121]
   >= 0
 R4370: sout[65] + sout[81] + sout[89] - sout[97] + sout[121] + xorout[121]
   >= 0
 R4371: sout[65] + sout[81] - sout[89] + sout[97] + sout[121] + xorout[121]
   >= 0
 R4372: sout[65] - sout[81] + sout[89] + sout[97] + sout[121] + xorout[121]
   >= 0
 R4373: - sout[65] + sout[81] + sout[89] + sout[97] + sout[121]
   + xorout[121] >= 0
 R4374: sout[65] + sout[81] + sout[89] + sout[97] + sout[121] - xorout[121]
   >= 0
 R4375: sout[65] + sout[81] - sout[89] - sout[97] - sout[121] + xorout[121]
   >= -2
 R4376: sout[65] - sout[81] + sout[89] - sout[97] - sout[121] + xorout[121]
   >= -2
 R4377: - sout[65] + sout[81] + sout[89] - sout[97] - sout[121]
   + xorout[121] >= -2
 R4378: sout[65] - sout[81] - sout[89] + sout[97] - sout[121] + xorout[121]
   >= -2
 R4379: - sout[65] + sout[81] - sout[89] + sout[97] - sout[121]
   + xorout[121] >= -2
 R4380: - sout[65] - sout[81] + sout[89] + sout[97] - sout[121]
   + xorout[121] >= -2
 R4381: sout[65] - sout[81] - sout[89] - sout[97] + sout[121] + xorout[121]
   >= -2
 R4382: - sout[65] + sout[81] - sout[89] - sout[97] + sout[121]
   + xorout[121] >= -2
 R4383: - sout[65] - sout[81] + sout[89] - sout[97] + sout[121]
   + xorout[121] >= -2
 R4384: - sout[65] - sout[81] - sout[89] + sout[97] + sout[121]
   + xorout[121] >= -2
 R4385: - sout[65] - sout[81] - sout[89] - sout[97] - sout[121]
   + xorout[121] >= -4
 R4386: - sout[120] + xorout[120] >= 0
 R4387: - sout[96] + xorout[120] >= 0
 R4388: - sout[88] + xorout[120] >= 0
 R4389: - sout[80] + xorout[120] >= 0
 R4390: - sout[64] + xorout[120] >= 0
 R4391: - sout[64] + sout[65] >= 0
 R4392: - sout[80] + sout[81] >= 0
 R4393: - sout[88] + sout[89] >= 0
 R4394: - sout[96] + sout[97] >= 0
 R4395: - sout[120] + sout[121] >= 0
 R4396: sout[64] + sout[80] + sout[88] + sout[96] + sout[120] - xorout[120]
   >= 0
 R4397: - xorout[120] + xorout[121] >= 0
 R4398: sout[65] - sout[81] - sout[89] - sout[97] - sout[121] + xorout[120]
   - xorout[121] >= -4
 R4399: - sout[65] + sout[81] - sout[89] - sout[97] - sout[121]
   + xorout[120] - xorout[121] >= -4
 R4400: - sout[65] - sout[81] + sout[89] - sout[97] - sout[121]
   + xorout[120] - xorout[121] >= -4
 R4401: - sout[65] - sout[81] - sout[89] + sout[97] - sout[121]
   + xorout[120] - xorout[121] >= -4
 R4402: - sout[65] - sout[81] - sout[89] - sout[97] + sout[121]
   + xorout[120] - xorout[121] >= -4
 R4403: sout[65] + sout[81] + sout[89] - sout[97] - sout[121] + xorout[120]
   - xorout[121] >= -2
 R4404: sout[65] + sout[81] - sout[89] + sout[97] - sout[121] + xorout[120]
   - xorout[121] >= -2
 R4405: sout[65] - sout[81] + sout[89] + sout[97] - sout[121] + xorout[120]
   - xorout[121] >= -2
 R4406: - sout[65] + sout[81] + sout[89] + sout[97] - sout[121]
   + xorout[120] - xorout[121] >= -2
 R4407: sout[65] + sout[81] - sout[89] - sout[97] + sout[121] + xorout[120]
   - xorout[121] >= -2
 R4408: sout[65] - sout[81] + sout[89] - sout[97] + sout[121] + xorout[120]
   - xorout[121] >= -2
 R4409: - sout[65] + sout[81] + sout[89] - sout[97] + sout[121]
   + xorout[120] - xorout[121] >= -2
 R4410: sout[65] - sout[81] - sout[89] + sout[97] + sout[121] + xorout[120]
   - xorout[121] >= -2
 R4411: - sout[65] + sout[81] - sout[89] + sout[97] + sout[121]
   + xorout[120] - xorout[121] >= -2
 R4412: - sout[65] - sout[81] + sout[89] + sout[97] + sout[121]
   + xorout[120] - xorout[121] >= -2
 R4413: sout[67] + sout[83] + sout[91] + sout[99] - sout[123] + xorout[123]
   >= 0
 R4414: sout[67] + sout[83] + sout[91] - sout[99] + sout[123] + xorout[123]
   >= 0
 R4415: sout[67] + sout[83] - sout[91] + sout[99] + sout[123] + xorout[123]
   >= 0
 R4416: sout[67] - sout[83] + sout[91] + sout[99] + sout[123] + xorout[123]
   >= 0
 R4417: - sout[67] + sout[83] + sout[91] + sout[99] + sout[123]
   + xorout[123] >= 0
 R4418: sout[67] + sout[83] + sout[91] + sout[99] + sout[123] - xorout[123]
   >= 0
 R4419: sout[67] + sout[83] - sout[91] - sout[99] - sout[123] + xorout[123]
   >= -2
 R4420: sout[67] - sout[83] + sout[91] - sout[99] - sout[123] + xorout[123]
   >= -2
 R4421: - sout[67] + sout[83] + sout[91] - sout[99] - sout[123]
   + xorout[123] >= -2
 R4422: sout[67] - sout[83] - sout[91] + sout[99] - sout[123] + xorout[123]
   >= -2
 R4423: - sout[67] + sout[83] - sout[91] + sout[99] - sout[123]
   + xorout[123] >= -2
 R4424: - sout[67] - sout[83] + sout[91] + sout[99] - sout[123]
   + xorout[123] >= -2
 R4425: sout[67] - sout[83] - sout[91] - sout[99] + sout[123] + xorout[123]
   >= -2
 R4426: - sout[67] + sout[83] - sout[91] - sout[99] + sout[123]
   + xorout[123] >= -2
 R4427: - sout[67] - sout[83] + sout[91] - sout[99] + sout[123]
   + xorout[123] >= -2
 R4428: - sout[67] - sout[83] - sout[91] + sout[99] + sout[123]
   + xorout[123] >= -2
 R4429: - sout[67] - sout[83] - sout[91] - sout[99] - sout[123]
   + xorout[123] >= -4
 R4430: - sout[122] + xorout[122] >= 0
 R4431: - sout[98] + xorout[122] >= 0
 R4432: - sout[90] + xorout[122] >= 0
 R4433: - sout[82] + xorout[122] >= 0
 R4434: - sout[66] + xorout[122] >= 0
 R4435: - sout[66] + sout[67] >= 0
 R4436: - sout[82] + sout[83] >= 0
 R4437: - sout[90] + sout[91] >= 0
 R4438: - sout[98] + sout[99] >= 0
 R4439: - sout[122] + sout[123] >= 0
 R4440: sout[66] + sout[82] + sout[90] + sout[98] + sout[122] - xorout[122]
   >= 0
 R4441: - xorout[122] + xorout[123] >= 0
 R4442: sout[67] - sout[83] - sout[91] - sout[99] - sout[123] + xorout[122]
   - xorout[123] >= -4
 R4443: - sout[67] + sout[83] - sout[91] - sout[99] - sout[123]
   + xorout[122] - xorout[123] >= -4
 R4444: - sout[67] - sout[83] + sout[91] - sout[99] - sout[123]
   + xorout[122] - xorout[123] >= -4
 R4445: - sout[67] - sout[83] - sout[91] + sout[99] - sout[123]
   + xorout[122] - xorout[123] >= -4
 R4446: - sout[67] - sout[83] - sout[91] - sout[99] + sout[123]
   + xorout[122] - xorout[123] >= -4
 R4447: sout[67] + sout[83] + sout[91] - sout[99] - sout[123] + xorout[122]
   - xorout[123] >= -2
 R4448: sout[67] + sout[83] - sout[91] + sout[99] - sout[123] + xorout[122]
   - xorout[123] >= -2
 R4449: sout[67] - sout[83] + sout[91] + sout[99] - sout[123] + xorout[122]
   - xorout[123] >= -2
 R4450: - sout[67] + sout[83] + sout[91] + sout[99] - sout[123]
   + xorout[122] - xorout[123] >= -2
 R4451: sout[67] + sout[83] - sout[91] - sout[99] + sout[123] + xorout[122]
   - xorout[123] >= -2
 R4452: sout[67] - sout[83] + sout[91] - sout[99] + sout[123] + xorout[122]
   - xorout[123] >= -2
 R4453: - sout[67] + sout[83] + sout[91] - sout[99] + sout[123]
   + xorout[122] - xorout[123] >= -2
 R4454: sout[67] - sout[83] - sout[91] + sout[99] + sout[123] + xorout[122]
   - xorout[123] >= -2
 R4455: - sout[67] + sout[83] - sout[91] + sout[99] + sout[123]
   + xorout[122] - xorout[123] >= -2
 R4456: - sout[67] - sout[83] + sout[91] + sout[99] + sout[123]
   + xorout[122] - xorout[123] >= -2
 R4457: sout[69] + sout[85] + sout[93] + sout[101] - sout[125]
   + xorout[125] >= 0
 R4458: sout[69] + sout[85] + sout[93] - sout[101] + sout[125]
   + xorout[125] >= 0
 R4459: sout[69] + sout[85] - sout[93] + sout[101] + sout[125]
   + xorout[125] >= 0
 R4460: sout[69] - sout[85] + sout[93] + sout[101] + sout[125]
   + xorout[125] >= 0
 R4461: - sout[69] + sout[85] + sout[93] + sout[101] + sout[125]
   + xorout[125] >= 0
 R4462: sout[69] + sout[85] + sout[93] + sout[101] + sout[125]
   - xorout[125] >= 0
 R4463: sout[69] + sout[85] - sout[93] - sout[101] - sout[125]
   + xorout[125] >= -2
 R4464: sout[69] - sout[85] + sout[93] - sout[101] - sout[125]
   + xorout[125] >= -2
 R4465: - sout[69] + sout[85] + sout[93] - sout[101] - sout[125]
   + xorout[125] >= -2
 R4466: sout[69] - sout[85] - sout[93] + sout[101] - sout[125]
   + xorout[125] >= -2
 R4467: - sout[69] + sout[85] - sout[93] + sout[101] - sout[125]
   + xorout[125] >= -2
 R4468: - sout[69] - sout[85] + sout[93] + sout[101] - sout[125]
   + xorout[125] >= -2
 R4469: sout[69] - sout[85] - sout[93] - sout[101] + sout[125]
   + xorout[125] >= -2
 R4470: - sout[69] + sout[85] - sout[93] - sout[101] + sout[125]
   + xorout[125] >= -2
 R4471: - sout[69] - sout[85] + sout[93] - sout[101] + sout[125]
   + xorout[125] >= -2
 R4472: - sout[69] - sout[85] - sout[93] + sout[101] + sout[125]
   + xorout[125] >= -2
 R4473: - sout[69] - sout[85] - sout[93] - sout[101] - sout[125]
   + xorout[125] >= -4
 R4474: - sout[124] + xorout[124] >= 0
 R4475: - sout[100] + xorout[124] >= 0
 R4476: - sout[92] + xorout[124] >= 0
 R4477: - sout[84] + xorout[124] >= 0
 R4478: - sout[68] + xorout[124] >= 0
 R4479: - sout[68] + sout[69] >= 0
 R4480: - sout[84] + sout[85] >= 0
 R4481: - sout[92] + sout[93] >= 0
 R4482: - sout[100] + sout[101] >= 0
 R4483: - sout[124] + sout[125] >= 0
 R4484: sout[68] + sout[84] + sout[92] + sout[100] + sout[124]
   - xorout[124] >= 0
 R4485: - xorout[124] + xorout[125] >= 0
 R4486: sout[69] - sout[85] - sout[93] - sout[101] - sout[125]
   + xorout[124] - xorout[125] >= -4
 R4487: - sout[69] + sout[85] - sout[93] - sout[101] - sout[125]
   + xorout[124] - xorout[125] >= -4
 R4488: - sout[69] - sout[85] + sout[93] - sout[101] - sout[125]
   + xorout[124] - xorout[125] >= -4
 R4489: - sout[69] - sout[85] - sout[93] + sout[101] - sout[125]
   + xorout[124] - xorout[125] >= -4
 R4490: - sout[69] - sout[85] - sout[93] - sout[101] + sout[125]
   + xorout[124] - xorout[125] >= -4
 R4491: sout[69] + sout[85] + sout[93] - sout[101] - sout[125]
   + xorout[124] - xorout[125] >= -2
 R4492: sout[69] + sout[85] - sout[93] + sout[101] - sout[125]
   + xorout[124] - xorout[125] >= -2
 R4493: sout[69] - sout[85] + sout[93] + sout[101] - sout[125]
   + xorout[124] - xorout[125] >= -2
 R4494: - sout[69] + sout[85] + sout[93] + sout[101] - sout[125]
   + xorout[124] - xorout[125] >= -2
 R4495: sout[69] + sout[85] - sout[93] - sout[101] + sout[125]
   + xorout[124] - xorout[125] >= -2
 R4496: sout[69] - sout[85] + sout[93] - sout[101] + sout[125]
   + xorout[124] - xorout[125] >= -2
 R4497: - sout[69] + sout[85] + sout[93] - sout[101] + sout[125]
   + xorout[124] - xorout[125] >= -2
 R4498: sout[69] - sout[85] - sout[93] + sout[101] + sout[125]
   + xorout[124] - xorout[125] >= -2
 R4499: - sout[69] + sout[85] - sout[93] + sout[101] + sout[125]
   + xorout[124] - xorout[125] >= -2
 R4500: - sout[69] - sout[85] + sout[93] + sout[101] + sout[125]
   + xorout[124] - xorout[125] >= -2
 R4501: sout[71] + sout[87] + sout[95] + sout[103] - sout[127]
   + xorout[127] >= 0
 R4502: sout[71] + sout[87] + sout[95] - sout[103] + sout[127]
   + xorout[127] >= 0
 R4503: sout[71] + sout[87] - sout[95] + sout[103] + sout[127]
   + xorout[127] >= 0
 R4504: sout[71] - sout[87] + sout[95] + sout[103] + sout[127]
   + xorout[127] >= 0
 R4505: - sout[71] + sout[87] + sout[95] + sout[103] + sout[127]
   + xorout[127] >= 0
 R4506: sout[71] + sout[87] + sout[95] + sout[103] + sout[127]
   - xorout[127] >= 0
 R4507: sout[71] + sout[87] - sout[95] - sout[103] - sout[127]
   + xorout[127] >= -2
 R4508: sout[71] - sout[87] + sout[95] - sout[103] - sout[127]
   + xorout[127] >= -2
 R4509: - sout[71] + sout[87] + sout[95] - sout[103] - sout[127]
   + xorout[127] >= -2
 R4510: sout[71] - sout[87] - sout[95] + sout[103] - sout[127]
   + xorout[127] >= -2
 R4511: - sout[71] + sout[87] - sout[95] + sout[103] - sout[127]
   + xorout[127] >= -2
 R4512: - sout[71] - sout[87] + sout[95] + sout[103] - sout[127]
   + xorout[127] >= -2
 R4513: sout[71] - sout[87] - sout[95] - sout[103] + sout[127]
   + xorout[127] >= -2
 R4514: - sout[71] + sout[87] - sout[95] - sout[103] + sout[127]
   + xorout[127] >= -2
 R4515: - sout[71] - sout[87] + sout[95] - sout[103] + sout[127]
   + xorout[127] >= -2
 R4516: - sout[71] - sout[87] - sout[95] + sout[103] + sout[127]
   + xorout[127] >= -2
 R4517: - sout[71] - sout[87] - sout[95] - sout[103] - sout[127]
   + xorout[127] >= -4
 R4518: - sout[126] + xorout[126] >= 0
 R4519: - sout[102] + xorout[126] >= 0
 R4520: - sout[94] + xorout[126] >= 0
 R4521: - sout[86] + xorout[126] >= 0
 R4522: - sout[70] + xorout[126] >= 0
 R4523: - sout[70] + sout[71] >= 0
 R4524: - sout[86] + sout[87] >= 0
 R4525: - sout[94] + sout[95] >= 0
 R4526: - sout[102] + sout[103] >= 0
 R4527: - sout[126] + sout[127] >= 0
 R4528: sout[70] + sout[86] + sout[94] + sout[102] + sout[126]
   - xorout[126] >= 0
 R4529: - xorout[126] + xorout[127] >= 0
 R4530: sout[71] - sout[87] - sout[95] - sout[103] - sout[127]
   + xorout[126] - xorout[127] >= -4
 R4531: - sout[71] + sout[87] - sout[95] - sout[103] - sout[127]
   + xorout[126] - xorout[127] >= -4
 R4532: - sout[71] - sout[87] + sout[95] - sout[103] - sout[127]
   + xorout[126] - xorout[127] >= -4
 R4533: - sout[71] - sout[87] - sout[95] + sout[103] - sout[127]
   + xorout[126] - xorout[127] >= -4
 R4534: - sout[71] - sout[87] - sout[95] - sout[103] + sout[127]
   + xorout[126] - xorout[127] >= -4
 R4535: sout[71] + sout[87] + sout[95] - sout[103] - sout[127]
   + xorout[126] - xorout[127] >= -2
 R4536: sout[71] + sout[87] - sout[95] + sout[103] - sout[127]
   + xorout[126] - xorout[127] >= -2
 R4537: sout[71] - sout[87] + sout[95] + sout[103] - sout[127]
   + xorout[126] - xorout[127] >= -2
 R4538: - sout[71] + sout[87] + sout[95] + sout[103] - sout[127]
   + xorout[126] - xorout[127] >= -2
 R4539: sout[71] + sout[87] - sout[95] - sout[103] + sout[127]
   + xorout[126] - xorout[127] >= -2
 R4540: sout[71] - sout[87] + sout[95] - sout[103] + sout[127]
   + xorout[126] - xorout[127] >= -2
 R4541: - sout[71] + sout[87] + sout[95] - sout[103] + sout[127]
   + xorout[126] - xorout[127] >= -2
 R4542: sout[71] - sout[87] - sout[95] + sout[103] + sout[127]
   + xorout[126] - xorout[127] >= -2
 R4543: - sout[71] + sout[87] - sout[95] + sout[103] + sout[127]
   + xorout[126] - xorout[127] >= -2
 R4544: - sout[71] - sout[87] + sout[95] + sout[103] + sout[127]
   + xorout[126] - xorout[127] >= -2
 R4545: - x[193] + x[257] + xorout[65] >= 0
 R4546: x[193] + x[257] - xorout[65] >= 0
 R4547: - x[192] + x[256] >= 0
 R4548: x[256] - xorout[64] >= 0
 R4549: x[193] - x[257] + xorout[65] >= 0
 R4550: - xorout[64] + xorout[65] >= 0
 R4551: - x[192] + x[193] >= 0
 R4552: x[192] - x[256] + xorout[64] >= 0
 R4553: - x[256] + x[257] >= 0
 R4554: - x[193] + x[256] - x[257] - xorout[65] >= -2
 R4555: - x[195] + x[259] + xorout[67] >= 0
 R4556: x[195] + x[259] - xorout[67] >= 0
 R4557: - x[194] + x[258] >= 0
 R4558: x[258] - xorout[66] >= 0
 R4559: x[195] - x[259] + xorout[67] >= 0
 R4560: - xorout[66] + xorout[67] >= 0
 R4561: - x[194] + x[195] >= 0
 R4562: x[194] - x[258] + xorout[66] >= 0
 R4563: - x[258] + x[259] >= 0
 R4564: - x[195] + x[258] - x[259] - xorout[67] >= -2
 R4565: - x[197] + x[261] + xorout[69] >= 0
 R4566: x[197] + x[261] - xorout[69] >= 0
 R4567: - x[196] + x[260] >= 0
 R4568: x[260] - xorout[68] >= 0
 R4569: x[197] - x[261] + xorout[69] >= 0
 R4570: - xorout[68] + xorout[69] >= 0
 R4571: - x[196] + x[197] >= 0
 R4572: x[196] - x[260] + xorout[68] >= 0
 R4573: - x[260] + x[261] >= 0
 R4574: - x[197] + x[260] - x[261] - xorout[69] >= -2
 R4575: - x[199] + x[263] + xorout[71] >= 0
 R4576: x[199] + x[263] - xorout[71] >= 0
 R4577: - x[198] + x[262] >= 0
 R4578: x[262] - xorout[70] >= 0
 R4579: x[199] - x[263] + xorout[71] >= 0
 R4580: - xorout[70] + xorout[71] >= 0
 R4581: - x[198] + x[199] >= 0
 R4582: x[198] - x[262] + xorout[70] >= 0
 R4583: - x[262] + x[263] >= 0
 R4584: - x[199] + x[262] - x[263] - xorout[71] >= -2
 R4585: - x[201] + x[265] + xorout[73] >= 0
 R4586: x[201] + x[265] - xorout[73] >= 0
 R4587: - x[200] + x[264] >= 0
 R4588: x[264] - xorout[72] >= 0
 R4589: x[201] - x[265] + xorout[73] >= 0
 R4590: - xorout[72] + xorout[73] >= 0
 R4591: - x[200] + x[201] >= 0
 R4592: x[200] - x[264] + xorout[72] >= 0
 R4593: - x[264] + x[265] >= 0
 R4594: - x[201] + x[264] - x[265] - xorout[73] >= -2
 R4595: - x[203] + x[267] + xorout[75] >= 0
 R4596: x[203] + x[267] - xorout[75] >= 0
 R4597: - x[202] + x[266] >= 0
 R4598: x[266] - xorout[74] >= 0
 R4599: x[203] - x[267] + xorout[75] >= 0
 R4600: - xorout[74] + xorout[75] >= 0
 R4601: - x[202] + x[203] >= 0
 R4602: x[202] - x[266] + xorout[74] >= 0
 R4603: - x[266] + x[267] >= 0
 R4604: - x[203] + x[266] - x[267] - xorout[75] >= -2
 R4605: - x[205] + x[269] + xorout[77] >= 0
 R4606: x[205] + x[269] - xorout[77] >= 0
 R4607: - x[204] + x[268] >= 0
 R4608: x[268] - xorout[76] >= 0
 R4609: x[205] - x[269] + xorout[77] >= 0
 R4610: - xorout[76] + xorout[77] >= 0
 R4611: - x[204] + x[205] >= 0
 R4612: x[204] - x[268] + xorout[76] >= 0
 R4613: - x[268] + x[269] >= 0
 R4614: - x[205] + x[268] - x[269] - xorout[77] >= -2
 R4615: - x[207] + x[271] + xorout[79] >= 0
 R4616: x[207] + x[271] - xorout[79] >= 0
 R4617: - x[206] + x[270] >= 0
 R4618: x[270] - xorout[78] >= 0
 R4619: x[207] - x[271] + xorout[79] >= 0
 R4620: - xorout[78] + xorout[79] >= 0
 R4621: - x[206] + x[207] >= 0
 R4622: x[206] - x[270] + xorout[78] >= 0
 R4623: - x[270] + x[271] >= 0
 R4624: - x[207] + x[270] - x[271] - xorout[79] >= -2
 R4625: - x[209] + x[273] + xorout[81] >= 0
 R4626: x[209] + x[273] - xorout[81] >= 0
 R4627: - x[208] + x[272] >= 0
 R4628: x[272] - xorout[80] >= 0
 R4629: x[209] - x[273] + xorout[81] >= 0
 R4630: - xorout[80] + xorout[81] >= 0
 R4631: - x[208] + x[209] >= 0
 R4632: x[208] - x[272] + xorout[80] >= 0
 R4633: - x[272] + x[273] >= 0
 R4634: - x[209] + x[272] - x[273] - xorout[81] >= -2
 R4635: - x[211] + x[275] + xorout[83] >= 0
 R4636: x[211] + x[275] - xorout[83] >= 0
 R4637: - x[210] + x[274] >= 0
 R4638: x[274] - xorout[82] >= 0
 R4639: x[211] - x[275] + xorout[83] >= 0
 R4640: - xorout[82] + xorout[83] >= 0
 R4641: - x[210] + x[211] >= 0
 R4642: x[210] - x[274] + xorout[82] >= 0
 R4643: - x[274] + x[275] >= 0
 R4644: - x[211] + x[274] - x[275] - xorout[83] >= -2
 R4645: - x[213] + x[277] + xorout[85] >= 0
 R4646: x[213] + x[277] - xorout[85] >= 0
 R4647: - x[212] + x[276] >= 0
 R4648: x[276] - xorout[84] >= 0
 R4649: x[213] - x[277] + xorout[85] >= 0
 R4650: - xorout[84] + xorout[85] >= 0
 R4651: - x[212] + x[213] >= 0
 R4652: x[212] - x[276] + xorout[84] >= 0
 R4653: - x[276] + x[277] >= 0
 R4654: - x[213] + x[276] - x[277] - xorout[85] >= -2
 R4655: - x[215] + x[279] + xorout[87] >= 0
 R4656: x[215] + x[279] - xorout[87] >= 0
 R4657: - x[214] + x[278] >= 0
 R4658: x[278] - xorout[86] >= 0
 R4659: x[215] - x[279] + xorout[87] >= 0
 R4660: - xorout[86] + xorout[87] >= 0
 R4661: - x[214] + x[215] >= 0
 R4662: x[214] - x[278] + xorout[86] >= 0
 R4663: - x[278] + x[279] >= 0
 R4664: - x[215] + x[278] - x[279] - xorout[87] >= -2
 R4665: - x[217] + x[281] + xorout[89] >= 0
 R4666: x[217] + x[281] - xorout[89] >= 0
 R4667: - x[216] + x[280] >= 0
 R4668: x[280] - xorout[88] >= 0
 R4669: x[217] - x[281] + xorout[89] >= 0
 R4670: - xorout[88] + xorout[89] >= 0
 R4671: - x[216] + x[217] >= 0
 R4672: x[216] - x[280] + xorout[88] >= 0
 R4673: - x[280] + x[281] >= 0
 R4674: - x[217] + x[280] - x[281] - xorout[89] >= -2
 R4675: - x[219] + x[283] + xorout[91] >= 0
 R4676: x[219] + x[283] - xorout[91] >= 0
 R4677: - x[218] + x[282] >= 0
 R4678: x[282] - xorout[90] >= 0
 R4679: x[219] - x[283] + xorout[91] >= 0
 R4680: - xorout[90] + xorout[91] >= 0
 R4681: - x[218] + x[219] >= 0
 R4682: x[218] - x[282] + xorout[90] >= 0
 R4683: - x[282] + x[283] >= 0
 R4684: - x[219] + x[282] - x[283] - xorout[91] >= -2
 R4685: - x[221] + x[285] + xorout[93] >= 0
 R4686: x[221] + x[285] - xorout[93] >= 0
 R4687: - x[220] + x[284] >= 0
 R4688: x[284] - xorout[92] >= 0
 R4689: x[221] - x[285] + xorout[93] >= 0
 R4690: - xorout[92] + xorout[93] >= 0
 R4691: - x[220] + x[221] >= 0
 R4692: x[220] - x[284] + xorout[92] >= 0
 R4693: - x[284] + x[285] >= 0
 R4694: - x[221] + x[284] - x[285] - xorout[93] >= -2
 R4695: - x[223] + x[287] + xorout[95] >= 0
 R4696: x[223] + x[287] - xorout[95] >= 0
 R4697: - x[222] + x[286] >= 0
 R4698: x[286] - xorout[94] >= 0
 R4699: x[223] - x[287] + xorout[95] >= 0
 R4700: - xorout[94] + xorout[95] >= 0
 R4701: - x[222] + x[223] >= 0
 R4702: x[222] - x[286] + xorout[94] >= 0
 R4703: - x[286] + x[287] >= 0
 R4704: - x[223] + x[286] - x[287] - xorout[95] >= -2
 R4705: - x[225] + x[289] + xorout[97] >= 0
 R4706: x[225] + x[289] - xorout[97] >= 0
 R4707: - x[224] + x[288] >= 0
 R4708: x[288] - xorout[96] >= 0
 R4709: x[225] - x[289] + xorout[97] >= 0
 R4710: - xorout[96] + xorout[97] >= 0
 R4711: - x[224] + x[225] >= 0
 R4712: x[224] - x[288] + xorout[96] >= 0
 R4713: - x[288] + x[289] >= 0
 R4714: - x[225] + x[288] - x[289] - xorout[97] >= -2
 R4715: - x[227] + x[291] + xorout[99] >= 0
 R4716: x[227] + x[291] - xorout[99] >= 0
 R4717: - x[226] + x[290] >= 0
 R4718: x[290] - xorout[98] >= 0
 R4719: x[227] - x[291] + xorout[99] >= 0
 R4720: - xorout[98] + xorout[99] >= 0
 R4721: - x[226] + x[227] >= 0
 R4722: x[226] - x[290] + xorout[98] >= 0
 R4723: - x[290] + x[291] >= 0
 R4724: - x[227] + x[290] - x[291] - xorout[99] >= -2
 R4725: - x[229] + x[293] + xorout[101] >= 0
 R4726: x[229] + x[293] - xorout[101] >= 0
 R4727: - x[228] + x[292] >= 0
 R4728: x[292] - xorout[100] >= 0
 R4729: x[229] - x[293] + xorout[101] >= 0
 R4730: - xorout[100] + xorout[101] >= 0
 R4731: - x[228] + x[229] >= 0
 R4732: x[228] - x[292] + xorout[100] >= 0
 R4733: - x[292] + x[293] >= 0
 R4734: - x[229] + x[292] - x[293] - xorout[101] >= -2
 R4735: - x[231] + x[295] + xorout[103] >= 0
 R4736: x[231] + x[295] - xorout[103] >= 0
 R4737: - x[230] + x[294] >= 0
 R4738: x[294] - xorout[102] >= 0
 R4739: x[231] - x[295] + xorout[103] >= 0
 R4740: - xorout[102] + xorout[103] >= 0
 R4741: - x[230] + x[231] >= 0
 R4742: x[230] - x[294] + xorout[102] >= 0
 R4743: - x[294] + x[295] >= 0
 R4744: - x[231] + x[294] - x[295] - xorout[103] >= -2
 R4745: - x[233] + x[297] + xorout[105] >= 0
 R4746: x[233] + x[297] - xorout[105] >= 0
 R4747: - x[232] + x[296] >= 0
 R4748: x[296] - xorout[104] >= 0
 R4749: x[233] - x[297] + xorout[105] >= 0
 R4750: - xorout[104] + xorout[105] >= 0
 R4751: - x[232] + x[233] >= 0
 R4752: x[232] - x[296] + xorout[104] >= 0
 R4753: - x[296] + x[297] >= 0
 R4754: - x[233] + x[296] - x[297] - xorout[105] >= -2
 R4755: - x[235] + x[299] + xorout[107] >= 0
 R4756: x[235] + x[299] - xorout[107] >= 0
 R4757: - x[234] + x[298] >= 0
 R4758: x[298] - xorout[106] >= 0
 R4759: x[235] - x[299] + xorout[107] >= 0
 R4760: - xorout[106] + xorout[107] >= 0
 R4761: - x[234] + x[235] >= 0
 R4762: x[234] - x[298] + xorout[106] >= 0
 R4763: - x[298] + x[299] >= 0
 R4764: - x[235] + x[298] - x[299] - xorout[107] >= -2
 R4765: - x[237] + x[301] + xorout[109] >= 0
 R4766: x[237] + x[301] - xorout[109] >= 0
 R4767: - x[236] + x[300] >= 0
 R4768: x[300] - xorout[108] >= 0
 R4769: x[237] - x[301] + xorout[109] >= 0
 R4770: - xorout[108] + xorout[109] >= 0
 R4771: - x[236] + x[237] >= 0
 R4772: x[236] - x[300] + xorout[108] >= 0
 R4773: - x[300] + x[301] >= 0
 R4774: - x[237] + x[300] - x[301] - xorout[109] >= -2
 R4775: - x[239] + x[303] + xorout[111] >= 0
 R4776: x[239] + x[303] - xorout[111] >= 0
 R4777: - x[238] + x[302] >= 0
 R4778: x[302] - xorout[110] >= 0
 R4779: x[239] - x[303] + xorout[111] >= 0
 R4780: - xorout[110] + xorout[111] >= 0
 R4781: - x[238] + x[239] >= 0
 R4782: x[238] - x[302] + xorout[110] >= 0
 R4783: - x[302] + x[303] >= 0
 R4784: - x[239] + x[302] - x[303] - xorout[111] >= -2
 R4785: - x[241] + x[305] + xorout[113] >= 0
 R4786: x[241] + x[305] - xorout[113] >= 0
 R4787: - x[240] + x[304] >= 0
 R4788: x[304] - xorout[112] >= 0
 R4789: x[241] - x[305] + xorout[113] >= 0
 R4790: - xorout[112] + xorout[113] >= 0
 R4791: - x[240] + x[241] >= 0
 R4792: x[240] - x[304] + xorout[112] >= 0
 R4793: - x[304] + x[305] >= 0
 R4794: - x[241] + x[304] - x[305] - xorout[113] >= -2
 R4795: - x[243] + x[307] + xorout[115] >= 0
 R4796: x[243] + x[307] - xorout[115] >= 0
 R4797: - x[242] + x[306] >= 0
 R4798: x[306] - xorout[114] >= 0
 R4799: x[243] - x[307] + xorout[115] >= 0
 R4800: - xorout[114] + xorout[115] >= 0
 R4801: - x[242] + x[243] >= 0
 R4802: x[242] - x[306] + xorout[114] >= 0
 R4803: - x[306] + x[307] >= 0
 R4804: - x[243] + x[306] - x[307] - xorout[115] >= -2
 R4805: - x[245] + x[309] + xorout[117] >= 0
 R4806: x[245] + x[309] - xorout[117] >= 0
 R4807: - x[244] + x[308] >= 0
 R4808: x[308] - xorout[116] >= 0
 R4809: x[245] - x[309] + xorout[117] >= 0
 R4810: - xorout[116] + xorout[117] >= 0
 R4811: - x[244] + x[245] >= 0
 R4812: x[244] - x[308] + xorout[116] >= 0
 R4813: - x[308] + x[309] >= 0
 R4814: - x[245] + x[308] - x[309] - xorout[117] >= -2
 R4815: - x[247] + x[311] + xorout[119] >= 0
 R4816: x[247] + x[311] - xorout[119] >= 0
 R4817: - x[246] + x[310] >= 0
 R4818: x[310] - xorout[118] >= 0
 R4819: x[247] - x[311] + xorout[119] >= 0
 R4820: - xorout[118] + xorout[119] >= 0
 R4821: - x[246] + x[247] >= 0
 R4822: x[246] - x[310] + xorout[118] >= 0
 R4823: - x[310] + x[311] >= 0
 R4824: - x[247] + x[310] - x[311] - xorout[119] >= -2
 R4825: - x[249] + x[313] + xorout[121] >= 0
 R4826: x[249] + x[313] - xorout[121] >= 0
 R4827: - x[248] + x[312] >= 0
 R4828: x[312] - xorout[120] >= 0
 R4829: x[249] - x[313] + xorout[121] >= 0
 R4830: - xorout[120] + xorout[121] >= 0
 R4831: - x[248] + x[249] >= 0
 R4832: x[248] - x[312] + xorout[120] >= 0
 R4833: - x[312] + x[313] >= 0
 R4834: - x[249] + x[312] - x[313] - xorout[121] >= -2
 R4835: - x[251] + x[315] + xorout[123] >= 0
 R4836: x[251] + x[315] - xorout[123] >= 0
 R4837: - x[250] + x[314] >= 0
 R4838: x[314] - xorout[122] >= 0
 R4839: x[251] - x[315] + xorout[123] >= 0
 R4840: - xorout[122] + xorout[123] >= 0
 R4841: - x[250] + x[251] >= 0
 R4842: x[250] - x[314] + xorout[122] >= 0
 R4843: - x[314] + x[315] >= 0
 R4844: - x[251] + x[314] - x[315] - xorout[123] >= -2
 R4845: - x[253] + x[317] + xorout[125] >= 0
 R4846: x[253] + x[317] - xorout[125] >= 0
 R4847: - x[252] + x[316] >= 0
 R4848: x[316] - xorout[124] >= 0
 R4849: x[253] - x[317] + xorout[125] >= 0
 R4850: - xorout[124] + xorout[125] >= 0
 R4851: - x[252] + x[253] >= 0
 R4852: x[252] - x[316] + xorout[124] >= 0
 R4853: - x[316] + x[317] >= 0
 R4854: - x[253] + x[316] - x[317] - xorout[125] >= -2
 R4855: - x[255] + x[319] + xorout[127] >= 0
 R4856: x[255] + x[319] - xorout[127] >= 0
 R4857: - x[254] + x[318] >= 0
 R4858: x[318] - xorout[126] >= 0
 R4859: x[255] - x[319] + xorout[127] >= 0
 R4860: - xorout[126] + xorout[127] >= 0
 R4861: - x[254] + x[255] >= 0
 R4862: x[254] - x[318] + xorout[126] >= 0
 R4863: - x[318] + x[319] >= 0
 R4864: - x[255] + x[318] - x[319] - xorout[127] >= -2
 R4865: x[128] - x[320] = 0
 R4866: x[129] - x[321] = 0
 R4867: x[130] - x[322] = 0
 R4868: x[131] - x[323] = 0
 R4869: x[132] - x[324] = 0
 R4870: x[133] - x[325] = 0
 R4871: x[134] - x[326] = 0
 R4872: x[135] - x[327] = 0
 R4873: x[136] - x[328] = 0
 R4874: x[137] - x[329] = 0
 R4875: x[138] - x[330] = 0
 R4876: x[139] - x[331] = 0
 R4877: x[140] - x[332] = 0
 R4878: x[141] - x[333] = 0
 R4879: x[142] - x[334] = 0
 R4880: x[143] - x[335] = 0
 R4881: x[144] - x[336] = 0
 R4882: x[145] - x[337] = 0
 R4883: x[146] - x[338] = 0
 R4884: x[147] - x[339] = 0
 R4885: x[148] - x[340] = 0
 R4886: x[149] - x[341] = 0
 R4887: x[150] - x[342] = 0
 R4888: x[151] - x[343] = 0
 R4889: x[152] - x[344] = 0
 R4890: x[153] - x[345] = 0
 R4891: x[154] - x[346] = 0
 R4892: x[155] - x[347] = 0
 R4893: x[156] - x[348] = 0
 R4894: x[157] - x[349] = 0
 R4895: x[158] - x[350] = 0
 R4896: x[159] - x[351] = 0
 R4897: x[160] - x[352] = 0
 R4898: x[161] - x[353] = 0
 R4899: x[162] - x[354] = 0
 R4900: x[163] - x[355] = 0
 R4901: x[164] - x[356] = 0
 R4902: x[165] - x[357] = 0
 R4903: x[166] - x[358] = 0
 R4904: x[167] - x[359] = 0
 R4905: x[168] - x[360] = 0
 R4906: x[169] - x[361] = 0
 R4907: x[170] - x[362] = 0
 R4908: x[171] - x[363] = 0
 R4909: x[172] - x[364] = 0
 R4910: x[173] - x[365] = 0
 R4911: x[174] - x[366] = 0
 R4912: x[175] - x[367] = 0
 R4913: x[176] - x[368] = 0
 R4914: x[177] - x[369] = 0
 R4915: x[178] - x[370] = 0
 R4916: x[179] - x[371] = 0
 R4917: x[180] - x[372] = 0
 R4918: x[181] - x[373] = 0
 R4919: x[182] - x[374] = 0
 R4920: x[183] - x[375] = 0
 R4921: x[184] - x[376] = 0
 R4922: x[185] - x[377] = 0
 R4923: x[186] - x[378] = 0
 R4924: x[187] - x[379] = 0
 R4925: x[188] - x[380] = 0
 R4926: x[189] - x[381] = 0
 R4927: x[190] - x[382] = 0
 R4928: x[191] - x[383] = 0
 R4929: - x[256] + x[257] >= 0
 R4930: - x[258] + x[259] >= 0
 R4931: - x[260] + x[261] >= 0
 R4932: - x[262] + x[263] >= 0
 R4933: - x[263] + sout[135] >= 0
 R4934: sout[128] - sout[129] >= 0
 R4935: sout[129] - sout[130] >= 0
 R4936: sout[130] - sout[131] >= 0
 R4937: sout[131] - sout[132] >= 0
 R4938: sout[132] - sout[133] >= 0
 R4939: sout[133] - sout[134] >= 0
 R4940: sout[134] - sout[135] >= 0
 R4941: - x[261] + sout[135] >= 0
 R4942: - x[259] + sout[135] >= 0
 R4943: - x[257] + sout[135] >= 0
 R4944: x[257] + x[259] + x[261] + x[263] - sout[128] >= 0
 R4945: - x[264] + x[265] >= 0
 R4946: - x[266] + x[267] >= 0
 R4947: - x[268] + x[269] >= 0
 R4948: - x[270] + x[271] >= 0
 R4949: - x[271] + sout[143] >= 0
 R4950: sout[136] - sout[137] >= 0
 R4951: sout[137] - sout[138] >= 0
 R4952: sout[138] - sout[139] >= 0
 R4953: sout[139] - sout[140] >= 0
 R4954: sout[140] - sout[141] >= 0
 R4955: sout[141] - sout[142] >= 0
 R4956: sout[142] - sout[143] >= 0
 R4957: - x[269] + sout[143] >= 0
 R4958: - x[267] + sout[143] >= 0
 R4959: - x[265] + sout[143] >= 0
 R4960: x[265] + x[267] + x[269] + x[271] - sout[136] >= 0
 R4961: - x[272] + x[273] >= 0
 R4962: - x[274] + x[275] >= 0
 R4963: - x[276] + x[277] >= 0
 R4964: - x[278] + x[279] >= 0
 R4965: - x[279] + sout[151] >= 0
 R4966: sout[144] - sout[145] >= 0
 R4967: sout[145] - sout[146] >= 0
 R4968: sout[146] - sout[147] >= 0
 R4969: sout[147] - sout[148] >= 0
 R4970: sout[148] - sout[149] >= 0
 R4971: sout[149] - sout[150] >= 0
 R4972: sout[150] - sout[151] >= 0
 R4973: - x[277] + sout[151] >= 0
 R4974: - x[275] + sout[151] >= 0
 R4975: - x[273] + sout[151] >= 0
 R4976: x[273] + x[275] + x[277] + x[279] - sout[144] >= 0
 R4977: - x[280] + x[281] >= 0
 R4978: - x[282] + x[283] >= 0
 R4979: - x[284] + x[285] >= 0
 R4980: - x[286] + x[287] >= 0
 R4981: - x[287] + sout[159] >= 0
 R4982: sout[152] - sout[153] >= 0
 R4983: sout[153] - sout[154] >= 0
 R4984: sout[154] - sout[155] >= 0
 R4985: sout[155] - sout[156] >= 0
 R4986: sout[156] - sout[157] >= 0
 R4987: sout[157] - sout[158] >= 0
 R4988: sout[158] - sout[159] >= 0
 R4989: - x[285] + sout[159] >= 0
 R4990: - x[283] + sout[159] >= 0
 R4991: - x[281] + sout[159] >= 0
 R4992: x[281] + x[283] + x[285] + x[287] - sout[152] >= 0
 R4993: - x[288] + x[289] >= 0
 R4994: - x[290] + x[291] >= 0
 R4995: - x[292] + x[293] >= 0
 R4996: - x[294] + x[295] >= 0
 R4997: - x[295] + sout[167] >= 0
 R4998: sout[160] - sout[161] >= 0
 R4999: sout[161] - sout[162] >= 0
 R5000: sout[162] - sout[163] >= 0
 R5001: sout[163] - sout[164] >= 0
 R5002: sout[164] - sout[165] >= 0
 R5003: sout[165] - sout[166] >= 0
 R5004: sout[166] - sout[167] >= 0
 R5005: - x[293] + sout[167] >= 0
 R5006: - x[291] + sout[167] >= 0
 R5007: - x[289] + sout[167] >= 0
 R5008: x[289] + x[291] + x[293] + x[295] - sout[160] >= 0
 R5009: - x[296] + x[297] >= 0
 R5010: - x[298] + x[299] >= 0
 R5011: - x[300] + x[301] >= 0
 R5012: - x[302] + x[303] >= 0
 R5013: - x[303] + sout[175] >= 0
 R5014: sout[168] - sout[169] >= 0
 R5015: sout[169] - sout[170] >= 0
 R5016: sout[170] - sout[171] >= 0
 R5017: sout[171] - sout[172] >= 0
 R5018: sout[172] - sout[173] >= 0
 R5019: sout[173] - sout[174] >= 0
 R5020: sout[174] - sout[175] >= 0
 R5021: - x[301] + sout[175] >= 0
 R5022: - x[299] + sout[175] >= 0
 R5023: - x[297] + sout[175] >= 0
 R5024: x[297] + x[299] + x[301] + x[303] - sout[168] >= 0
 R5025: - x[304] + x[305] >= 0
 R5026: - x[306] + x[307] >= 0
 R5027: - x[308] + x[309] >= 0
 R5028: - x[310] + x[311] >= 0
 R5029: - x[311] + sout[183] >= 0
 R5030: sout[176] - sout[177] >= 0
 R5031: sout[177] - sout[178] >= 0
 R5032: sout[178] - sout[179] >= 0
 R5033: sout[179] - sout[180] >= 0
 R5034: sout[180] - sout[181] >= 0
 R5035: sout[181] - sout[182] >= 0
 R5036: sout[182] - sout[183] >= 0
 R5037: - x[309] + sout[183] >= 0
 R5038: - x[307] + sout[183] >= 0
 R5039: - x[305] + sout[183] >= 0
 R5040: x[305] + x[307] + x[309] + x[311] - sout[176] >= 0
 R5041: - x[312] + x[313] >= 0
 R5042: - x[314] + x[315] >= 0
 R5043: - x[316] + x[317] >= 0
 R5044: - x[318] + x[319] >= 0
 R5045: - x[319] + sout[191] >= 0
 R5046: sout[184] - sout[185] >= 0
 R5047: sout[185] - sout[186] >= 0
 R5048: sout[186] - sout[187] >= 0
 R5049: sout[187] - sout[188] >= 0
 R5050: sout[188] - sout[189] >= 0
 R5051: sout[189] - sout[190] >= 0
 R5052: sout[190] - sout[191] >= 0
 R5053: - x[317] + sout[191] >= 0
 R5054: - x[315] + sout[191] >= 0
 R5055: - x[313] + sout[191] >= 0
 R5056: x[313] + x[315] + x[317] + x[319] - sout[184] >= 0
 R5057: sout[137] + sout[145] + sout[153] + sout[161] + sout[169]
   - sout[177] + xorout[129] >= 0
 R5058: sout[137] + sout[145] + sout[153] + sout[161] - sout[169]
   + sout[177] + xorout[129] >= 0
 R5059: sout[137] + sout[145] + sout[153] - sout[161] + sout[169]
   + sout[177] + xorout[129] >= 0
 R5060: sout[137] + sout[145] - sout[153] + sout[161] + sout[169]
   + sout[177] + xorout[129] >= 0
 R5061: sout[137] - sout[145] + sout[153] + sout[161] + sout[169]
   + sout[177] + xorout[129] >= 0
 R5062: - sout[137] + sout[145] + sout[153] + sout[161] + sout[169]
   + sout[177] + xorout[129] >= 0
 R5063: sout[137] + sout[145] + sout[153] + sout[161] + sout[169]
   + sout[177] - xorout[129] >= 0
 R5064: sout[137] + sout[145] + sout[153] - sout[161] - sout[169]
   - sout[177] + xorout[129] >= -2
 R5065: sout[137] + sout[145] - sout[153] + sout[161] - sout[169]
   - sout[177] + xorout[129] >= -2
 R5066: sout[137] - sout[145] + sout[153] + sout[161] - sout[169]
   - sout[177] + xorout[129] >= -2
 R5067: - sout[137] + sout[145] + sout[153] + sout[161] - sout[169]
   - sout[177] + xorout[129] >= -2
 R5068: sout[137] + sout[145] - sout[153] - sout[161] + sout[169]
   - sout[177] + xorout[129] >= -2
 R5069: sout[137] - sout[145] + sout[153] - sout[161] + sout[169]
   - sout[177] + xorout[129] >= -2
 R5070: - sout[137] + sout[145] + sout[153] - sout[161] + sout[169]
   - sout[177] + xorout[129] >= -2
 R5071: sout[137] - sout[145] - sout[153] + sout[161] + sout[169]
   - sout[177] + xorout[129] >= -2
 R5072: - sout[137] + sout[145] - sout[153] + sout[161] + sout[169]
   - sout[177] + xorout[129] >= -2
 R5073: - sout[137] - sout[145] + sout[153] + sout[161] + sout[169]
   - sout[177] + xorout[129] >= -2
 R5074: sout[137] + sout[145] - sout[153] - sout[161] - sout[169]
   + sout[177] + xorout[129] >= -2
 R5075: sout[137] - sout[145] + sout[153] - sout[161] - sout[169]
   + sout[177] + xorout[129] >= -2
 R5076: - sout[137] + sout[145] + sout[153] - sout[161] - sout[169]
   + sout[177] + xorout[129] >= -2
 R5077: sout[137] - sout[145] - sout[153] + sout[161] - sout[169]
   + sout[177] + xorout[129] >= -2
 R5078: - sout[137] + sout[145] - sout[153] + sout[161] - sout[169]
   + sout[177] + xorout[129] >= -2
 R5079: - sout[137] - sout[145] + sout[153] + sout[161] - sout[169]
   + sout[177] + xorout[129] >= -2
 R5080: sout[137] - sout[145] - sout[153] - sout[161] + sout[169]
   + sout[177] + xorout[129] >= -2
 R5081: - sout[137] + sout[145] - sout[153] - sout[161] + sout[169]
   + sout[177] + xorout[129] >= -2
 R5082: - sout[137] - sout[145] + sout[153] - sout[161] + sout[169]
   + sout[177] + xorout[129] >= -2
 R5083: - sout[137] - sout[145] - sout[153] + sout[161] + sout[169]
   + sout[177] + xorout[129] >= -2
 R5084: sout[137] - sout[145] - sout[153] - sout[161] - sout[169]
   - sout[177] + xorout[129] >= -4
 R5085: - sout[137] + sout[145] - sout[153] - sout[161] - sout[169]
   - sout[177] + xorout[129] >= -4
 R5086: - sout[137] - sout[145] + sout[153] - sout[161] - sout[169]
   - sout[177] + xorout[129] >= -4
 R5087: - sout[137] - sout[145] - sout[153] + sout[161] - sout[169]
   - sout[177] + xorout[129] >= -4
 R5088: - sout[137] - sout[145] - sout[153] - sout[161] + sout[169]
   - sout[177] + xorout[129] >= -4
 R5089: - sout[137] - sout[145] - sout[153] - sout[161] - sout[169]
   + sout[177] + xorout[129] >= -4
 R5090: - sout[176] + xorout[128] >= 0
 R5091: - sout[168] + xorout[128] >= 0
 R5092: - sout[160] + xorout[128] >= 0
 R5093: - sout[152] + xorout[128] >= 0
 R5094: - sout[144] + xorout[128] >= 0
 R5095: - sout[136] + xorout[128] >= 0
 R5096: - sout[136] + sout[137] >= 0
 R5097: - sout[144] + sout[145] >= 0
 R5098: - sout[152] + sout[153] >= 0
 R5099: - sout[160] + sout[161] >= 0
 R5100: - sout[168] + sout[169] >= 0
 R5101: - sout[176] + sout[177] >= 0
 R5102: sout[136] + sout[144] + sout[152] + sout[160] + sout[168]
   + sout[176] - xorout[128] >= 0
 R5103: - xorout[128] + xorout[129] >= 0
 R5104: - sout[137] - sout[145] - sout[153] - sout[161] - sout[169]
   - sout[177] + xorout[128] - xorout[129] >= -6
 R5105: sout[137] + sout[145] - sout[153] - sout[161] - sout[169]
   - sout[177] + xorout[128] - xorout[129] >= -4
 R5106: sout[137] - sout[145] + sout[153] - sout[161] - sout[169]
   - sout[177] + xorout[128] - xorout[129] >= -4
 R5107: - sout[137] + sout[145] + sout[153] - sout[161] - sout[169]
   - sout[177] + xorout[128] - xorout[129] >= -4
 R5108: sout[137] - sout[145] - sout[153] + sout[161] - sout[169]
   - sout[177] + xorout[128] - xorout[129] >= -4
 R5109: - sout[137] + sout[145] - sout[153] + sout[161] - sout[169]
   - sout[177] + xorout[128] - xorout[129] >= -4
 R5110: - sout[137] - sout[145] + sout[153] + sout[161] - sout[169]
   - sout[177] + xorout[128] - xorout[129] >= -4
 R5111: sout[137] - sout[145] - sout[153] - sout[161] + sout[169]
   - sout[177] + xorout[128] - xorout[129] >= -4
 R5112: - sout[137] + sout[145] - sout[153] - sout[161] + sout[169]
   - sout[177] + xorout[128] - xorout[129] >= -4
 R5113: - sout[137] - sout[145] + sout[153] - sout[161] + sout[169]
   - sout[177] + xorout[128] - xorout[129] >= -4
 R5114: - sout[137] - sout[145] - sout[153] + sout[161] + sout[169]
   - sout[177] + xorout[128] - xorout[129] >= -4
 R5115: sout[137] - sout[145] - sout[153] - sout[161] - sout[169]
   + sout[177] + xorout[128] - xorout[129] >= -4
 R5116: - sout[137] + sout[145] - sout[153] - sout[161] - sout[169]
   + sout[177] + xorout[128] - xorout[129] >= -4
 R5117: - sout[137] - sout[145] + sout[153] - sout[161] - sout[169]
   + sout[177] + xorout[128] - xorout[129] >= -4
 R5118: - sout[137] - sout[145] - sout[153] + sout[161] - sout[169]
   + sout[177] + xorout[128] - xorout[129] >= -4
 R5119: - sout[137] - sout[145] - sout[153] - sout[161] + sout[169]
   + sout[177] + xorout[128] - xorout[129] >= -4
 R5120: sout[137] + sout[145] + sout[153] + sout[161] - sout[169]
   - sout[177] + xorout[128] - xorout[129] >= -2
 R5121: sout[137] + sout[145] + sout[153] - sout[161] + sout[169]
   - sout[177] + xorout[128] - xorout[129] >= -2
 R5122: sout[137] + sout[145] - sout[153] + sout[161] + sout[169]
   - sout[177] + xorout[128] - xorout[129] >= -2
 R5123: sout[137] - sout[145] + sout[153] + sout[161] + sout[169]
   - sout[177] + xorout[128] - xorout[129] >= -2
 R5124: - sout[137] + sout[145] + sout[153] + sout[161] + sout[169]
   - sout[177] + xorout[128] - xorout[129] >= -2
 R5125: sout[137] + sout[145] + sout[153] - sout[161] - sout[169]
   + sout[177] + xorout[128] - xorout[129] >= -2
 R5126: sout[137] + sout[145] - sout[153] + sout[161] - sout[169]
   + sout[177] + xorout[128] - xorout[129] >= -2
 R5127: sout[137] - sout[145] + sout[153] + sout[161] - sout[169]
   + sout[177] + xorout[128] - xorout[129] >= -2
 R5128: - sout[137] + sout[145] + sout[153] + sout[161] - sout[169]
   + sout[177] + xorout[128] - xorout[129] >= -2
 R5129: sout[137] + sout[145] - sout[153] - sout[161] + sout[169]
   + sout[177] + xorout[128] - xorout[129] >= -2
 R5130: sout[137] - sout[145] + sout[153] - sout[161] + sout[169]
   + sout[177] + xorout[128] - xorout[129] >= -2
 R5131: - sout[137] + sout[145] + sout[153] - sout[161] + sout[169]
   + sout[177] + xorout[128] - xorout[129] >= -2
 R5132: sout[137] - sout[145] - sout[153] + sout[161] + sout[169]
   + sout[177] + xorout[128] - xorout[129] >= -2
 R5133: - sout[137] + sout[145] - sout[153] + sout[161] + sout[169]
   + sout[177] + xorout[128] - xorout[129] >= -2
 R5134: - sout[137] - sout[145] + sout[153] + sout[161] + sout[169]
   + sout[177] + xorout[128] - xorout[129] >= -2
 R5135: sout[139] + sout[147] + sout[155] + sout[163] + sout[171]
   - sout[179] + xorout[131] >= 0
 R5136: sout[139] + sout[147] + sout[155] + sout[163] - sout[171]
   + sout[179] + xorout[131] >= 0
 R5137: sout[139] + sout[147] + sout[155] - sout[163] + sout[171]
   + sout[179] + xorout[131] >= 0
 R5138: sout[139] + sout[147] - sout[155] + sout[163] + sout[171]
   + sout[179] + xorout[131] >= 0
 R5139: sout[139] - sout[147] + sout[155] + sout[163] + sout[171]
   + sout[179] + xorout[131] >= 0
 R5140: - sout[139] + sout[147] + sout[155] + sout[163] + sout[171]
   + sout[179] + xorout[131] >= 0
 R5141: sout[139] + sout[147] + sout[155] + sout[163] + sout[171]
   + sout[179] - xorout[131] >= 0
 R5142: sout[139] + sout[147] + sout[155] - sout[163] - sout[171]
   - sout[179] + xorout[131] >= -2
 R5143: sout[139] + sout[147] - sout[155] + sout[163] - sout[171]
   - sout[179] + xorout[131] >= -2
 R5144: sout[139] - sout[147] + sout[155] + sout[163] - sout[171]
   - sout[179] + xorout[131] >= -2
 R5145: - sout[139] + sout[147] + sout[155] + sout[163] - sout[171]
   - sout[179] + xorout[131] >= -2
 R5146: sout[139] + sout[147] - sout[155] - sout[163] + sout[171]
   - sout[179] + xorout[131] >= -2
 R5147: sout[139] - sout[147] + sout[155] - sout[163] + sout[171]
   - sout[179] + xorout[131] >= -2
 R5148: - sout[139] + sout[147] + sout[155] - sout[163] + sout[171]
   - sout[179] + xorout[131] >= -2
 R5149: sout[139] - sout[147] - sout[155] + sout[163] + sout[171]
   - sout[179] + xorout[131] >= -2
 R5150: - sout[139] + sout[147] - sout[155] + sout[163] + sout[171]
   - sout[179] + xorout[131] >= -2
 R5151: - sout[139] - sout[147] + sout[155] + sout[163] + sout[171]
   - sout[179] + xorout[131] >= -2
 R5152: sout[139] + sout[147] - sout[155] - sout[163] - sout[171]
   + sout[179] + xorout[131] >= -2
 R5153: sout[139] - sout[147] + sout[155] - sout[163] - sout[171]
   + sout[179] + xorout[131] >= -2
 R5154: - sout[139] + sout[147] + sout[155] - sout[163] - sout[171]
   + sout[179] + xorout[131] >= -2
 R5155: sout[139] - sout[147] - sout[155] + sout[163] - sout[171]
   + sout[179] + xorout[131] >= -2
 R5156: - sout[139] + sout[147] - sout[155] + sout[163] - sout[171]
   + sout[179] + xorout[131] >= -2
 R5157: - sout[139] - sout[147] + sout[155] + sout[163] - sout[171]
   + sout[179] + xorout[131] >= -2
 R5158: sout[139] - sout[147] - sout[155] - sout[163] + sout[171]
   + sout[179] + xorout[131] >= -2
 R5159: - sout[139] + sout[147] - sout[155] - sout[163] + sout[171]
   + sout[179] + xorout[131] >= -2
 R5160: - sout[139] - sout[147] + sout[155] - sout[163] + sout[171]
   + sout[179] + xorout[131] >= -2
 R5161: - sout[139] - sout[147] - sout[155] + sout[163] + sout[171]
   + sout[179] + xorout[131] >= -2
 R5162: sout[139] - sout[147] - sout[155] - sout[163] - sout[171]
   - sout[179] + xorout[131] >= -4
 R5163: - sout[139] + sout[147] - sout[155] - sout[163] - sout[171]
   - sout[179] + xorout[131] >= -4
 R5164: - sout[139] - sout[147] + sout[155] - sout[163] - sout[171]
   - sout[179] + xorout[131] >= -4
 R5165: - sout[139] - sout[147] - sout[155] + sout[163] - sout[171]
   - sout[179] + xorout[131] >= -4
 R5166: - sout[139] - sout[147] - sout[155] - sout[163] + sout[171]
   - sout[179] + xorout[131] >= -4
 R5167: - sout[139] - sout[147] - sout[155] - sout[163] - sout[171]
   + sout[179] + xorout[131] >= -4
 R5168: - sout[178] + xorout[130] >= 0
 R5169: - sout[170] + xorout[130] >= 0
 R5170: - sout[162] + xorout[130] >= 0
 R5171: - sout[154] + xorout[130] >= 0
 R5172: - sout[146] + xorout[130] >= 0
 R5173: - sout[138] + xorout[130] >= 0
 R5174: - sout[138] + sout[139] >= 0
 R5175: - sout[146] + sout[147] >= 0
 R5176: - sout[154] + sout[155] >= 0
 R5177: - sout[162] + sout[163] >= 0
 R5178: - sout[170] + sout[171] >= 0
 R5179: - sout[178] + sout[179] >= 0
 R5180: sout[138] + sout[146] + sout[154] + sout[162] + sout[170]
   + sout[178] - xorout[130] >= 0
 R5181: - xorout[130] + xorout[131] >= 0
 R5182: - sout[139] - sout[147] - sout[155] - sout[163] - sout[171]
   - sout[179] + xorout[130] - xorout[131] >= -6
 R5183: sout[139] + sout[147] - sout[155] - sout[163] - sout[171]
   - sout[179] + xorout[130] - xorout[131] >= -4
 R5184: sout[139] - sout[147] + sout[155] - sout[163] - sout[171]
   - sout[179] + xorout[130] - xorout[131] >= -4
 R5185: - sout[139] + sout[147] + sout[155] - sout[163] - sout[171]
   - sout[179] + xorout[130] - xorout[131] >= -4
 R5186: sout[139] - sout[147] - sout[155] + sout[163] - sout[171]
   - sout[179] + xorout[130] - xorout[131] >= -4
 R5187: - sout[139] + sout[147] - sout[155] + sout[163] - sout[171]
   - sout[179] + xorout[130] - xorout[131] >= -4
 R5188: - sout[139] - sout[147] + sout[155] + sout[163] - sout[171]
   - sout[179] + xorout[130] - xorout[131] >= -4
 R5189: sout[139] - sout[147] - sout[155] - sout[163] + sout[171]
   - sout[179] + xorout[130] - xorout[131] >= -4
 R5190: - sout[139] + sout[147] - sout[155] - sout[163] + sout[171]
   - sout[179] + xorout[130] - xorout[131] >= -4
 R5191: - sout[139] - sout[147] + sout[155] - sout[163] + sout[171]
   - sout[179] + xorout[130] - xorout[131] >= -4
 R5192: - sout[139] - sout[147] - sout[155] + sout[163] + sout[171]
   - sout[179] + xorout[130] - xorout[131] >= -4
 R5193: sout[139] - sout[147] - sout[155] - sout[163] - sout[171]
   + sout[179] + xorout[130] - xorout[131] >= -4
 R5194: - sout[139] + sout[147] - sout[155] - sout[163] - sout[171]
   + sout[179] + xorout[130] - xorout[131] >= -4
 R5195: - sout[139] - sout[147] + sout[155] - sout[163] - sout[171]
   + sout[179] + xorout[130] - xorout[131] >= -4
 R5196: - sout[139] - sout[147] - sout[155] + sout[163] - sout[171]
   + sout[179] + xorout[130] - xorout[131] >= -4
 R5197: - sout[139] - sout[147] - sout[155] - sout[163] + sout[171]
   + sout[179] + xorout[130] - xorout[131] >= -4
 R5198: sout[139] + sout[147] + sout[155] + sout[163] - sout[171]
   - sout[179] + xorout[130] - xorout[131] >= -2
 R5199: sout[139] + sout[147] + sout[155] - sout[163] + sout[171]
   - sout[179] + xorout[130] - xorout[131] >= -2
 R5200: sout[139] + sout[147] - sout[155] + sout[163] + sout[171]
   - sout[179] + xorout[130] - xorout[131] >= -2
 R5201: sout[139] - sout[147] + sout[155] + sout[163] + sout[171]
   - sout[179] + xorout[130] - xorout[131] >= -2
 R5202: - sout[139] + sout[147] + sout[155] + sout[163] + sout[171]
   - sout[179] + xorout[130] - xorout[131] >= -2
 R5203: sout[139] + sout[147] + sout[155] - sout[163] - sout[171]
   + sout[179] + xorout[130] - xorout[131] >= -2
 R5204: sout[139] + sout[147] - sout[155] + sout[163] - sout[171]
   + sout[179] + xorout[130] - xorout[131] >= -2
 R5205: sout[139] - sout[147] + sout[155] + sout[163] - sout[171]
   + sout[179] + xorout[130] - xorout[131] >= -2
 R5206: - sout[139] + sout[147] + sout[155] + sout[163] - sout[171]
   + sout[179] + xorout[130] - xorout[131] >= -2
 R5207: sout[139] + sout[147] - sout[155] - sout[163] + sout[171]
   + sout[179] + xorout[130] - xorout[131] >= -2
 R5208: sout[139] - sout[147] + sout[155] - sout[163] + sout[171]
   + sout[179] + xorout[130] - xorout[131] >= -2
 R5209: - sout[139] + sout[147] + sout[155] - sout[163] + sout[171]
   + sout[179] + xorout[130] - xorout[131] >= -2
 R5210: sout[139] - sout[147] - sout[155] + sout[163] + sout[171]
   + sout[179] + xorout[130] - xorout[131] >= -2
 R5211: - sout[139] + sout[147] - sout[155] + sout[163] + sout[171]
   + sout[179] + xorout[130] - xorout[131] >= -2
 R5212: - sout[139] - sout[147] + sout[155] + sout[163] + sout[171]
   + sout[179] + xorout[130] - xorout[131] >= -2
 R5213: sout[141] + sout[149] + sout[157] + sout[165] + sout[173]
   - sout[181] + xorout[133] >= 0
 R5214: sout[141] + sout[149] + sout[157] + sout[165] - sout[173]
   + sout[181] + xorout[133] >= 0
 R5215: sout[141] + sout[149] + sout[157] - sout[165] + sout[173]
   + sout[181] + xorout[133] >= 0
 R5216: sout[141] + sout[149] - sout[157] + sout[165] + sout[173]
   + sout[181] + xorout[133] >= 0
 R5217: sout[141] - sout[149] + sout[157] + sout[165] + sout[173]
   + sout[181] + xorout[133] >= 0
 R5218: - sout[141] + sout[149] + sout[157] + sout[165] + sout[173]
   + sout[181] + xorout[133] >= 0
 R5219: sout[141] + sout[149] + sout[157] + sout[165] + sout[173]
   + sout[181] - xorout[133] >= 0
 R5220: sout[141] + sout[149] + sout[157] - sout[165] - sout[173]
   - sout[181] + xorout[133] >= -2
 R5221: sout[141] + sout[149] - sout[157] + sout[165] - sout[173]
   - sout[181] + xorout[133] >= -2
 R5222: sout[141] - sout[149] + sout[157] + sout[165] - sout[173]
   - sout[181] + xorout[133] >= -2
 R5223: - sout[141] + sout[149] + sout[157] + sout[165] - sout[173]
   - sout[181] + xorout[133] >= -2
 R5224: sout[141] + sout[149] - sout[157] - sout[165] + sout[173]
   - sout[181] + xorout[133] >= -2
 R5225: sout[141] - sout[149] + sout[157] - sout[165] + sout[173]
   - sout[181] + xorout[133] >= -2
 R5226: - sout[141] + sout[149] + sout[157] - sout[165] + sout[173]
   - sout[181] + xorout[133] >= -2
 R5227: sout[141] - sout[149] - sout[157] + sout[165] + sout[173]
   - sout[181] + xorout[133] >= -2
 R5228: - sout[141] + sout[149] - sout[157] + sout[165] + sout[173]
   - sout[181] + xorout[133] >= -2
 R5229: - sout[141] - sout[149] + sout[157] + sout[165] + sout[173]
   - sout[181] + xorout[133] >= -2
 R5230: sout[141] + sout[149] - sout[157] - sout[165] - sout[173]
   + sout[181] + xorout[133] >= -2
 R5231: sout[141] - sout[149] + sout[157] - sout[165] - sout[173]
   + sout[181] + xorout[133] >= -2
 R5232: - sout[141] + sout[149] + sout[157] - sout[165] - sout[173]
   + sout[181] + xorout[133] >= -2
 R5233: sout[141] - sout[149] - sout[157] + sout[165] - sout[173]
   + sout[181] + xorout[133] >= -2
 R5234: - sout[141] + sout[149] - sout[157] + sout[165] - sout[173]
   + sout[181] + xorout[133] >= -2
 R5235: - sout[141] - sout[149] + sout[157] + sout[165] - sout[173]
   + sout[181] + xorout[133] >= -2
 R5236: sout[141] - sout[149] - sout[157] - sout[165] + sout[173]
   + sout[181] + xorout[133] >= -2
 R5237: - sout[141] + sout[149] - sout[157] - sout[165] + sout[173]
   + sout[181] + xorout[133] >= -2
 R5238: - sout[141] - sout[149] + sout[157] - sout[165] + sout[173]
   + sout[181] + xorout[133] >= -2
 R5239: - sout[141] - sout[149] - sout[157] + sout[165] + sout[173]
   + sout[181] + xorout[133] >= -2
 R5240: sout[141] - sout[149] - sout[157] - sout[165] - sout[173]
   - sout[181] + xorout[133] >= -4
 R5241: - sout[141] + sout[149] - sout[157] - sout[165] - sout[173]
   - sout[181] + xorout[133] >= -4
 R5242: - sout[141] - sout[149] + sout[157] - sout[165] - sout[173]
   - sout[181] + xorout[133] >= -4
 R5243: - sout[141] - sout[149] - sout[157] + sout[165] - sout[173]
   - sout[181] + xorout[133] >= -4
 R5244: - sout[141] - sout[149] - sout[157] - sout[165] + sout[173]
   - sout[181] + xorout[133] >= -4
 R5245: - sout[141] - sout[149] - sout[157] - sout[165] - sout[173]
   + sout[181] + xorout[133] >= -4
 R5246: - sout[180] + xorout[132] >= 0
 R5247: - sout[172] + xorout[132] >= 0
 R5248: - sout[164] + xorout[132] >= 0
 R5249: - sout[156] + xorout[132] >= 0
 R5250: - sout[148] + xorout[132] >= 0
 R5251: - sout[140] + xorout[132] >= 0
 R5252: - sout[140] + sout[141] >= 0
 R5253: - sout[148] + sout[149] >= 0
 R5254: - sout[156] + sout[157] >= 0
 R5255: - sout[164] + sout[165] >= 0
 R5256: - sout[172] + sout[173] >= 0
 R5257: - sout[180] + sout[181] >= 0
 R5258: sout[140] + sout[148] + sout[156] + sout[164] + sout[172]
   + sout[180] - xorout[132] >= 0
 R5259: - xorout[132] + xorout[133] >= 0
 R5260: - sout[141] - sout[149] - sout[157] - sout[165] - sout[173]
   - sout[181] + xorout[132] - xorout[133] >= -6
 R5261: sout[141] + sout[149] - sout[157] - sout[165] - sout[173]
   - sout[181] + xorout[132] - xorout[133] >= -4
 R5262: sout[141] - sout[149] + sout[157] - sout[165] - sout[173]
   - sout[181] + xorout[132] - xorout[133] >= -4
 R5263: - sout[141] + sout[149] + sout[157] - sout[165] - sout[173]
   - sout[181] + xorout[132] - xorout[133] >= -4
 R5264: sout[141] - sout[149] - sout[157] + sout[165] - sout[173]
   - sout[181] + xorout[132] - xorout[133] >= -4
 R5265: - sout[141] + sout[149] - sout[157] + sout[165] - sout[173]
   - sout[181] + xorout[132] - xorout[133] >= -4
 R5266: - sout[141] - sout[149] + sout[157] + sout[165] - sout[173]
   - sout[181] + xorout[132] - xorout[133] >= -4
 R5267: sout[141] - sout[149] - sout[157] - sout[165] + sout[173]
   - sout[181] + xorout[132] - xorout[133] >= -4
 R5268: - sout[141] + sout[149] - sout[157] - sout[165] + sout[173]
   - sout[181] + xorout[132] - xorout[133] >= -4
 R5269: - sout[141] - sout[149] + sout[157] - sout[165] + sout[173]
   - sout[181] + xorout[132] - xorout[133] >= -4
 R5270: - sout[141] - sout[149] - sout[157] + sout[165] + sout[173]
   - sout[181] + xorout[132] - xorout[133] >= -4
 R5271: sout[141] - sout[149] - sout[157] - sout[165] - sout[173]
   + sout[181] + xorout[132] - xorout[133] >= -4
 R5272: - sout[141] + sout[149] - sout[157] - sout[165] - sout[173]
   + sout[181] + xorout[132] - xorout[133] >= -4
 R5273: - sout[141] - sout[149] + sout[157] - sout[165] - sout[173]
   + sout[181] + xorout[132] - xorout[133] >= -4
 R5274: - sout[141] - sout[149] - sout[157] + sout[165] - sout[173]
   + sout[181] + xorout[132] - xorout[133] >= -4
 R5275: - sout[141] - sout[149] - sout[157] - sout[165] + sout[173]
   + sout[181] + xorout[132] - xorout[133] >= -4
 R5276: sout[141] + sout[149] + sout[157] + sout[165] - sout[173]
   - sout[181] + xorout[132] - xorout[133] >= -2
 R5277: sout[141] + sout[149] + sout[157] - sout[165] + sout[173]
   - sout[181] + xorout[132] - xorout[133] >= -2
 R5278: sout[141] + sout[149] - sout[157] + sout[165] + sout[173]
   - sout[181] + xorout[132] - xorout[133] >= -2
 R5279: sout[141] - sout[149] + sout[157] + sout[165] + sout[173]
   - sout[181] + xorout[132] - xorout[133] >= -2
 R5280: - sout[141] + sout[149] + sout[157] + sout[165] + sout[173]
   - sout[181] + xorout[132] - xorout[133] >= -2
 R5281: sout[141] + sout[149] + sout[157] - sout[165] - sout[173]
   + sout[181] + xorout[132] - xorout[133] >= -2
 R5282: sout[141] + sout[149] - sout[157] + sout[165] - sout[173]
   + sout[181] + xorout[132] - xorout[133] >= -2
 R5283: sout[141] - sout[149] + sout[157] + sout[165] - sout[173]
   + sout[181] + xorout[132] - xorout[133] >= -2
 R5284: - sout[141] + sout[149] + sout[157] + sout[165] - sout[173]
   + sout[181] + xorout[132] - xorout[133] >= -2
 R5285: sout[141] + sout[149] - sout[157] - sout[165] + sout[173]
   + sout[181] + xorout[132] - xorout[133] >= -2
 R5286: sout[141] - sout[149] + sout[157] - sout[165] + sout[173]
   + sout[181] + xorout[132] - xorout[133] >= -2
 R5287: - sout[141] + sout[149] + sout[157] - sout[165] + sout[173]
   + sout[181] + xorout[132] - xorout[133] >= -2
 R5288: sout[141] - sout[149] - sout[157] + sout[165] + sout[173]
   + sout[181] + xorout[132] - xorout[133] >= -2
 R5289: - sout[141] + sout[149] - sout[157] + sout[165] + sout[173]
   + sout[181] + xorout[132] - xorout[133] >= -2
 R5290: - sout[141] - sout[149] + sout[157] + sout[165] + sout[173]
   + sout[181] + xorout[132] - xorout[133] >= -2
 R5291: sout[143] + sout[151] + sout[159] + sout[167] + sout[175]
   - sout[183] + xorout[135] >= 0
 R5292: sout[143] + sout[151] + sout[159] + sout[167] - sout[175]
   + sout[183] + xorout[135] >= 0
 R5293: sout[143] + sout[151] + sout[159] - sout[167] + sout[175]
   + sout[183] + xorout[135] >= 0
 R5294: sout[143] + sout[151] - sout[159] + sout[167] + sout[175]
   + sout[183] + xorout[135] >= 0
 R5295: sout[143] - sout[151] + sout[159] + sout[167] + sout[175]
   + sout[183] + xorout[135] >= 0
 R5296: - sout[143] + sout[151] + sout[159] + sout[167] + sout[175]
   + sout[183] + xorout[135] >= 0
 R5297: sout[143] + sout[151] + sout[159] + sout[167] + sout[175]
   + sout[183] - xorout[135] >= 0
 R5298: sout[143] + sout[151] + sout[159] - sout[167] - sout[175]
   - sout[183] + xorout[135] >= -2
 R5299: sout[143] + sout[151] - sout[159] + sout[167] - sout[175]
   - sout[183] + xorout[135] >= -2
 R5300: sout[143] - sout[151] + sout[159] + sout[167] - sout[175]
   - sout[183] + xorout[135] >= -2
 R5301: - sout[143] + sout[151] + sout[159] + sout[167] - sout[175]
   - sout[183] + xorout[135] >= -2
 R5302: sout[143] + sout[151] - sout[159] - sout[167] + sout[175]
   - sout[183] + xorout[135] >= -2
 R5303: sout[143] - sout[151] + sout[159] - sout[167] + sout[175]
   - sout[183] + xorout[135] >= -2
 R5304: - sout[143] + sout[151] + sout[159] - sout[167] + sout[175]
   - sout[183] + xorout[135] >= -2
 R5305: sout[143] - sout[151] - sout[159] + sout[167] + sout[175]
   - sout[183] + xorout[135] >= -2
 R5306: - sout[143] + sout[151] - sout[159] + sout[167] + sout[175]
   - sout[183] + xorout[135] >= -2
 R5307: - sout[143] - sout[151] + sout[159] + sout[167] + sout[175]
   - sout[183] + xorout[135] >= -2
 R5308: sout[143] + sout[151] - sout[159] - sout[167] - sout[175]
   + sout[183] + xorout[135] >= -2
 R5309: sout[143] - sout[151] + sout[159] - sout[167] - sout[175]
   + sout[183] + xorout[135] >= -2
 R5310: - sout[143] + sout[151] + sout[159] - sout[167] - sout[175]
   + sout[183] + xorout[135] >= -2
 R5311: sout[143] - sout[151] - sout[159] + sout[167] - sout[175]
   + sout[183] + xorout[135] >= -2
 R5312: - sout[143] + sout[151] - sout[159] + sout[167] - sout[175]
   + sout[183] + xorout[135] >= -2
 R5313: - sout[143] - sout[151] + sout[159] + sout[167] - sout[175]
   + sout[183] + xorout[135] >= -2
 R5314: sout[143] - sout[151] - sout[159] - sout[167] + sout[175]
   + sout[183] + xorout[135] >= -2
 R5315: - sout[143] + sout[151] - sout[159] - sout[167] + sout[175]
   + sout[183] + xorout[135] >= -2
 R5316: - sout[143] - sout[151] + sout[159] - sout[167] + sout[175]
   + sout[183] + xorout[135] >= -2
 R5317: - sout[143] - sout[151] - sout[159] + sout[167] + sout[175]
   + sout[183] + xorout[135] >= -2
 R5318: sout[143] - sout[151] - sout[159] - sout[167] - sout[175]
   - sout[183] + xorout[135] >= -4
 R5319: - sout[143] + sout[151] - sout[159] - sout[167] - sout[175]
   - sout[183] + xorout[135] >= -4
 R5320: - sout[143] - sout[151] + sout[159] - sout[167] - sout[175]
   - sout[183] + xorout[135] >= -4
 R5321: - sout[143] - sout[151] - sout[159] + sout[167] - sout[175]
   - sout[183] + xorout[135] >= -4
 R5322: - sout[143] - sout[151] - sout[159] - sout[167] + sout[175]
   - sout[183] + xorout[135] >= -4
 R5323: - sout[143] - sout[151] - sout[159] - sout[167] - sout[175]
   + sout[183] + xorout[135] >= -4
 R5324: - sout[182] + xorout[134] >= 0
 R5325: - sout[174] + xorout[134] >= 0
 R5326: - sout[166] + xorout[134] >= 0
 R5327: - sout[158] + xorout[134] >= 0
 R5328: - sout[150] + xorout[134] >= 0
 R5329: - sout[142] + xorout[134] >= 0
 R5330: - sout[142] + sout[143] >= 0
 R5331: - sout[150] + sout[151] >= 0
 R5332: - sout[158] + sout[159] >= 0
 R5333: - sout[166] + sout[167] >= 0
 R5334: - sout[174] + sout[175] >= 0
 R5335: - sout[182] + sout[183] >= 0
 R5336: sout[142] + sout[150] + sout[158] + sout[166] + sout[174]
   + sout[182] - xorout[134] >= 0
 R5337: - xorout[134] + xorout[135] >= 0
 R5338: - sout[143] - sout[151] - sout[159] - sout[167] - sout[175]
   - sout[183] + xorout[134] - xorout[135] >= -6
 R5339: sout[143] + sout[151] - sout[159] - sout[167] - sout[175]
   - sout[183] + xorout[134] - xorout[135] >= -4
 R5340: sout[143] - sout[151] + sout[159] - sout[167] - sout[175]
   - sout[183] + xorout[134] - xorout[135] >= -4
 R5341: - sout[143] + sout[151] + sout[159] - sout[167] - sout[175]
   - sout[183] + xorout[134] - xorout[135] >= -4
 R5342: sout[143] - sout[151] - sout[159] + sout[167] - sout[175]
   - sout[183] + xorout[134] - xorout[135] >= -4
 R5343: - sout[143] + sout[151] - sout[159] + sout[167] - sout[175]
   - sout[183] + xorout[134] - xorout[135] >= -4
 R5344: - sout[143] - sout[151] + sout[159] + sout[167] - sout[175]
   - sout[183] + xorout[134] - xorout[135] >= -4
 R5345: sout[143] - sout[151] - sout[159] - sout[167] + sout[175]
   - sout[183] + xorout[134] - xorout[135] >= -4
 R5346: - sout[143] + sout[151] - sout[159] - sout[167] + sout[175]
   - sout[183] + xorout[134] - xorout[135] >= -4
 R5347: - sout[143] - sout[151] + sout[159] - sout[167] + sout[175]
   - sout[183] + xorout[134] - xorout[135] >= -4
 R5348: - sout[143] - sout[151] - sout[159] + sout[167] + sout[175]
   - sout[183] + xorout[134] - xorout[135] >= -4
 R5349: sout[143] - sout[151] - sout[159] - sout[167] - sout[175]
   + sout[183] + xorout[134] - xorout[135] >= -4
 R5350: - sout[143] + sout[151] - sout[159] - sout[167] - sout[175]
   + sout[183] + xorout[134] - xorout[135] >= -4
 R5351: - sout[143] - sout[151] + sout[159] - sout[167] - sout[175]
   + sout[183] + xorout[134] - xorout[135] >= -4
 R5352: - sout[143] - sout[151] - sout[159] + sout[167] - sout[175]
   + sout[183] + xorout[134] - xorout[135] >= -4
 R5353: - sout[143] - sout[151] - sout[159] - sout[167] + sout[175]
   + sout[183] + xorout[134] - xorout[135] >= -4
 R5354: sout[143] + sout[151] + sout[159] + sout[167] - sout[175]
   - sout[183] + xorout[134] - xorout[135] >= -2
 R5355: sout[143] + sout[151] + sout[159] - sout[167] + sout[175]
   - sout[183] + xorout[134] - xorout[135] >= -2
 R5356: sout[143] + sout[151] - sout[159] + sout[167] + sout[175]
   - sout[183] + xorout[134] - xorout[135] >= -2
 R5357: sout[143] - sout[151] + sout[159] + sout[167] + sout[175]
   - sout[183] + xorout[134] - xorout[135] >= -2
 R5358: - sout[143] + sout[151] + sout[159] + sout[167] + sout[175]
   - sout[183] + xorout[134] - xorout[135] >= -2
 R5359: sout[143] + sout[151] + sout[159] - sout[167] - sout[175]
   + sout[183] + xorout[134] - xorout[135] >= -2
 R5360: sout[143] + sout[151] - sout[159] + sout[167] - sout[175]
   + sout[183] + xorout[134] - xorout[135] >= -2
 R5361: sout[143] - sout[151] + sout[159] + sout[167] - sout[175]
   + sout[183] + xorout[134] - xorout[135] >= -2
 R5362: - sout[143] + sout[151] + sout[159] + sout[167] - sout[175]
   + sout[183] + xorout[134] - xorout[135] >= -2
 R5363: sout[143] + sout[151] - sout[159] - sout[167] + sout[175]
   + sout[183] + xorout[134] - xorout[135] >= -2
 R5364: sout[143] - sout[151] + sout[159] - sout[167] + sout[175]
   + sout[183] + xorout[134] - xorout[135] >= -2
 R5365: - sout[143] + sout[151] + sout[159] - sout[167] + sout[175]
   + sout[183] + xorout[134] - xorout[135] >= -2
 R5366: sout[143] - sout[151] - sout[159] + sout[167] + sout[175]
   + sout[183] + xorout[134] - xorout[135] >= -2
 R5367: - sout[143] + sout[151] - sout[159] + sout[167] + sout[175]
   + sout[183] + xorout[134] - xorout[135] >= -2
 R5368: - sout[143] - sout[151] + sout[159] + sout[167] + sout[175]
   + sout[183] + xorout[134] - xorout[135] >= -2
 R5369: sout[129] + sout[145] + sout[153] + sout[169] + sout[177]
   - sout[185] + xorout[137] >= 0
 R5370: sout[129] + sout[145] + sout[153] + sout[169] - sout[177]
   + sout[185] + xorout[137] >= 0
 R5371: sout[129] + sout[145] + sout[153] - sout[169] + sout[177]
   + sout[185] + xorout[137] >= 0
 R5372: sout[129] + sout[145] - sout[153] + sout[169] + sout[177]
   + sout[185] + xorout[137] >= 0
 R5373: sout[129] - sout[145] + sout[153] + sout[169] + sout[177]
   + sout[185] + xorout[137] >= 0
 R5374: - sout[129] + sout[145] + sout[153] + sout[169] + sout[177]
   + sout[185] + xorout[137] >= 0
 R5375: sout[129] + sout[145] + sout[153] + sout[169] + sout[177]
   + sout[185] - xorout[137] >= 0
 R5376: sout[129] + sout[145] + sout[153] - sout[169] - sout[177]
   - sout[185] + xorout[137] >= -2
 R5377: sout[129] + sout[145] - sout[153] + sout[169] - sout[177]
   - sout[185] + xorout[137] >= -2
 R5378: sout[129] - sout[145] + sout[153] + sout[169] - sout[177]
   - sout[185] + xorout[137] >= -2
 R5379: - sout[129] + sout[145] + sout[153] + sout[169] - sout[177]
   - sout[185] + xorout[137] >= -2
 R5380: sout[129] + sout[145] - sout[153] - sout[169] + sout[177]
   - sout[185] + xorout[137] >= -2
 R5381: sout[129] - sout[145] + sout[153] - sout[169] + sout[177]
   - sout[185] + xorout[137] >= -2
 R5382: - sout[129] + sout[145] + sout[153] - sout[169] + sout[177]
   - sout[185] + xorout[137] >= -2
 R5383: sout[129] - sout[145] - sout[153] + sout[169] + sout[177]
   - sout[185] + xorout[137] >= -2
 R5384: - sout[129] + sout[145] - sout[153] + sout[169] + sout[177]
   - sout[185] + xorout[137] >= -2
 R5385: - sout[129] - sout[145] + sout[153] + sout[169] + sout[177]
   - sout[185] + xorout[137] >= -2
 R5386: sout[129] + sout[145] - sout[153] - sout[169] - sout[177]
   + sout[185] + xorout[137] >= -2
 R5387: sout[129] - sout[145] + sout[153] - sout[169] - sout[177]
   + sout[185] + xorout[137] >= -2
 R5388: - sout[129] + sout[145] + sout[153] - sout[169] - sout[177]
   + sout[185] + xorout[137] >= -2
 R5389: sout[129] - sout[145] - sout[153] + sout[169] - sout[177]
   + sout[185] + xorout[137] >= -2
 R5390: - sout[129] + sout[145] - sout[153] + sout[169] - sout[177]
   + sout[185] + xorout[137] >= -2
 R5391: - sout[129] - sout[145] + sout[153] + sout[169] - sout[177]
   + sout[185] + xorout[137] >= -2
 R5392: sout[129] - sout[145] - sout[153] - sout[169] + sout[177]
   + sout[185] + xorout[137] >= -2
 R5393: - sout[129] + sout[145] - sout[153] - sout[169] + sout[177]
   + sout[185] + xorout[137] >= -2
 R5394: - sout[129] - sout[145] + sout[153] - sout[169] + sout[177]
   + sout[185] + xorout[137] >= -2
 R5395: - sout[129] - sout[145] - sout[153] + sout[169] + sout[177]
   + sout[185] + xorout[137] >= -2
 R5396: sout[129] - sout[145] - sout[153] - sout[169] - sout[177]
   - sout[185] + xorout[137] >= -4
 R5397: - sout[129] + sout[145] - sout[153] - sout[169] - sout[177]
   - sout[185] + xorout[137] >= -4
 R5398: - sout[129] - sout[145] + sout[153] - sout[169] - sout[177]
   - sout[185] + xorout[137] >= -4
 R5399: - sout[129] - sout[145] - sout[153] + sout[169] - sout[177]
   - sout[185] + xorout[137] >= -4
 R5400: - sout[129] - sout[145] - sout[153] - sout[169] + sout[177]
   - sout[185] + xorout[137] >= -4
 R5401: - sout[129] - sout[145] - sout[153] - sout[169] - sout[177]
   + sout[185] + xorout[137] >= -4
 R5402: - sout[184] + xorout[136] >= 0
 R5403: - sout[176] + xorout[136] >= 0
 R5404: - sout[168] + xorout[136] >= 0
 R5405: - sout[152] + xorout[136] >= 0
 R5406: - sout[144] + xorout[136] >= 0
 R5407: - sout[128] + xorout[136] >= 0
 R5408: - sout[128] + sout[129] >= 0
 R5409: - sout[144] + sout[145] >= 0
 R5410: - sout[152] + sout[153] >= 0
 R5411: - sout[168] + sout[169] >= 0
 R5412: - sout[176] + sout[177] >= 0
 R5413: - sout[184] + sout[185] >= 0
 R5414: sout[128] + sout[144] + sout[152] + sout[168] + sout[176]
   + sout[184] - xorout[136] >= 0
 R5415: - xorout[136] + xorout[137] >= 0
 R5416: - sout[129] - sout[145] - sout[153] - sout[169] - sout[177]
   - sout[185] + xorout[136] - xorout[137] >= -6
 R5417: sout[129] + sout[145] - sout[153] - sout[169] - sout[177]
   - sout[185] + xorout[136] - xorout[137] >= -4
 R5418: sout[129] - sout[145] + sout[153] - sout[169] - sout[177]
   - sout[185] + xorout[136] - xorout[137] >= -4
 R5419: - sout[129] + sout[145] + sout[153] - sout[169] - sout[177]
   - sout[185] + xorout[136] - xorout[137] >= -4
 R5420: sout[129] - sout[145] - sout[153] + sout[169] - sout[177]
   - sout[185] + xorout[136] - xorout[137] >= -4
 R5421: - sout[129] + sout[145] - sout[153] + sout[169] - sout[177]
   - sout[185] + xorout[136] - xorout[137] >= -4
 R5422: - sout[129] - sout[145] + sout[153] + sout[169] - sout[177]
   - sout[185] + xorout[136] - xorout[137] >= -4
 R5423: sout[129] - sout[145] - sout[153] - sout[169] + sout[177]
   - sout[185] + xorout[136] - xorout[137] >= -4
 R5424: - sout[129] + sout[145] - sout[153] - sout[169] + sout[177]
   - sout[185] + xorout[136] - xorout[137] >= -4
 R5425: - sout[129] - sout[145] + sout[153] - sout[169] + sout[177]
   - sout[185] + xorout[136] - xorout[137] >= -4
 R5426: - sout[129] - sout[145] - sout[153] + sout[169] + sout[177]
   - sout[185] + xorout[136] - xorout[137] >= -4
 R5427: sout[129] - sout[145] - sout[153] - sout[169] - sout[177]
   + sout[185] + xorout[136] - xorout[137] >= -4
 R5428: - sout[129] + sout[145] - sout[153] - sout[169] - sout[177]
   + sout[185] + xorout[136] - xorout[137] >= -4
 R5429: - sout[129] - sout[145] + sout[153] - sout[169] - sout[177]
   + sout[185] + xorout[136] - xorout[137] >= -4
 R5430: - sout[129] - sout[145] - sout[153] + sout[169] - sout[177]
   + sout[185] + xorout[136] - xorout[137] >= -4
 R5431: - sout[129] - sout[145] - sout[153] - sout[169] + sout[177]
   + sout[185] + xorout[136] - xorout[137] >= -4
 R5432: sout[129] + sout[145] + sout[153] + sout[169] - sout[177]
   - sout[185] + xorout[136] - xorout[137] >= -2
 R5433: sout[129] + sout[145] + sout[153] - sout[169] + sout[177]
   - sout[185] + xorout[136] - xorout[137] >= -2
 R5434: sout[129] + sout[145] - sout[153] + sout[169] + sout[177]
   - sout[185] + xorout[136] - xorout[137] >= -2
 R5435: sout[129] - sout[145] + sout[153] + sout[169] + sout[177]
   - sout[185] + xorout[136] - xorout[137] >= -2
 R5436: - sout[129] + sout[145] + sout[153] + sout[169] + sout[177]
   - sout[185] + xorout[136] - xorout[137] >= -2
 R5437: sout[129] + sout[145] + sout[153] - sout[169] - sout[177]
   + sout[185] + xorout[136] - xorout[137] >= -2
 R5438: sout[129] + sout[145] - sout[153] + sout[169] - sout[177]
   + sout[185] + xorout[136] - xorout[137] >= -2
 R5439: sout[129] - sout[145] + sout[153] + sout[169] - sout[177]
   + sout[185] + xorout[136] - xorout[137] >= -2
 R5440: - sout[129] + sout[145] + sout[153] + sout[169] - sout[177]
   + sout[185] + xorout[136] - xorout[137] >= -2
 R5441: sout[129] + sout[145] - sout[153] - sout[169] + sout[177]
   + sout[185] + xorout[136] - xorout[137] >= -2
 R5442: sout[129] - sout[145] + sout[153] - sout[169] + sout[177]
   + sout[185] + xorout[136] - xorout[137] >= -2
 R5443: - sout[129] + sout[145] + sout[153] - sout[169] + sout[177]
   + sout[185] + xorout[136] - xorout[137] >= -2
 R5444: sout[129] - sout[145] - sout[153] + sout[169] + sout[177]
   + sout[185] + xorout[136] - xorout[137] >= -2
 R5445: - sout[129] + sout[145] - sout[153] + sout[169] + sout[177]
   + sout[185] + xorout[136] - xorout[137] >= -2
 R5446: - sout[129] - sout[145] + sout[153] + sout[169] + sout[177]
   + sout[185] + xorout[136] - xorout[137] >= -2
 R5447: sout[131] + sout[147] + sout[155] + sout[171] + sout[179]
   - sout[187] + xorout[139] >= 0
 R5448: sout[131] + sout[147] + sout[155] + sout[171] - sout[179]
   + sout[187] + xorout[139] >= 0
 R5449: sout[131] + sout[147] + sout[155] - sout[171] + sout[179]
   + sout[187] + xorout[139] >= 0
 R5450: sout[131] + sout[147] - sout[155] + sout[171] + sout[179]
   + sout[187] + xorout[139] >= 0
 R5451: sout[131] - sout[147] + sout[155] + sout[171] + sout[179]
   + sout[187] + xorout[139] >= 0
 R5452: - sout[131] + sout[147] + sout[155] + sout[171] + sout[179]
   + sout[187] + xorout[139] >= 0
 R5453: sout[131] + sout[147] + sout[155] + sout[171] + sout[179]
   + sout[187] - xorout[139] >= 0
 R5454: sout[131] + sout[147] + sout[155] - sout[171] - sout[179]
   - sout[187] + xorout[139] >= -2
 R5455: sout[131] + sout[147] - sout[155] + sout[171] - sout[179]
   - sout[187] + xorout[139] >= -2
 R5456: sout[131] - sout[147] + sout[155] + sout[171] - sout[179]
   - sout[187] + xorout[139] >= -2
 R5457: - sout[131] + sout[147] + sout[155] + sout[171] - sout[179]
   - sout[187] + xorout[139] >= -2
 R5458: sout[131] + sout[147] - sout[155] - sout[171] + sout[179]
   - sout[187] + xorout[139] >= -2
 R5459: sout[131] - sout[147] + sout[155] - sout[171] + sout[179]
   - sout[187] + xorout[139] >= -2
 R5460: - sout[131] + sout[147] + sout[155] - sout[171] + sout[179]
   - sout[187] + xorout[139] >= -2
 R5461: sout[131] - sout[147] - sout[155] + sout[171] + sout[179]
   - sout[187] + xorout[139] >= -2
 R5462: - sout[131] + sout[147] - sout[155] + sout[171] + sout[179]
   - sout[187] + xorout[139] >= -2
 R5463: - sout[131] - sout[147] + sout[155] + sout[171] + sout[179]
   - sout[187] + xorout[139] >= -2
 R5464: sout[131] + sout[147] - sout[155] - sout[171] - sout[179]
   + sout[187] + xorout[139] >= -2
 R5465: sout[131] - sout[147] + sout[155] - sout[171] - sout[179]
   + sout[187] + xorout[139] >= -2
 R5466: - sout[131] + sout[147] + sout[155] - sout[171] - sout[179]
   + sout[187] + xorout[139] >= -2
 R5467: sout[131] - sout[147] - sout[155] + sout[171] - sout[179]
   + sout[187] + xorout[139] >= -2
 R5468: - sout[131] + sout[147] - sout[155] + sout[171] - sout[179]
   + sout[187] + xorout[139] >= -2
 R5469: - sout[131] - sout[147] + sout[155] + sout[171] - sout[179]
   + sout[187] + xorout[139] >= -2
 R5470: sout[131] - sout[147] - sout[155] - sout[171] + sout[179]
   + sout[187] + xorout[139] >= -2
 R5471: - sout[131] + sout[147] - sout[155] - sout[171] + sout[179]
   + sout[187] + xorout[139] >= -2
 R5472: - sout[131] - sout[147] + sout[155] - sout[171] + sout[179]
   + sout[187] + xorout[139] >= -2
 R5473: - sout[131] - sout[147] - sout[155] + sout[171] + sout[179]
   + sout[187] + xorout[139] >= -2
 R5474: sout[131] - sout[147] - sout[155] - sout[171] - sout[179]
   - sout[187] + xorout[139] >= -4
 R5475: - sout[131] + sout[147] - sout[155] - sout[171] - sout[179]
   - sout[187] + xorout[139] >= -4
 R5476: - sout[131] - sout[147] + sout[155] - sout[171] - sout[179]
   - sout[187] + xorout[139] >= -4
 R5477: - sout[131] - sout[147] - sout[155] + sout[171] - sout[179]
   - sout[187] + xorout[139] >= -4
 R5478: - sout[131] - sout[147] - sout[155] - sout[171] + sout[179]
   - sout[187] + xorout[139] >= -4
 R5479: - sout[131] - sout[147] - sout[155] - sout[171] - sout[179]
   + sout[187] + xorout[139] >= -4
 R5480: - sout[186] + xorout[138] >= 0
 R5481: - sout[178] + xorout[138] >= 0
 R5482: - sout[170] + xorout[138] >= 0
 R5483: - sout[154] + xorout[138] >= 0
 R5484: - sout[146] + xorout[138] >= 0
 R5485: - sout[130] + xorout[138] >= 0
 R5486: - sout[130] + sout[131] >= 0
 R5487: - sout[146] + sout[147] >= 0
 R5488: - sout[154] + sout[155] >= 0
 R5489: - sout[170] + sout[171] >= 0
 R5490: - sout[178] + sout[179] >= 0
 R5491: - sout[186] + sout[187] >= 0
 R5492: sout[130] + sout[146] + sout[154] + sout[170] + sout[178]
   + sout[186] - xorout[138] >= 0
 R5493: - xorout[138] + xorout[139] >= 0
 R5494: - sout[131] - sout[147] - sout[155] - sout[171] - sout[179]
   - sout[187] + xorout[138] - xorout[139] >= -6
 R5495: sout[131] + sout[147] - sout[155] - sout[171] - sout[179]
   - sout[187] + xorout[138] - xorout[139] >= -4
 R5496: sout[131] - sout[147] + sout[155] - sout[171] - sout[179]
   - sout[187] + xorout[138] - xorout[139] >= -4
 R5497: - sout[131] + sout[147] + sout[155] - sout[171] - sout[179]
   - sout[187] + xorout[138] - xorout[139] >= -4
 R5498: sout[131] - sout[147] - sout[155] + sout[171] - sout[179]
   - sout[187] + xorout[138] - xorout[139] >= -4
 R5499: - sout[131] + sout[147] - sout[155] + sout[171] - sout[179]
   - sout[187] + xorout[138] - xorout[139] >= -4
 R5500: - sout[131] - sout[147] + sout[155] + sout[171] - sout[179]
   - sout[187] + xorout[138] - xorout[139] >= -4
 R5501: sout[131] - sout[147] - sout[155] - sout[171] + sout[179]
   - sout[187] + xorout[138] - xorout[139] >= -4
 R5502: - sout[131] + sout[147] - sout[155] - sout[171] + sout[179]
   - sout[187] + xorout[138] - xorout[139] >= -4
 R5503: - sout[131] - sout[147] + sout[155] - sout[171] + sout[179]
   - sout[187] + xorout[138] - xorout[139] >= -4
 R5504: - sout[131] - sout[147] - sout[155] + sout[171] + sout[179]
   - sout[187] + xorout[138] - xorout[139] >= -4
 R5505: sout[131] - sout[147] - sout[155] - sout[171] - sout[179]
   + sout[187] + xorout[138] - xorout[139] >= -4
 R5506: - sout[131] + sout[147] - sout[155] - sout[171] - sout[179]
   + sout[187] + xorout[138] - xorout[139] >= -4
 R5507: - sout[131] - sout[147] + sout[155] - sout[171] - sout[179]
   + sout[187] + xorout[138] - xorout[139] >= -4
 R5508: - sout[131] - sout[147] - sout[155] + sout[171] - sout[179]
   + sout[187] + xorout[138] - xorout[139] >= -4
 R5509: - sout[131] - sout[147] - sout[155] - sout[171] + sout[179]
   + sout[187] + xorout[138] - xorout[139] >= -4
 R5510: sout[131] + sout[147] + sout[155] + sout[171] - sout[179]
   - sout[187] + xorout[138] - xorout[139] >= -2
 R5511: sout[131] + sout[147] + sout[155] - sout[171] + sout[179]
   - sout[187] + xorout[138] - xorout[139] >= -2
 R5512: sout[131] + sout[147] - sout[155] + sout[171] + sout[179]
   - sout[187] + xorout[138] - xorout[139] >= -2
 R5513: sout[131] - sout[147] + sout[155] + sout[171] + sout[179]
   - sout[187] + xorout[138] - xorout[139] >= -2
 R5514: - sout[131] + sout[147] + sout[155] + sout[171] + sout[179]
   - sout[187] + xorout[138] - xorout[139] >= -2
 R5515: sout[131] + sout[147] + sout[155] - sout[171] - sout[179]
   + sout[187] + xorout[138] - xorout[139] >= -2
 R5516: sout[131] + sout[147] - sout[155] + sout[171] - sout[179]
   + sout[187] + xorout[138] - xorout[139] >= -2
 R5517: sout[131] - sout[147] + sout[155] + sout[171] - sout[179]
   + sout[187] + xorout[138] - xorout[139] >= -2
 R5518: - sout[131] + sout[147] + sout[155] + sout[171] - sout[179]
   + sout[187] + xorout[138] - xorout[139] >= -2
 R5519: sout[131] + sout[147] - sout[155] - sout[171] + sout[179]
   + sout[187] + xorout[138] - xorout[139] >= -2
 R5520: sout[131] - sout[147] + sout[155] - sout[171] + sout[179]
   + sout[187] + xorout[138] - xorout[139] >= -2
 R5521: - sout[131] + sout[147] + sout[155] - sout[171] + sout[179]
   + sout[187] + xorout[138] - xorout[139] >= -2
 R5522: sout[131] - sout[147] - sout[155] + sout[171] + sout[179]
   + sout[187] + xorout[138] - xorout[139] >= -2
 R5523: - sout[131] + sout[147] - sout[155] + sout[171] + sout[179]
   + sout[187] + xorout[138] - xorout[139] >= -2
 R5524: - sout[131] - sout[147] + sout[155] + sout[171] + sout[179]
   + sout[187] + xorout[138] - xorout[139] >= -2
 R5525: sout[133] + sout[149] + sout[157] + sout[173] + sout[181]
   - sout[189] + xorout[141] >= 0
 R5526: sout[133] + sout[149] + sout[157] + sout[173] - sout[181]
   + sout[189] + xorout[141] >= 0
 R5527: sout[133] + sout[149] + sout[157] - sout[173] + sout[181]
   + sout[189] + xorout[141] >= 0
 R5528: sout[133] + sout[149] - sout[157] + sout[173] + sout[181]
   + sout[189] + xorout[141] >= 0
 R5529: sout[133] - sout[149] + sout[157] + sout[173] + sout[181]
   + sout[189] + xorout[141] >= 0
 R5530: - sout[133] + sout[149] + sout[157] + sout[173] + sout[181]
   + sout[189] + xorout[141] >= 0
 R5531: sout[133] + sout[149] + sout[157] + sout[173] + sout[181]
   + sout[189] - xorout[141] >= 0
 R5532: sout[133] + sout[149] + sout[157] - sout[173] - sout[181]
   - sout[189] + xorout[141] >= -2
 R5533: sout[133] + sout[149] - sout[157] + sout[173] - sout[181]
   - sout[189] + xorout[141] >= -2
 R5534: sout[133] - sout[149] + sout[157] + sout[173] - sout[181]
   - sout[189] + xorout[141] >= -2
 R5535: - sout[133] + sout[149] + sout[157] + sout[173] - sout[181]
   - sout[189] + xorout[141] >= -2
 R5536: sout[133] + sout[149] - sout[157] - sout[173] + sout[181]
   - sout[189] + xorout[141] >= -2
 R5537: sout[133] - sout[149] + sout[157] - sout[173] + sout[181]
   - sout[189] + xorout[141] >= -2
 R5538: - sout[133] + sout[149] + sout[157] - sout[173] + sout[181]
   - sout[189] + xorout[141] >= -2
 R5539: sout[133] - sout[149] - sout[157] + sout[173] + sout[181]
   - sout[189] + xorout[141] >= -2
 R5540: - sout[133] + sout[149] - sout[157] + sout[173] + sout[181]
   - sout[189] + xorout[141] >= -2
 R5541: - sout[133] - sout[149] + sout[157] + sout[173] + sout[181]
   - sout[189] + xorout[141] >= -2
 R5542: sout[133] + sout[149] - sout[157] - sout[173] - sout[181]
   + sout[189] + xorout[141] >= -2
 R5543: sout[133] - sout[149] + sout[157] - sout[173] - sout[181]
   + sout[189] + xorout[141] >= -2
 R5544: - sout[133] + sout[149] + sout[157] - sout[173] - sout[181]
   + sout[189] + xorout[141] >= -2
 R5545: sout[133] - sout[149] - sout[157] + sout[173] - sout[181]
   + sout[189] + xorout[141] >= -2
 R5546: - sout[133] + sout[149] - sout[157] + sout[173] - sout[181]
   + sout[189] + xorout[141] >= -2
 R5547: - sout[133] - sout[149] + sout[157] + sout[173] - sout[181]
   + sout[189] + xorout[141] >= -2
 R5548: sout[133] - sout[149] - sout[157] - sout[173] + sout[181]
   + sout[189] + xorout[141] >= -2
 R5549: - sout[133] + sout[149] - sout[157] - sout[173] + sout[181]
   + sout[189] + xorout[141] >= -2
 R5550: - sout[133] - sout[149] + sout[157] - sout[173] + sout[181]
   + sout[189] + xorout[141] >= -2
 R5551: - sout[133] - sout[149] - sout[157] + sout[173] + sout[181]
   + sout[189] + xorout[141] >= -2
 R5552: sout[133] - sout[149] - sout[157] - sout[173] - sout[181]
   - sout[189] + xorout[141] >= -4
 R5553: - sout[133] + sout[149] - sout[157] - sout[173] - sout[181]
   - sout[189] + xorout[141] >= -4
 R5554: - sout[133] - sout[149] + sout[157] - sout[173] - sout[181]
   - sout[189] + xorout[141] >= -4
 R5555: - sout[133] - sout[149] - sout[157] + sout[173] - sout[181]
   - sout[189] + xorout[141] >= -4
 R5556: - sout[133] - sout[149] - sout[157] - sout[173] + sout[181]
   - sout[189] + xorout[141] >= -4
 R5557: - sout[133] - sout[149] - sout[157] - sout[173] - sout[181]
   + sout[189] + xorout[141] >= -4
 R5558: - sout[188] + xorout[140] >= 0
 R5559: - sout[180] + xorout[140] >= 0
 R5560: - sout[172] + xorout[140] >= 0
 R5561: - sout[156] + xorout[140] >= 0
 R5562: - sout[148] + xorout[140] >= 0
 R5563: - sout[132] + xorout[140] >= 0
 R5564: - sout[132] + sout[133] >= 0
 R5565: - sout[148] + sout[149] >= 0
 R5566: - sout[156] + sout[157] >= 0
 R5567: - sout[172] + sout[173] >= 0
 R5568: - sout[180] + sout[181] >= 0
 R5569: - sout[188] + sout[189] >= 0
 R5570: sout[132] + sout[148] + sout[156] + sout[172] + sout[180]
   + sout[188] - xorout[140] >= 0
 R5571: - xorout[140] + xorout[141] >= 0
 R5572: - sout[133] - sout[149] - sout[157] - sout[173] - sout[181]
   - sout[189] + xorout[140] - xorout[141] >= -6
 R5573: sout[133] + sout[149] - sout[157] - sout[173] - sout[181]
   - sout[189] + xorout[140] - xorout[141] >= -4
 R5574: sout[133] - sout[149] + sout[157] - sout[173] - sout[181]
   - sout[189] + xorout[140] - xorout[141] >= -4
 R5575: - sout[133] + sout[149] + sout[157] - sout[173] - sout[181]
   - sout[189] + xorout[140] - xorout[141] >= -4
 R5576: sout[133] - sout[149] - sout[157] + sout[173] - sout[181]
   - sout[189] + xorout[140] - xorout[141] >= -4
 R5577: - sout[133] + sout[149] - sout[157] + sout[173] - sout[181]
   - sout[189] + xorout[140] - xorout[141] >= -4
 R5578: - sout[133] - sout[149] + sout[157] + sout[173] - sout[181]
   - sout[189] + xorout[140] - xorout[141] >= -4
 R5579: sout[133] - sout[149] - sout[157] - sout[173] + sout[181]
   - sout[189] + xorout[140] - xorout[141] >= -4
 R5580: - sout[133] + sout[149] - sout[157] - sout[173] + sout[181]
   - sout[189] + xorout[140] - xorout[141] >= -4
 R5581: - sout[133] - sout[149] + sout[157] - sout[173] + sout[181]
   - sout[189] + xorout[140] - xorout[141] >= -4
 R5582: - sout[133] - sout[149] - sout[157] + sout[173] + sout[181]
   - sout[189] + xorout[140] - xorout[141] >= -4
 R5583: sout[133] - sout[149] - sout[157] - sout[173] - sout[181]
   + sout[189] + xorout[140] - xorout[141] >= -4
 R5584: - sout[133] + sout[149] - sout[157] - sout[173] - sout[181]
   + sout[189] + xorout[140] - xorout[141] >= -4
 R5585: - sout[133] - sout[149] + sout[157] - sout[173] - sout[181]
   + sout[189] + xorout[140] - xorout[141] >= -4
 R5586: - sout[133] - sout[149] - sout[157] + sout[173] - sout[181]
   + sout[189] + xorout[140] - xorout[141] >= -4
 R5587: - sout[133] - sout[149] - sout[157] - sout[173] + sout[181]
   + sout[189] + xorout[140] - xorout[141] >= -4
 R5588: sout[133] + sout[149] + sout[157] + sout[173] - sout[181]
   - sout[189] + xorout[140] - xorout[141] >= -2
 R5589: sout[133] + sout[149] + sout[157] - sout[173] + sout[181]
   - sout[189] + xorout[140] - xorout[141] >= -2
 R5590: sout[133] + sout[149] - sout[157] + sout[173] + sout[181]
   - sout[189] + xorout[140] - xorout[141] >= -2
 R5591: sout[133] - sout[149] + sout[157] + sout[173] + sout[181]
   - sout[189] + xorout[140] - xorout[141] >= -2
 R5592: - sout[133] + sout[149] + sout[157] + sout[173] + sout[181]
   - sout[189] + xorout[140] - xorout[141] >= -2
 R5593: sout[133] + sout[149] + sout[157] - sout[173] - sout[181]
   + sout[189] + xorout[140] - xorout[141] >= -2
 R5594: sout[133] + sout[149] - sout[157] + sout[173] - sout[181]
   + sout[189] + xorout[140] - xorout[141] >= -2
 R5595: sout[133] - sout[149] + sout[157] + sout[173] - sout[181]
   + sout[189] + xorout[140] - xorout[141] >= -2
 R5596: - sout[133] + sout[149] + sout[157] + sout[173] - sout[181]
   + sout[189] + xorout[140] - xorout[141] >= -2
 R5597: sout[133] + sout[149] - sout[157] - sout[173] + sout[181]
   + sout[189] + xorout[140] - xorout[141] >= -2
 R5598: sout[133] - sout[149] + sout[157] - sout[173] + sout[181]
   + sout[189] + xorout[140] - xorout[141] >= -2
 R5599: - sout[133] + sout[149] + sout[157] - sout[173] + sout[181]
   + sout[189] + xorout[140] - xorout[141] >= -2
 R5600: sout[133] - sout[149] - sout[157] + sout[173] + sout[181]
   + sout[189] + xorout[140] - xorout[141] >= -2
 R5601: - sout[133] + sout[149] - sout[157] + sout[173] + sout[181]
   + sout[189] + xorout[140] - xorout[141] >= -2
 R5602: - sout[133] - sout[149] + sout[157] + sout[173] + sout[181]
   + sout[189] + xorout[140] - xorout[141] >= -2
 R5603: sout[135] + sout[151] + sout[159] + sout[175] + sout[183]
   - sout[191] + xorout[143] >= 0
 R5604: sout[135] + sout[151] + sout[159] + sout[175] - sout[183]
   + sout[191] + xorout[143] >= 0
 R5605: sout[135] + sout[151] + sout[159] - sout[175] + sout[183]
   + sout[191] + xorout[143] >= 0
 R5606: sout[135] + sout[151] - sout[159] + sout[175] + sout[183]
   + sout[191] + xorout[143] >= 0
 R5607: sout[135] - sout[151] + sout[159] + sout[175] + sout[183]
   + sout[191] + xorout[143] >= 0
 R5608: - sout[135] + sout[151] + sout[159] + sout[175] + sout[183]
   + sout[191] + xorout[143] >= 0
 R5609: sout[135] + sout[151] + sout[159] + sout[175] + sout[183]
   + sout[191] - xorout[143] >= 0
 R5610: sout[135] + sout[151] + sout[159] - sout[175] - sout[183]
   - sout[191] + xorout[143] >= -2
 R5611: sout[135] + sout[151] - sout[159] + sout[175] - sout[183]
   - sout[191] + xorout[143] >= -2
 R5612: sout[135] - sout[151] + sout[159] + sout[175] - sout[183]
   - sout[191] + xorout[143] >= -2
 R5613: - sout[135] + sout[151] + sout[159] + sout[175] - sout[183]
   - sout[191] + xorout[143] >= -2
 R5614: sout[135] + sout[151] - sout[159] - sout[175] + sout[183]
   - sout[191] + xorout[143] >= -2
 R5615: sout[135] - sout[151] + sout[159] - sout[175] + sout[183]
   - sout[191] + xorout[143] >= -2
 R5616: - sout[135] + sout[151] + sout[159] - sout[175] + sout[183]
   - sout[191] + xorout[143] >= -2
 R5617: sout[135] - sout[151] - sout[159] + sout[175] + sout[183]
   - sout[191] + xorout[143] >= -2
 R5618: - sout[135] + sout[151] - sout[159] + sout[175] + sout[183]
   - sout[191] + xorout[143] >= -2
 R5619: - sout[135] - sout[151] + sout[159] + sout[175] + sout[183]
   - sout[191] + xorout[143] >= -2
 R5620: sout[135] + sout[151] - sout[159] - sout[175] - sout[183]
   + sout[191] + xorout[143] >= -2
 R5621: sout[135] - sout[151] + sout[159] - sout[175] - sout[183]
   + sout[191] + xorout[143] >= -2
 R5622: - sout[135] + sout[151] + sout[159] - sout[175] - sout[183]
   + sout[191] + xorout[143] >= -2
 R5623: sout[135] - sout[151] - sout[159] + sout[175] - sout[183]
   + sout[191] + xorout[143] >= -2
 R5624: - sout[135] + sout[151] - sout[159] + sout[175] - sout[183]
   + sout[191] + xorout[143] >= -2
 R5625: - sout[135] - sout[151] + sout[159] + sout[175] - sout[183]
   + sout[191] + xorout[143] >= -2
 R5626: sout[135] - sout[151] - sout[159] - sout[175] + sout[183]
   + sout[191] + xorout[143] >= -2
 R5627: - sout[135] + sout[151] - sout[159] - sout[175] + sout[183]
   + sout[191] + xorout[143] >= -2
 R5628: - sout[135] - sout[151] + sout[159] - sout[175] + sout[183]
   + sout[191] + xorout[143] >= -2
 R5629: - sout[135] - sout[151] - sout[159] + sout[175] + sout[183]
   + sout[191] + xorout[143] >= -2
 R5630: sout[135] - sout[151] - sout[159] - sout[175] - sout[183]
   - sout[191] + xorout[143] >= -4
 R5631: - sout[135] + sout[151] - sout[159] - sout[175] - sout[183]
   - sout[191] + xorout[143] >= -4
 R5632: - sout[135] - sout[151] + sout[159] - sout[175] - sout[183]
   - sout[191] + xorout[143] >= -4
 R5633: - sout[135] - sout[151] - sout[159] + sout[175] - sout[183]
   - sout[191] + xorout[143] >= -4
 R5634: - sout[135] - sout[151] - sout[159] - sout[175] + sout[183]
   - sout[191] + xorout[143] >= -4
 R5635: - sout[135] - sout[151] - sout[159] - sout[175] - sout[183]
   + sout[191] + xorout[143] >= -4
 R5636: - sout[190] + xorout[142] >= 0
 R5637: - sout[182] + xorout[142] >= 0
 R5638: - sout[174] + xorout[142] >= 0
 R5639: - sout[158] + xorout[142] >= 0
 R5640: - sout[150] + xorout[142] >= 0
 R5641: - sout[134] + xorout[142] >= 0
 R5642: - sout[134] + sout[135] >= 0
 R5643: - sout[150] + sout[151] >= 0
 R5644: - sout[158] + sout[159] >= 0
 R5645: - sout[174] + sout[175] >= 0
 R5646: - sout[182] + sout[183] >= 0
 R5647: - sout[190] + sout[191] >= 0
 R5648: sout[134] + sout[150] + sout[158] + sout[174] + sout[182]
   + sout[190] - xorout[142] >= 0
 R5649: - xorout[142] + xorout[143] >= 0
 R5650: - sout[135] - sout[151] - sout[159] - sout[175] - sout[183]
   - sout[191] + xorout[142] - xorout[143] >= -6
 R5651: sout[135] + sout[151] - sout[159] - sout[175] - sout[183]
   - sout[191] + xorout[142] - xorout[143] >= -4
 R5652: sout[135] - sout[151] + sout[159] - sout[175] - sout[183]
   - sout[191] + xorout[142] - xorout[143] >= -4
 R5653: - sout[135] + sout[151] + sout[159] - sout[175] - sout[183]
   - sout[191] + xorout[142] - xorout[143] >= -4
 R5654: sout[135] - sout[151] - sout[159] + sout[175] - sout[183]
   - sout[191] + xorout[142] - xorout[143] >= -4
 R5655: - sout[135] + sout[151] - sout[159] + sout[175] - sout[183]
   - sout[191] + xorout[142] - xorout[143] >= -4
 R5656: - sout[135] - sout[151] + sout[159] + sout[175] - sout[183]
   - sout[191] + xorout[142] - xorout[143] >= -4
 R5657: sout[135] - sout[151] - sout[159] - sout[175] + sout[183]
   - sout[191] + xorout[142] - xorout[143] >= -4
 R5658: - sout[135] + sout[151] - sout[159] - sout[175] + sout[183]
   - sout[191] + xorout[142] - xorout[143] >= -4
 R5659: - sout[135] - sout[151] + sout[159] - sout[175] + sout[183]
   - sout[191] + xorout[142] - xorout[143] >= -4
 R5660: - sout[135] - sout[151] - sout[159] + sout[175] + sout[183]
   - sout[191] + xorout[142] - xorout[143] >= -4
 R5661: sout[135] - sout[151] - sout[159] - sout[175] - sout[183]
   + sout[191] + xorout[142] - xorout[143] >= -4
 R5662: - sout[135] + sout[151] - sout[159] - sout[175] - sout[183]
   + sout[191] + xorout[142] - xorout[143] >= -4
 R5663: - sout[135] - sout[151] + sout[159] - sout[175] - sout[183]
   + sout[191] + xorout[142] - xorout[143] >= -4
 R5664: - sout[135] - sout[151] - sout[159] + sout[175] - sout[183]
   + sout[191] + xorout[142] - xorout[143] >= -4
 R5665: - sout[135] - sout[151] - sout[159] - sout[175] + sout[183]
   + sout[191] + xorout[142] - xorout[143] >= -4
 R5666: sout[135] + sout[151] + sout[159] + sout[175] - sout[183]
   - sout[191] + xorout[142] - xorout[143] >= -2
 R5667: sout[135] + sout[151] + sout[159] - sout[175] + sout[183]
   - sout[191] + xorout[142] - xorout[143] >= -2
 R5668: sout[135] + sout[151] - sout[159] + sout[175] + sout[183]
   - sout[191] + xorout[142] - xorout[143] >= -2
 R5669: sout[135] - sout[151] + sout[159] + sout[175] + sout[183]
   - sout[191] + xorout[142] - xorout[143] >= -2
 R5670: - sout[135] + sout[151] + sout[159] + sout[175] + sout[183]
   - sout[191] + xorout[142] - xorout[143] >= -2
 R5671: sout[135] + sout[151] + sout[159] - sout[175] - sout[183]
   + sout[191] + xorout[142] - xorout[143] >= -2
 R5672: sout[135] + sout[151] - sout[159] + sout[175] - sout[183]
   + sout[191] + xorout[142] - xorout[143] >= -2
 R5673: sout[135] - sout[151] + sout[159] + sout[175] - sout[183]
   + sout[191] + xorout[142] - xorout[143] >= -2
 R5674: - sout[135] + sout[151] + sout[159] + sout[175] - sout[183]
   + sout[191] + xorout[142] - xorout[143] >= -2
 R5675: sout[135] + sout[151] - sout[159] - sout[175] + sout[183]
   + sout[191] + xorout[142] - xorout[143] >= -2
 R5676: sout[135] - sout[151] + sout[159] - sout[175] + sout[183]
   + sout[191] + xorout[142] - xorout[143] >= -2
 R5677: - sout[135] + sout[151] + sout[159] - sout[175] + sout[183]
   + sout[191] + xorout[142] - xorout[143] >= -2
 R5678: sout[135] - sout[151] - sout[159] + sout[175] + sout[183]
   + sout[191] + xorout[142] - xorout[143] >= -2
 R5679: - sout[135] + sout[151] - sout[159] + sout[175] + sout[183]
   + sout[191] + xorout[142] - xorout[143] >= -2
 R5680: - sout[135] - sout[151] + sout[159] + sout[175] + sout[183]
   + sout[191] + xorout[142] - xorout[143] >= -2
 R5681: sout[129] + sout[137] + sout[153] + sout[161] + sout[177]
   - sout[185] + xorout[145] >= 0
 R5682: sout[129] + sout[137] + sout[153] + sout[161] - sout[177]
   + sout[185] + xorout[145] >= 0
 R5683: sout[129] + sout[137] + sout[153] - sout[161] + sout[177]
   + sout[185] + xorout[145] >= 0
 R5684: sout[129] + sout[137] - sout[153] + sout[161] + sout[177]
   + sout[185] + xorout[145] >= 0
 R5685: sout[129] - sout[137] + sout[153] + sout[161] + sout[177]
   + sout[185] + xorout[145] >= 0
 R5686: - sout[129] + sout[137] + sout[153] + sout[161] + sout[177]
   + sout[185] + xorout[145] >= 0
 R5687: sout[129] + sout[137] + sout[153] + sout[161] + sout[177]
   + sout[185] - xorout[145] >= 0
 R5688: sout[129] + sout[137] + sout[153] - sout[161] - sout[177]
   - sout[185] + xorout[145] >= -2
 R5689: sout[129] + sout[137] - sout[153] + sout[161] - sout[177]
   - sout[185] + xorout[145] >= -2
 R5690: sout[129] - sout[137] + sout[153] + sout[161] - sout[177]
   - sout[185] + xorout[145] >= -2
 R5691: - sout[129] + sout[137] + sout[153] + sout[161] - sout[177]
   - sout[185] + xorout[145] >= -2
 R5692: sout[129] + sout[137] - sout[153] - sout[161] + sout[177]
   - sout[185] + xorout[145] >= -2
 R5693: sout[129] - sout[137] + sout[153] - sout[161] + sout[177]
   - sout[185] + xorout[145] >= -2
 R5694: - sout[129] + sout[137] + sout[153] - sout[161] + sout[177]
   - sout[185] + xorout[145] >= -2
 R5695: sout[129] - sout[137] - sout[153] + sout[161] + sout[177]
   - sout[185] + xorout[145] >= -2
 R5696: - sout[129] + sout[137] - sout[153] + sout[161] + sout[177]
   - sout[185] + xorout[145] >= -2
 R5697: - sout[129] - sout[137] + sout[153] + sout[161] + sout[177]
   - sout[185] + xorout[145] >= -2
 R5698: sout[129] + sout[137] - sout[153] - sout[161] - sout[177]
   + sout[185] + xorout[145] >= -2
 R5699: sout[129] - sout[137] + sout[153] - sout[161] - sout[177]
   + sout[185] + xorout[145] >= -2
 R5700: - sout[129] + sout[137] + sout[153] - sout[161] - sout[177]
   + sout[185] + xorout[145] >= -2
 R5701: sout[129] - sout[137] - sout[153] + sout[161] - sout[177]
   + sout[185] + xorout[145] >= -2
 R5702: - sout[129] + sout[137] - sout[153] + sout[161] - sout[177]
   + sout[185] + xorout[145] >= -2
 R5703: - sout[129] - sout[137] + sout[153] + sout[161] - sout[177]
   + sout[185] + xorout[145] >= -2
 R5704: sout[129] - sout[137] - sout[153] - sout[161] + sout[177]
   + sout[185] + xorout[145] >= -2
 R5705: - sout[129] + sout[137] - sout[153] - sout[161] + sout[177]
   + sout[185] + xorout[145] >= -2
 R5706: - sout[129] - sout[137] + sout[153] - sout[161] + sout[177]
   + sout[185] + xorout[145] >= -2
 R5707: - sout[129] - sout[137] - sout[153] + sout[161] + sout[177]
   + sout[185] + xorout[145] >= -2
 R5708: sout[129] - sout[137] - sout[153] - sout[161] - sout[177]
   - sout[185] + xorout[145] >= -4
 R5709: - sout[129] + sout[137] - sout[153] - sout[161] - sout[177]
   - sout[185] + xorout[145] >= -4
 R5710: - sout[129] - sout[137] + sout[153] - sout[161] - sout[177]
   - sout[185] + xorout[145] >= -4
 R5711: - sout[129] - sout[137] - sout[153] + sout[161] - sout[177]
   - sout[185] + xorout[145] >= -4
 R5712: - sout[129] - sout[137] - sout[153] - sout[161] + sout[177]
   - sout[185] + xorout[145] >= -4
 R5713: - sout[129] - sout[137] - sout[153] - sout[161] - sout[177]
   + sout[185] + xorout[145] >= -4
 R5714: - sout[184] + xorout[144] >= 0
 R5715: - sout[176] + xorout[144] >= 0
 R5716: - sout[160] + xorout[144] >= 0
 R5717: - sout[152] + xorout[144] >= 0
 R5718: - sout[136] + xorout[144] >= 0
 R5719: - sout[128] + xorout[144] >= 0
 R5720: - sout[128] + sout[129] >= 0
 R5721: - sout[136] + sout[137] >= 0
 R5722: - sout[152] + sout[153] >= 0
 R5723: - sout[160] + sout[161] >= 0
 R5724: - sout[176] + sout[177] >= 0
 R5725: - sout[184] + sout[185] >= 0
 R5726: sout[128] + sout[136] + sout[152] + sout[160] + sout[176]
   + sout[184] - xorout[144] >= 0
 R5727: - xorout[144] + xorout[145] >= 0
 R5728: - sout[129] - sout[137] - sout[153] - sout[161] - sout[177]
   - sout[185] + xorout[144] - xorout[145] >= -6
 R5729: sout[129] + sout[137] - sout[153] - sout[161] - sout[177]
   - sout[185] + xorout[144] - xorout[145] >= -4
 R5730: sout[129] - sout[137] + sout[153] - sout[161] - sout[177]
   - sout[185] + xorout[144] - xorout[145] >= -4
 R5731: - sout[129] + sout[137] + sout[153] - sout[161] - sout[177]
   - sout[185] + xorout[144] - xorout[145] >= -4
 R5732: sout[129] - sout[137] - sout[153] + sout[161] - sout[177]
   - sout[185] + xorout[144] - xorout[145] >= -4
 R5733: - sout[129] + sout[137] - sout[153] + sout[161] - sout[177]
   - sout[185] + xorout[144] - xorout[145] >= -4
 R5734: - sout[129] - sout[137] + sout[153] + sout[161] - sout[177]
   - sout[185] + xorout[144] - xorout[145] >= -4
 R5735: sout[129] - sout[137] - sout[153] - sout[161] + sout[177]
   - sout[185] + xorout[144] - xorout[145] >= -4
 R5736: - sout[129] + sout[137] - sout[153] - sout[161] + sout[177]
   - sout[185] + xorout[144] - xorout[145] >= -4
 R5737: - sout[129] - sout[137] + sout[153] - sout[161] + sout[177]
   - sout[185] + xorout[144] - xorout[145] >= -4
 R5738: - sout[129] - sout[137] - sout[153] + sout[161] + sout[177]
   - sout[185] + xorout[144] - xorout[145] >= -4
 R5739: sout[129] - sout[137] - sout[153] - sout[161] - sout[177]
   + sout[185] + xorout[144] - xorout[145] >= -4
 R5740: - sout[129] + sout[137] - sout[153] - sout[161] - sout[177]
   + sout[185] + xorout[144] - xorout[145] >= -4
 R5741: - sout[129] - sout[137] + sout[153] - sout[161] - sout[177]
   + sout[185] + xorout[144] - xorout[145] >= -4
 R5742: - sout[129] - sout[137] - sout[153] + sout[161] - sout[177]
   + sout[185] + xorout[144] - xorout[145] >= -4
 R5743: - sout[129] - sout[137] - sout[153] - sout[161] + sout[177]
   + sout[185] + xorout[144] - xorout[145] >= -4
 R5744: sout[129] + sout[137] + sout[153] + sout[161] - sout[177]
   - sout[185] + xorout[144] - xorout[145] >= -2
 R5745: sout[129] + sout[137] + sout[153] - sout[161] + sout[177]
   - sout[185] + xorout[144] - xorout[145] >= -2
 R5746: sout[129] + sout[137] - sout[153] + sout[161] + sout[177]
   - sout[185] + xorout[144] - xorout[145] >= -2
 R5747: sout[129] - sout[137] + sout[153] + sout[161] + sout[177]
   - sout[185] + xorout[144] - xorout[145] >= -2
 R5748: - sout[129] + sout[137] + sout[153] + sout[161] + sout[177]
   - sout[185] + xorout[144] - xorout[145] >= -2
 R5749: sout[129] + sout[137] + sout[153] - sout[161] - sout[177]
   + sout[185] + xorout[144] - xorout[145] >= -2
 R5750: sout[129] + sout[137] - sout[153] + sout[161] - sout[177]
   + sout[185] + xorout[144] - xorout[145] >= -2
 R5751: sout[129] - sout[137] + sout[153] + sout[161] - sout[177]
   + sout[185] + xorout[144] - xorout[145] >= -2
 R5752: - sout[129] + sout[137] + sout[153] + sout[161] - sout[177]
   + sout[185] + xorout[144] - xorout[145] >= -2
 R5753: sout[129] + sout[137] - sout[153] - sout[161] + sout[177]
   + sout[185] + xorout[144] - xorout[145] >= -2
 R5754: sout[129] - sout[137] + sout[153] - sout[161] + sout[177]
   + sout[185] + xorout[144] - xorout[145] >= -2
 R5755: - sout[129] + sout[137] + sout[153] - sout[161] + sout[177]
   + sout[185] + xorout[144] - xorout[145] >= -2
 R5756: sout[129] - sout[137] - sout[153] + sout[161] + sout[177]
   + sout[185] + xorout[144] - xorout[145] >= -2
 R5757: - sout[129] + sout[137] - sout[153] + sout[161] + sout[177]
   + sout[185] + xorout[144] - xorout[145] >= -2
 R5758: - sout[129] - sout[137] + sout[153] + sout[161] + sout[177]
   + sout[185] + xorout[144] - xorout[145] >= -2
 R5759: sout[131] + sout[139] + sout[155] + sout[163] + sout[179]
   - sout[187] + xorout[147] >= 0
 R5760: sout[131] + sout[139] + sout[155] + sout[163] - sout[179]
   + sout[187] + xorout[147] >= 0
 R5761: sout[131] + sout[139] + sout[155] - sout[163] + sout[179]
   + sout[187] + xorout[147] >= 0
 R5762: sout[131] + sout[139] - sout[155] + sout[163] + sout[179]
   + sout[187] + xorout[147] >= 0
 R5763: sout[131] - sout[139] + sout[155] + sout[163] + sout[179]
   + sout[187] + xorout[147] >= 0
 R5764: - sout[131] + sout[139] + sout[155] + sout[163] + sout[179]
   + sout[187] + xorout[147] >= 0
 R5765: sout[131] + sout[139] + sout[155] + sout[163] + sout[179]
   + sout[187] - xorout[147] >= 0
 R5766: sout[131] + sout[139] + sout[155] - sout[163] - sout[179]
   - sout[187] + xorout[147] >= -2
 R5767: sout[131] + sout[139] - sout[155] + sout[163] - sout[179]
   - sout[187] + xorout[147] >= -2
 R5768: sout[131] - sout[139] + sout[155] + sout[163] - sout[179]
   - sout[187] + xorout[147] >= -2
 R5769: - sout[131] + sout[139] + sout[155] + sout[163] - sout[179]
   - sout[187] + xorout[147] >= -2
 R5770: sout[131] + sout[139] - sout[155] - sout[163] + sout[179]
   - sout[187] + xorout[147] >= -2
 R5771: sout[131] - sout[139] + sout[155] - sout[163] + sout[179]
   - sout[187] + xorout[147] >= -2
 R5772: - sout[131] + sout[139] + sout[155] - sout[163] + sout[179]
   - sout[187] + xorout[147] >= -2
 R5773: sout[131] - sout[139] - sout[155] + sout[163] + sout[179]
   - sout[187] + xorout[147] >= -2
 R5774: - sout[131] + sout[139] - sout[155] + sout[163] + sout[179]
   - sout[187] + xorout[147] >= -2
 R5775: - sout[131] - sout[139] + sout[155] + sout[163] + sout[179]
   - sout[187] + xorout[147] >= -2
 R5776: sout[131] + sout[139] - sout[155] - sout[163] - sout[179]
   + sout[187] + xorout[147] >= -2
 R5777: sout[131] - sout[139] + sout[155] - sout[163] - sout[179]
   + sout[187] + xorout[147] >= -2
 R5778: - sout[131] + sout[139] + sout[155] - sout[163] - sout[179]
   + sout[187] + xorout[147] >= -2
 R5779: sout[131] - sout[139] - sout[155] + sout[163] - sout[179]
   + sout[187] + xorout[147] >= -2
 R5780: - sout[131] + sout[139] - sout[155] + sout[163] - sout[179]
   + sout[187] + xorout[147] >= -2
 R5781: - sout[131] - sout[139] + sout[155] + sout[163] - sout[179]
   + sout[187] + xorout[147] >= -2
 R5782: sout[131] - sout[139] - sout[155] - sout[163] + sout[179]
   + sout[187] + xorout[147] >= -2
 R5783: - sout[131] + sout[139] - sout[155] - sout[163] + sout[179]
   + sout[187] + xorout[147] >= -2
 R5784: - sout[131] - sout[139] + sout[155] - sout[163] + sout[179]
   + sout[187] + xorout[147] >= -2
 R5785: - sout[131] - sout[139] - sout[155] + sout[163] + sout[179]
   + sout[187] + xorout[147] >= -2
 R5786: sout[131] - sout[139] - sout[155] - sout[163] - sout[179]
   - sout[187] + xorout[147] >= -4
 R5787: - sout[131] + sout[139] - sout[155] - sout[163] - sout[179]
   - sout[187] + xorout[147] >= -4
 R5788: - sout[131] - sout[139] + sout[155] - sout[163] - sout[179]
   - sout[187] + xorout[147] >= -4
 R5789: - sout[131] - sout[139] - sout[155] + sout[163] - sout[179]
   - sout[187] + xorout[147] >= -4
 R5790: - sout[131] - sout[139] - sout[155] - sout[163] + sout[179]
   - sout[187] + xorout[147] >= -4
 R5791: - sout[131] - sout[139] - sout[155] - sout[163] - sout[179]
   + sout[187] + xorout[147] >= -4
 R5792: - sout[186] + xorout[146] >= 0
 R5793: - sout[178] + xorout[146] >= 0
 R5794: - sout[162] + xorout[146] >= 0
 R5795: - sout[154] + xorout[146] >= 0
 R5796: - sout[138] + xorout[146] >= 0
 R5797: - sout[130] + xorout[146] >= 0
 R5798: - sout[130] + sout[131] >= 0
 R5799: - sout[138] + sout[139] >= 0
 R5800: - sout[154] + sout[155] >= 0
 R5801: - sout[162] + sout[163] >= 0
 R5802: - sout[178] + sout[179] >= 0
 R5803: - sout[186] + sout[187] >= 0
 R5804: sout[130] + sout[138] + sout[154] + sout[162] + sout[178]
   + sout[186] - xorout[146] >= 0
 R5805: - xorout[146] + xorout[147] >= 0
 R5806: - sout[131] - sout[139] - sout[155] - sout[163] - sout[179]
   - sout[187] + xorout[146] - xorout[147] >= -6
 R5807: sout[131] + sout[139] - sout[155] - sout[163] - sout[179]
   - sout[187] + xorout[146] - xorout[147] >= -4
 R5808: sout[131] - sout[139] + sout[155] - sout[163] - sout[179]
   - sout[187] + xorout[146] - xorout[147] >= -4
 R5809: - sout[131] + sout[139] + sout[155] - sout[163] - sout[179]
   - sout[187] + xorout[146] - xorout[147] >= -4
 R5810: sout[131] - sout[139] - sout[155] + sout[163] - sout[179]
   - sout[187] + xorout[146] - xorout[147] >= -4
 R5811: - sout[131] + sout[139] - sout[155] + sout[163] - sout[179]
   - sout[187] + xorout[146] - xorout[147] >= -4
 R5812: - sout[131] - sout[139] + sout[155] + sout[163] - sout[179]
   - sout[187] + xorout[146] - xorout[147] >= -4
 R5813: sout[131] - sout[139] - sout[155] - sout[163] + sout[179]
   - sout[187] + xorout[146] - xorout[147] >= -4
 R5814: - sout[131] + sout[139] - sout[155] - sout[163] + sout[179]
   - sout[187] + xorout[146] - xorout[147] >= -4
 R5815: - sout[131] - sout[139] + sout[155] - sout[163] + sout[179]
   - sout[187] + xorout[146] - xorout[147] >= -4
 R5816: - sout[131] - sout[139] - sout[155] + sout[163] + sout[179]
   - sout[187] + xorout[146] - xorout[147] >= -4
 R5817: sout[131] - sout[139] - sout[155] - sout[163] - sout[179]
   + sout[187] + xorout[146] - xorout[147] >= -4
 R5818: - sout[131] + sout[139] - sout[155] - sout[163] - sout[179]
   + sout[187] + xorout[146] - xorout[147] >= -4
 R5819: - sout[131] - sout[139] + sout[155] - sout[163] - sout[179]
   + sout[187] + xorout[146] - xorout[147] >= -4
 R5820: - sout[131] - sout[139] - sout[155] + sout[163] - sout[179]
   + sout[187] + xorout[146] - xorout[147] >= -4
 R5821: - sout[131] - sout[139] - sout[155] - sout[163] + sout[179]
   + sout[187] + xorout[146] - xorout[147] >= -4
 R5822: sout[131] + sout[139] + sout[155] + sout[163] - sout[179]
   - sout[187] + xorout[146] - xorout[147] >= -2
 R5823: sout[131] + sout[139] + sout[155] - sout[163] + sout[179]
   - sout[187] + xorout[146] - xorout[147] >= -2
 R5824: sout[131] + sout[139] - sout[155] + sout[163] + sout[179]
   - sout[187] + xorout[146] - xorout[147] >= -2
 R5825: sout[131] - sout[139] + sout[155] + sout[163] + sout[179]
   - sout[187] + xorout[146] - xorout[147] >= -2
 R5826: - sout[131] + sout[139] + sout[155] + sout[163] + sout[179]
   - sout[187] + xorout[146] - xorout[147] >= -2
 R5827: sout[131] + sout[139] + sout[155] - sout[163] - sout[179]
   + sout[187] + xorout[146] - xorout[147] >= -2
 R5828: sout[131] + sout[139] - sout[155] + sout[163] - sout[179]
   + sout[187] + xorout[146] - xorout[147] >= -2
 R5829: sout[131] - sout[139] + sout[155] + sout[163] - sout[179]
   + sout[187] + xorout[146] - xorout[147] >= -2
 R5830: - sout[131] + sout[139] + sout[155] + sout[163] - sout[179]
   + sout[187] + xorout[146] - xorout[147] >= -2
 R5831: sout[131] + sout[139] - sout[155] - sout[163] + sout[179]
   + sout[187] + xorout[146] - xorout[147] >= -2
 R5832: sout[131] - sout[139] + sout[155] - sout[163] + sout[179]
   + sout[187] + xorout[146] - xorout[147] >= -2
 R5833: - sout[131] + sout[139] + sout[155] - sout[163] + sout[179]
   + sout[187] + xorout[146] - xorout[147] >= -2
 R5834: sout[131] - sout[139] - sout[155] + sout[163] + sout[179]
   + sout[187] + xorout[146] - xorout[147] >= -2
 R5835: - sout[131] + sout[139] - sout[155] + sout[163] + sout[179]
   + sout[187] + xorout[146] - xorout[147] >= -2
 R5836: - sout[131] - sout[139] + sout[155] + sout[163] + sout[179]
   + sout[187] + xorout[146] - xorout[147] >= -2
 R5837: sout[133] + sout[141] + sout[157] + sout[165] + sout[181]
   - sout[189] + xorout[149] >= 0
 R5838: sout[133] + sout[141] + sout[157] + sout[165] - sout[181]
   + sout[189] + xorout[149] >= 0
 R5839: sout[133] + sout[141] + sout[157] - sout[165] + sout[181]
   + sout[189] + xorout[149] >= 0
 R5840: sout[133] + sout[141] - sout[157] + sout[165] + sout[181]
   + sout[189] + xorout[149] >= 0
 R5841: sout[133] - sout[141] + sout[157] + sout[165] + sout[181]
   + sout[189] + xorout[149] >= 0
 R5842: - sout[133] + sout[141] + sout[157] + sout[165] + sout[181]
   + sout[189] + xorout[149] >= 0
 R5843: sout[133] + sout[141] + sout[157] + sout[165] + sout[181]
   + sout[189] - xorout[149] >= 0
 R5844: sout[133] + sout[141] + sout[157] - sout[165] - sout[181]
   - sout[189] + xorout[149] >= -2
 R5845: sout[133] + sout[141] - sout[157] + sout[165] - sout[181]
   - sout[189] + xorout[149] >= -2
 R5846: sout[133] - sout[141] + sout[157] + sout[165] - sout[181]
   - sout[189] + xorout[149] >= -2
 R5847: - sout[133] + sout[141] + sout[157] + sout[165] - sout[181]
   - sout[189] + xorout[149] >= -2
 R5848: sout[133] + sout[141] - sout[157] - sout[165] + sout[181]
   - sout[189] + xorout[149] >= -2
 R5849: sout[133] - sout[141] + sout[157] - sout[165] + sout[181]
   - sout[189] + xorout[149] >= -2
 R5850: - sout[133] + sout[141] + sout[157] - sout[165] + sout[181]
   - sout[189] + xorout[149] >= -2
 R5851: sout[133] - sout[141] - sout[157] + sout[165] + sout[181]
   - sout[189] + xorout[149] >= -2
 R5852: - sout[133] + sout[141] - sout[157] + sout[165] + sout[181]
   - sout[189] + xorout[149] >= -2
 R5853: - sout[133] - sout[141] + sout[157] + sout[165] + sout[181]
   - sout[189] + xorout[149] >= -2
 R5854: sout[133] + sout[141] - sout[157] - sout[165] - sout[181]
   + sout[189] + xorout[149] >= -2
 R5855: sout[133] - sout[141] + sout[157] - sout[165] - sout[181]
   + sout[189] + xorout[149] >= -2
 R5856: - sout[133] + sout[141] + sout[157] - sout[165] - sout[181]
   + sout[189] + xorout[149] >= -2
 R5857: sout[133] - sout[141] - sout[157] + sout[165] - sout[181]
   + sout[189] + xorout[149] >= -2
 R5858: - sout[133] + sout[141] - sout[157] + sout[165] - sout[181]
   + sout[189] + xorout[149] >= -2
 R5859: - sout[133] - sout[141] + sout[157] + sout[165] - sout[181]
   + sout[189] + xorout[149] >= -2
 R5860: sout[133] - sout[141] - sout[157] - sout[165] + sout[181]
   + sout[189] + xorout[149] >= -2
 R5861: - sout[133] + sout[141] - sout[157] - sout[165] + sout[181]
   + sout[189] + xorout[149] >= -2
 R5862: - sout[133] - sout[141] + sout[157] - sout[165] + sout[181]
   + sout[189] + xorout[149] >= -2
 R5863: - sout[133] - sout[141] - sout[157] + sout[165] + sout[181]
   + sout[189] + xorout[149] >= -2
 R5864: sout[133] - sout[141] - sout[157] - sout[165] - sout[181]
   - sout[189] + xorout[149] >= -4
 R5865: - sout[133] + sout[141] - sout[157] - sout[165] - sout[181]
   - sout[189] + xorout[149] >= -4
 R5866: - sout[133] - sout[141] + sout[157] - sout[165] - sout[181]
   - sout[189] + xorout[149] >= -4
 R5867: - sout[133] - sout[141] - sout[157] + sout[165] - sout[181]
   - sout[189] + xorout[149] >= -4
 R5868: - sout[133] - sout[141] - sout[157] - sout[165] + sout[181]
   - sout[189] + xorout[149] >= -4
 R5869: - sout[133] - sout[141] - sout[157] - sout[165] - sout[181]
   + sout[189] + xorout[149] >= -4
 R5870: - sout[188] + xorout[148] >= 0
 R5871: - sout[180] + xorout[148] >= 0
 R5872: - sout[164] + xorout[148] >= 0
 R5873: - sout[156] + xorout[148] >= 0
 R5874: - sout[140] + xorout[148] >= 0
 R5875: - sout[132] + xorout[148] >= 0
 R5876: - sout[132] + sout[133] >= 0
 R5877: - sout[140] + sout[141] >= 0
 R5878: - sout[156] + sout[157] >= 0
 R5879: - sout[164] + sout[165] >= 0
 R5880: - sout[180] + sout[181] >= 0
 R5881: - sout[188] + sout[189] >= 0
 R5882: sout[132] + sout[140] + sout[156] + sout[164] + sout[180]
   + sout[188] - xorout[148] >= 0
 R5883: - xorout[148] + xorout[149] >= 0
 R5884: - sout[133] - sout[141] - sout[157] - sout[165] - sout[181]
   - sout[189] + xorout[148] - xorout[149] >= -6
 R5885: sout[133] + sout[141] - sout[157] - sout[165] - sout[181]
   - sout[189] + xorout[148] - xorout[149] >= -4
 R5886: sout[133] - sout[141] + sout[157] - sout[165] - sout[181]
   - sout[189] + xorout[148] - xorout[149] >= -4
 R5887: - sout[133] + sout[141] + sout[157] - sout[165] - sout[181]
   - sout[189] + xorout[148] - xorout[149] >= -4
 R5888: sout[133] - sout[141] - sout[157] + sout[165] - sout[181]
   - sout[189] + xorout[148] - xorout[149] >= -4
 R5889: - sout[133] + sout[141] - sout[157] + sout[165] - sout[181]
   - sout[189] + xorout[148] - xorout[149] >= -4
 R5890: - sout[133] - sout[141] + sout[157] + sout[165] - sout[181]
   - sout[189] + xorout[148] - xorout[149] >= -4
 R5891: sout[133] - sout[141] - sout[157] - sout[165] + sout[181]
   - sout[189] + xorout[148] - xorout[149] >= -4
 R5892: - sout[133] + sout[141] - sout[157] - sout[165] + sout[181]
   - sout[189] + xorout[148] - xorout[149] >= -4
 R5893: - sout[133] - sout[141] + sout[157] - sout[165] + sout[181]
   - sout[189] + xorout[148] - xorout[149] >= -4
 R5894: - sout[133] - sout[141] - sout[157] + sout[165] + sout[181]
   - sout[189] + xorout[148] - xorout[149] >= -4
 R5895: sout[133] - sout[141] - sout[157] - sout[165] - sout[181]
   + sout[189] + xorout[148] - xorout[149] >= -4
 R5896: - sout[133] + sout[141] - sout[157] - sout[165] - sout[181]
   + sout[189] + xorout[148] - xorout[149] >= -4
 R5897: - sout[133] - sout[141] + sout[157] - sout[165] - sout[181]
   + sout[189] + xorout[148] - xorout[149] >= -4
 R5898: - sout[133] - sout[141] - sout[157] + sout[165] - sout[181]
   + sout[189] + xorout[148] - xorout[149] >= -4
 R5899: - sout[133] - sout[141] - sout[157] - sout[165] + sout[181]
   + sout[189] + xorout[148] - xorout[149] >= -4
 R5900: sout[133] + sout[141] + sout[157] + sout[165] - sout[181]
   - sout[189] + xorout[148] - xorout[149] >= -2
 R5901: sout[133] + sout[141] + sout[157] - sout[165] + sout[181]
   - sout[189] + xorout[148] - xorout[149] >= -2
 R5902: sout[133] + sout[141] - sout[157] + sout[165] + sout[181]
   - sout[189] + xorout[148] - xorout[149] >= -2
 R5903: sout[133] - sout[141] + sout[157] + sout[165] + sout[181]
   - sout[189] + xorout[148] - xorout[149] >= -2
 R5904: - sout[133] + sout[141] + sout[157] + sout[165] + sout[181]
   - sout[189] + xorout[148] - xorout[149] >= -2
 R5905: sout[133] + sout[141] + sout[157] - sout[165] - sout[181]
   + sout[189] + xorout[148] - xorout[149] >= -2
 R5906: sout[133] + sout[141] - sout[157] + sout[165] - sout[181]
   + sout[189] + xorout[148] - xorout[149] >= -2
 R5907: sout[133] - sout[141] + sout[157] + sout[165] - sout[181]
   + sout[189] + xorout[148] - xorout[149] >= -2
 R5908: - sout[133] + sout[141] + sout[157] + sout[165] - sout[181]
   + sout[189] + xorout[148] - xorout[149] >= -2
 R5909: sout[133] + sout[141] - sout[157] - sout[165] + sout[181]
   + sout[189] + xorout[148] - xorout[149] >= -2
 R5910: sout[133] - sout[141] + sout[157] - sout[165] + sout[181]
   + sout[189] + xorout[148] - xorout[149] >= -2
 R5911: - sout[133] + sout[141] + sout[157] - sout[165] + sout[181]
   + sout[189] + xorout[148] - xorout[149] >= -2
 R5912: sout[133] - sout[141] - sout[157] + sout[165] + sout[181]
   + sout[189] + xorout[148] - xorout[149] >= -2
 R5913: - sout[133] + sout[141] - sout[157] + sout[165] + sout[181]
   + sout[189] + xorout[148] - xorout[149] >= -2
 R5914: - sout[133] - sout[141] + sout[157] + sout[165] + sout[181]
   + sout[189] + xorout[148] - xorout[149] >= -2
 R5915: sout[135] + sout[143] + sout[159] + sout[167] + sout[183]
   - sout[191] + xorout[151] >= 0
 R5916: sout[135] + sout[143] + sout[159] + sout[167] - sout[183]
   + sout[191] + xorout[151] >= 0
 R5917: sout[135] + sout[143] + sout[159] - sout[167] + sout[183]
   + sout[191] + xorout[151] >= 0
 R5918: sout[135] + sout[143] - sout[159] + sout[167] + sout[183]
   + sout[191] + xorout[151] >= 0
 R5919: sout[135] - sout[143] + sout[159] + sout[167] + sout[183]
   + sout[191] + xorout[151] >= 0
 R5920: - sout[135] + sout[143] + sout[159] + sout[167] + sout[183]
   + sout[191] + xorout[151] >= 0
 R5921: sout[135] + sout[143] + sout[159] + sout[167] + sout[183]
   + sout[191] - xorout[151] >= 0
 R5922: sout[135] + sout[143] + sout[159] - sout[167] - sout[183]
   - sout[191] + xorout[151] >= -2
 R5923: sout[135] + sout[143] - sout[159] + sout[167] - sout[183]
   - sout[191] + xorout[151] >= -2
 R5924: sout[135] - sout[143] + sout[159] + sout[167] - sout[183]
   - sout[191] + xorout[151] >= -2
 R5925: - sout[135] + sout[143] + sout[159] + sout[167] - sout[183]
   - sout[191] + xorout[151] >= -2
 R5926: sout[135] + sout[143] - sout[159] - sout[167] + sout[183]
   - sout[191] + xorout[151] >= -2
 R5927: sout[135] - sout[143] + sout[159] - sout[167] + sout[183]
   - sout[191] + xorout[151] >= -2
 R5928: - sout[135] + sout[143] + sout[159] - sout[167] + sout[183]
   - sout[191] + xorout[151] >= -2
 R5929: sout[135] - sout[143] - sout[159] + sout[167] + sout[183]
   - sout[191] + xorout[151] >= -2
 R5930: - sout[135] + sout[143] - sout[159] + sout[167] + sout[183]
   - sout[191] + xorout[151] >= -2
 R5931: - sout[135] - sout[143] + sout[159] + sout[167] + sout[183]
   - sout[191] + xorout[151] >= -2
 R5932: sout[135] + sout[143] - sout[159] - sout[167] - sout[183]
   + sout[191] + xorout[151] >= -2
 R5933: sout[135] - sout[143] + sout[159] - sout[167] - sout[183]
   + sout[191] + xorout[151] >= -2
 R5934: - sout[135] + sout[143] + sout[159] - sout[167] - sout[183]
   + sout[191] + xorout[151] >= -2
 R5935: sout[135] - sout[143] - sout[159] + sout[167] - sout[183]
   + sout[191] + xorout[151] >= -2
 R5936: - sout[135] + sout[143] - sout[159] + sout[167] - sout[183]
   + sout[191] + xorout[151] >= -2
 R5937: - sout[135] - sout[143] + sout[159] + sout[167] - sout[183]
   + sout[191] + xorout[151] >= -2
 R5938: sout[135] - sout[143] - sout[159] - sout[167] + sout[183]
   + sout[191] + xorout[151] >= -2
 R5939: - sout[135] + sout[143] - sout[159] - sout[167] + sout[183]
   + sout[191] + xorout[151] >= -2
 R5940: - sout[135] - sout[143] + sout[159] - sout[167] + sout[183]
   + sout[191] + xorout[151] >= -2
 R5941: - sout[135] - sout[143] - sout[159] + sout[167] + sout[183]
   + sout[191] + xorout[151] >= -2
 R5942: sout[135] - sout[143] - sout[159] - sout[167] - sout[183]
   - sout[191] + xorout[151] >= -4
 R5943: - sout[135] + sout[143] - sout[159] - sout[167] - sout[183]
   - sout[191] + xorout[151] >= -4
 R5944: - sout[135] - sout[143] + sout[159] - sout[167] - sout[183]
   - sout[191] + xorout[151] >= -4
 R5945: - sout[135] - sout[143] - sout[159] + sout[167] - sout[183]
   - sout[191] + xorout[151] >= -4
 R5946: - sout[135] - sout[143] - sout[159] - sout[167] + sout[183]
   - sout[191] + xorout[151] >= -4
 R5947: - sout[135] - sout[143] - sout[159] - sout[167] - sout[183]
   + sout[191] + xorout[151] >= -4
 R5948: - sout[190] + xorout[150] >= 0
 R5949: - sout[182] + xorout[150] >= 0
 R5950: - sout[166] + xorout[150] >= 0
 R5951: - sout[158] + xorout[150] >= 0
 R5952: - sout[142] + xorout[150] >= 0
 R5953: - sout[134] + xorout[150] >= 0
 R5954: - sout[134] + sout[135] >= 0
 R5955: - sout[142] + sout[143] >= 0
 R5956: - sout[158] + sout[159] >= 0
 R5957: - sout[166] + sout[167] >= 0
 R5958: - sout[182] + sout[183] >= 0
 R5959: - sout[190] + sout[191] >= 0
 R5960: sout[134] + sout[142] + sout[158] + sout[166] + sout[182]
   + sout[190] - xorout[150] >= 0
 R5961: - xorout[150] + xorout[151] >= 0
 R5962: - sout[135] - sout[143] - sout[159] - sout[167] - sout[183]
   - sout[191] + xorout[150] - xorout[151] >= -6
 R5963: sout[135] + sout[143] - sout[159] - sout[167] - sout[183]
   - sout[191] + xorout[150] - xorout[151] >= -4
 R5964: sout[135] - sout[143] + sout[159] - sout[167] - sout[183]
   - sout[191] + xorout[150] - xorout[151] >= -4
 R5965: - sout[135] + sout[143] + sout[159] - sout[167] - sout[183]
   - sout[191] + xorout[150] - xorout[151] >= -4
 R5966: sout[135] - sout[143] - sout[159] + sout[167] - sout[183]
   - sout[191] + xorout[150] - xorout[151] >= -4
 R5967: - sout[135] + sout[143] - sout[159] + sout[167] - sout[183]
   - sout[191] + xorout[150] - xorout[151] >= -4
 R5968: - sout[135] - sout[143] + sout[159] + sout[167] - sout[183]
   - sout[191] + xorout[150] - xorout[151] >= -4
 R5969: sout[135] - sout[143] - sout[159] - sout[167] + sout[183]
   - sout[191] + xorout[150] - xorout[151] >= -4
 R5970: - sout[135] + sout[143] - sout[159] - sout[167] + sout[183]
   - sout[191] + xorout[150] - xorout[151] >= -4
 R5971: - sout[135] - sout[143] + sout[159] - sout[167] + sout[183]
   - sout[191] + xorout[150] - xorout[151] >= -4
 R5972: - sout[135] - sout[143] - sout[159] + sout[167] + sout[183]
   - sout[191] + xorout[150] - xorout[151] >= -4
 R5973: sout[135] - sout[143] - sout[159] - sout[167] - sout[183]
   + sout[191] + xorout[150] - xorout[151] >= -4
 R5974: - sout[135] + sout[143] - sout[159] - sout[167] - sout[183]
   + sout[191] + xorout[150] - xorout[151] >= -4
 R5975: - sout[135] - sout[143] + sout[159] - sout[167] - sout[183]
   + sout[191] + xorout[150] - xorout[151] >= -4
 R5976: - sout[135] - sout[143] - sout[159] + sout[167] - sout[183]
   + sout[191] + xorout[150] - xorout[151] >= -4
 R5977: - sout[135] - sout[143] - sout[159] - sout[167] + sout[183]
   + sout[191] + xorout[150] - xorout[151] >= -4
 R5978: sout[135] + sout[143] + sout[159] + sout[167] - sout[183]
   - sout[191] + xorout[150] - xorout[151] >= -2
 R5979: sout[135] + sout[143] + sout[159] - sout[167] + sout[183]
   - sout[191] + xorout[150] - xorout[151] >= -2
 R5980: sout[135] + sout[143] - sout[159] + sout[167] + sout[183]
   - sout[191] + xorout[150] - xorout[151] >= -2
 R5981: sout[135] - sout[143] + sout[159] + sout[167] + sout[183]
   - sout[191] + xorout[150] - xorout[151] >= -2
 R5982: - sout[135] + sout[143] + sout[159] + sout[167] + sout[183]
   - sout[191] + xorout[150] - xorout[151] >= -2
 R5983: sout[135] + sout[143] + sout[159] - sout[167] - sout[183]
   + sout[191] + xorout[150] - xorout[151] >= -2
 R5984: sout[135] + sout[143] - sout[159] + sout[167] - sout[183]
   + sout[191] + xorout[150] - xorout[151] >= -2
 R5985: sout[135] - sout[143] + sout[159] + sout[167] - sout[183]
   + sout[191] + xorout[150] - xorout[151] >= -2
 R5986: - sout[135] + sout[143] + sout[159] + sout[167] - sout[183]
   + sout[191] + xorout[150] - xorout[151] >= -2
 R5987: sout[135] + sout[143] - sout[159] - sout[167] + sout[183]
   + sout[191] + xorout[150] - xorout[151] >= -2
 R5988: sout[135] - sout[143] + sout[159] - sout[167] + sout[183]
   + sout[191] + xorout[150] - xorout[151] >= -2
 R5989: - sout[135] + sout[143] + sout[159] - sout[167] + sout[183]
   + sout[191] + xorout[150] - xorout[151] >= -2
 R5990: sout[135] - sout[143] - sout[159] + sout[167] + sout[183]
   + sout[191] + xorout[150] - xorout[151] >= -2
 R5991: - sout[135] + sout[143] - sout[159] + sout[167] + sout[183]
   + sout[191] + xorout[150] - xorout[151] >= -2
 R5992: - sout[135] - sout[143] + sout[159] + sout[167] + sout[183]
   + sout[191] + xorout[150] - xorout[151] >= -2
 R5993: sout[129] + sout[137] + sout[145] + sout[161] + sout[169]
   - sout[185] + xorout[153] >= 0
 R5994: sout[129] + sout[137] + sout[145] + sout[161] - sout[169]
   + sout[185] + xorout[153] >= 0
 R5995: sout[129] + sout[137] + sout[145] - sout[161] + sout[169]
   + sout[185] + xorout[153] >= 0
 R5996: sout[129] + sout[137] - sout[145] + sout[161] + sout[169]
   + sout[185] + xorout[153] >= 0
 R5997: sout[129] - sout[137] + sout[145] + sout[161] + sout[169]
   + sout[185] + xorout[153] >= 0
 R5998: - sout[129] + sout[137] + sout[145] + sout[161] + sout[169]
   + sout[185] + xorout[153] >= 0
 R5999: sout[129] + sout[137] + sout[145] + sout[161] + sout[169]
   + sout[185] - xorout[153] >= 0
 R6000: sout[129] + sout[137] + sout[145] - sout[161] - sout[169]
   - sout[185] + xorout[153] >= -2
 R6001: sout[129] + sout[137] - sout[145] + sout[161] - sout[169]
   - sout[185] + xorout[153] >= -2
 R6002: sout[129] - sout[137] + sout[145] + sout[161] - sout[169]
   - sout[185] + xorout[153] >= -2
 R6003: - sout[129] + sout[137] + sout[145] + sout[161] - sout[169]
   - sout[185] + xorout[153] >= -2
 R6004: sout[129] + sout[137] - sout[145] - sout[161] + sout[169]
   - sout[185] + xorout[153] >= -2
 R6005: sout[129] - sout[137] + sout[145] - sout[161] + sout[169]
   - sout[185] + xorout[153] >= -2
 R6006: - sout[129] + sout[137] + sout[145] - sout[161] + sout[169]
   - sout[185] + xorout[153] >= -2
 R6007: sout[129] - sout[137] - sout[145] + sout[161] + sout[169]
   - sout[185] + xorout[153] >= -2
 R6008: - sout[129] + sout[137] - sout[145] + sout[161] + sout[169]
   - sout[185] + xorout[153] >= -2
 R6009: - sout[129] - sout[137] + sout[145] + sout[161] + sout[169]
   - sout[185] + xorout[153] >= -2
 R6010: sout[129] + sout[137] - sout[145] - sout[161] - sout[169]
   + sout[185] + xorout[153] >= -2
 R6011: sout[129] - sout[137] + sout[145] - sout[161] - sout[169]
   + sout[185] + xorout[153] >= -2
 R6012: - sout[129] + sout[137] + sout[145] - sout[161] - sout[169]
   + sout[185] + xorout[153] >= -2
 R6013: sout[129] - sout[137] - sout[145] + sout[161] - sout[169]
   + sout[185] + xorout[153] >= -2
 R6014: - sout[129] + sout[137] - sout[145] + sout[161] - sout[169]
   + sout[185] + xorout[153] >= -2
 R6015: - sout[129] - sout[137] + sout[145] + sout[161] - sout[169]
   + sout[185] + xorout[153] >= -2
 R6016: sout[129] - sout[137] - sout[145] - sout[161] + sout[169]
   + sout[185] + xorout[153] >= -2
 R6017: - sout[129] + sout[137] - sout[145] - sout[161] + sout[169]
   + sout[185] + xorout[153] >= -2
 R6018: - sout[129] - sout[137] + sout[145] - sout[161] + sout[169]
   + sout[185] + xorout[153] >= -2
 R6019: - sout[129] - sout[137] - sout[145] + sout[161] + sout[169]
   + sout[185] + xorout[153] >= -2
 R6020: sout[129] - sout[137] - sout[145] - sout[161] - sout[169]
   - sout[185] + xorout[153] >= -4
 R6021: - sout[129] + sout[137] - sout[145] - sout[161] - sout[169]
   - sout[185] + xorout[153] >= -4
 R6022: - sout[129] - sout[137] + sout[145] - sout[161] - sout[169]
   - sout[185] + xorout[153] >= -4
 R6023: - sout[129] - sout[137] - sout[145] + sout[161] - sout[169]
   - sout[185] + xorout[153] >= -4
 R6024: - sout[129] - sout[137] - sout[145] - sout[161] + sout[169]
   - sout[185] + xorout[153] >= -4
 R6025: - sout[129] - sout[137] - sout[145] - sout[161] - sout[169]
   + sout[185] + xorout[153] >= -4
 R6026: - sout[184] + xorout[152] >= 0
 R6027: - sout[168] + xorout[152] >= 0
 R6028: - sout[160] + xorout[152] >= 0
 R6029: - sout[144] + xorout[152] >= 0
 R6030: - sout[136] + xorout[152] >= 0
 R6031: - sout[128] + xorout[152] >= 0
 R6032: - sout[128] + sout[129] >= 0
 R6033: - sout[136] + sout[137] >= 0
 R6034: - sout[144] + sout[145] >= 0
 R6035: - sout[160] + sout[161] >= 0
 R6036: - sout[168] + sout[169] >= 0
 R6037: - sout[184] + sout[185] >= 0
 R6038: sout[128] + sout[136] + sout[144] + sout[160] + sout[168]
   + sout[184] - xorout[152] >= 0
 R6039: - xorout[152] + xorout[153] >= 0
 R6040: - sout[129] - sout[137] - sout[145] - sout[161] - sout[169]
   - sout[185] + xorout[152] - xorout[153] >= -6
 R6041: sout[129] + sout[137] - sout[145] - sout[161] - sout[169]
   - sout[185] + xorout[152] - xorout[153] >= -4
 R6042: sout[129] - sout[137] + sout[145] - sout[161] - sout[169]
   - sout[185] + xorout[152] - xorout[153] >= -4
 R6043: - sout[129] + sout[137] + sout[145] - sout[161] - sout[169]
   - sout[185] + xorout[152] - xorout[153] >= -4
 R6044: sout[129] - sout[137] - sout[145] + sout[161] - sout[169]
   - sout[185] + xorout[152] - xorout[153] >= -4
 R6045: - sout[129] + sout[137] - sout[145] + sout[161] - sout[169]
   - sout[185] + xorout[152] - xorout[153] >= -4
 R6046: - sout[129] - sout[137] + sout[145] + sout[161] - sout[169]
   - sout[185] + xorout[152] - xorout[153] >= -4
 R6047: sout[129] - sout[137] - sout[145] - sout[161] + sout[169]
   - sout[185] + xorout[152] - xorout[153] >= -4
 R6048: - sout[129] + sout[137] - sout[145] - sout[161] + sout[169]
   - sout[185] + xorout[152] - xorout[153] >= -4
 R6049: - sout[129] - sout[137] + sout[145] - sout[161] + sout[169]
   - sout[185] + xorout[152] - xorout[153] >= -4
 R6050: - sout[129] - sout[137] - sout[145] + sout[161] + sout[169]
   - sout[185] + xorout[152] - xorout[153] >= -4
 R6051: sout[129] - sout[137] - sout[145] - sout[161] - sout[169]
   + sout[185] + xorout[152] - xorout[153] >= -4
 R6052: - sout[129] + sout[137] - sout[145] - sout[161] - sout[169]
   + sout[185] + xorout[152] - xorout[153] >= -4
 R6053: - sout[129] - sout[137] + sout[145] - sout[161] - sout[169]
   + sout[185] + xorout[152] - xorout[153] >= -4
 R6054: - sout[129] - sout[137] - sout[145] + sout[161] - sout[169]
   + sout[185] + xorout[152] - xorout[153] >= -4
 R6055: - sout[129] - sout[137] - sout[145] - sout[161] + sout[169]
   + sout[185] + xorout[152] - xorout[153] >= -4
 R6056: sout[129] + sout[137] + sout[145] + sout[161] - sout[169]
   - sout[185] + xorout[152] - xorout[153] >= -2
 R6057: sout[129] + sout[137] + sout[145] - sout[161] + sout[169]
   - sout[185] + xorout[152] - xorout[153] >= -2
 R6058: sout[129] + sout[137] - sout[145] + sout[161] + sout[169]
   - sout[185] + xorout[152] - xorout[153] >= -2
 R6059: sout[129] - sout[137] + sout[145] + sout[161] + sout[169]
   - sout[185] + xorout[152] - xorout[153] >= -2
 R6060: - sout[129] + sout[137] + sout[145] + sout[161] + sout[169]
   - sout[185] + xorout[152] - xorout[153] >= -2
 R6061: sout[129] + sout[137] + sout[145] - sout[161] - sout[169]
   + sout[185] + xorout[152] - xorout[153] >= -2
 R6062: sout[129] + sout[137] - sout[145] + sout[161] - sout[169]
   + sout[185] + xorout[152] - xorout[153] >= -2
 R6063: sout[129] - sout[137] + sout[145] + sout[161] - sout[169]
   + sout[185] + xorout[152] - xorout[153] >= -2
 R6064: - sout[129] + sout[137] + sout[145] + sout[161] - sout[169]
   + sout[185] + xorout[152] - xorout[153] >= -2
 R6065: sout[129] + sout[137] - sout[145] - sout[161] + sout[169]
   + sout[185] + xorout[152] - xorout[153] >= -2
 R6066: sout[129] - sout[137] + sout[145] - sout[161] + sout[169]
   + sout[185] + xorout[152] - xorout[153] >= -2
 R6067: - sout[129] + sout[137] + sout[145] - sout[161] + sout[169]
   + sout[185] + xorout[152] - xorout[153] >= -2
 R6068: sout[129] - sout[137] - sout[145] + sout[161] + sout[169]
   + sout[185] + xorout[152] - xorout[153] >= -2
 R6069: - sout[129] + sout[137] - sout[145] + sout[161] + sout[169]
   + sout[185] + xorout[152] - xorout[153] >= -2
 R6070: - sout[129] - sout[137] + sout[145] + sout[161] + sout[169]
   + sout[185] + xorout[152] - xorout[153] >= -2
 R6071: sout[131] + sout[139] + sout[147] + sout[163] + sout[171]
   - sout[187] + xorout[155] >= 0
 R6072: sout[131] + sout[139] + sout[147] + sout[163] - sout[171]
   + sout[187] + xorout[155] >= 0
 R6073: sout[131] + sout[139] + sout[147] - sout[163] + sout[171]
   + sout[187] + xorout[155] >= 0
 R6074: sout[131] + sout[139] - sout[147] + sout[163] + sout[171]
   + sout[187] + xorout[155] >= 0
 R6075: sout[131] - sout[139] + sout[147] + sout[163] + sout[171]
   + sout[187] + xorout[155] >= 0
 R6076: - sout[131] + sout[139] + sout[147] + sout[163] + sout[171]
   + sout[187] + xorout[155] >= 0
 R6077: sout[131] + sout[139] + sout[147] + sout[163] + sout[171]
   + sout[187] - xorout[155] >= 0
 R6078: sout[131] + sout[139] + sout[147] - sout[163] - sout[171]
   - sout[187] + xorout[155] >= -2
 R6079: sout[131] + sout[139] - sout[147] + sout[163] - sout[171]
   - sout[187] + xorout[155] >= -2
 R6080: sout[131] - sout[139] + sout[147] + sout[163] - sout[171]
   - sout[187] + xorout[155] >= -2
 R6081: - sout[131] + sout[139] + sout[147] + sout[163] - sout[171]
   - sout[187] + xorout[155] >= -2
 R6082: sout[131] + sout[139] - sout[147] - sout[163] + sout[171]
   - sout[187] + xorout[155] >= -2
 R6083: sout[131] - sout[139] + sout[147] - sout[163] + sout[171]
   - sout[187] + xorout[155] >= -2
 R6084: - sout[131] + sout[139] + sout[147] - sout[163] + sout[171]
   - sout[187] + xorout[155] >= -2
 R6085: sout[131] - sout[139] - sout[147] + sout[163] + sout[171]
   - sout[187] + xorout[155] >= -2
 R6086: - sout[131] + sout[139] - sout[147] + sout[163] + sout[171]
   - sout[187] + xorout[155] >= -2
 R6087: - sout[131] - sout[139] + sout[147] + sout[163] + sout[171]
   - sout[187] + xorout[155] >= -2
 R6088: sout[131] + sout[139] - sout[147] - sout[163] - sout[171]
   + sout[187] + xorout[155] >= -2
 R6089: sout[131] - sout[139] + sout[147] - sout[163] - sout[171]
   + sout[187] + xorout[155] >= -2
 R6090: - sout[131] + sout[139] + sout[147] - sout[163] - sout[171]
   + sout[187] + xorout[155] >= -2
 R6091: sout[131] - sout[139] - sout[147] + sout[163] - sout[171]
   + sout[187] + xorout[155] >= -2
 R6092: - sout[131] + sout[139] - sout[147] + sout[163] - sout[171]
   + sout[187] + xorout[155] >= -2
 R6093: - sout[131] - sout[139] + sout[147] + sout[163] - sout[171]
   + sout[187] + xorout[155] >= -2
 R6094: sout[131] - sout[139] - sout[147] - sout[163] + sout[171]
   + sout[187] + xorout[155] >= -2
 R6095: - sout[131] + sout[139] - sout[147] - sout[163] + sout[171]
   + sout[187] + xorout[155] >= -2
 R6096: - sout[131] - sout[139] + sout[147] - sout[163] + sout[171]
   + sout[187] + xorout[155] >= -2
 R6097: - sout[131] - sout[139] - sout[147] + sout[163] + sout[171]
   + sout[187] + xorout[155] >= -2
 R6098: sout[131] - sout[139] - sout[147] - sout[163] - sout[171]
   - sout[187] + xorout[155] >= -4
 R6099: - sout[131] + sout[139] - sout[147] - sout[163] - sout[171]
   - sout[187] + xorout[155] >= -4
 R6100: - sout[131] - sout[139] + sout[147] - sout[163] - sout[171]
   - sout[187] + xorout[155] >= -4
 R6101: - sout[131] - sout[139] - sout[147] + sout[163] - sout[171]
   - sout[187] + xorout[155] >= -4
 R6102: - sout[131] - sout[139] - sout[147] - sout[163] + sout[171]
   - sout[187] + xorout[155] >= -4
 R6103: - sout[131] - sout[139] - sout[147] - sout[163] - sout[171]
   + sout[187] + xorout[155] >= -4
 R6104: - sout[186] + xorout[154] >= 0
 R6105: - sout[170] + xorout[154] >= 0
 R6106: - sout[162] + xorout[154] >= 0
 R6107: - sout[146] + xorout[154] >= 0
 R6108: - sout[138] + xorout[154] >= 0
 R6109: - sout[130] + xorout[154] >= 0
 R6110: - sout[130] + sout[131] >= 0
 R6111: - sout[138] + sout[139] >= 0
 R6112: - sout[146] + sout[147] >= 0
 R6113: - sout[162] + sout[163] >= 0
 R6114: - sout[170] + sout[171] >= 0
 R6115: - sout[186] + sout[187] >= 0
 R6116: sout[130] + sout[138] + sout[146] + sout[162] + sout[170]
   + sout[186] - xorout[154] >= 0
 R6117: - xorout[154] + xorout[155] >= 0
 R6118: - sout[131] - sout[139] - sout[147] - sout[163] - sout[171]
   - sout[187] + xorout[154] - xorout[155] >= -6
 R6119: sout[131] + sout[139] - sout[147] - sout[163] - sout[171]
   - sout[187] + xorout[154] - xorout[155] >= -4
 R6120: sout[131] - sout[139] + sout[147] - sout[163] - sout[171]
   - sout[187] + xorout[154] - xorout[155] >= -4
 R6121: - sout[131] + sout[139] + sout[147] - sout[163] - sout[171]
   - sout[187] + xorout[154] - xorout[155] >= -4
 R6122: sout[131] - sout[139] - sout[147] + sout[163] - sout[171]
   - sout[187] + xorout[154] - xorout[155] >= -4
 R6123: - sout[131] + sout[139] - sout[147] + sout[163] - sout[171]
   - sout[187] + xorout[154] - xorout[155] >= -4
 R6124: - sout[131] - sout[139] + sout[147] + sout[163] - sout[171]
   - sout[187] + xorout[154] - xorout[155] >= -4
 R6125: sout[131] - sout[139] - sout[147] - sout[163] + sout[171]
   - sout[187] + xorout[154] - xorout[155] >= -4
 R6126: - sout[131] + sout[139] - sout[147] - sout[163] + sout[171]
   - sout[187] + xorout[154] - xorout[155] >= -4
 R6127: - sout[131] - sout[139] + sout[147] - sout[163] + sout[171]
   - sout[187] + xorout[154] - xorout[155] >= -4
 R6128: - sout[131] - sout[139] - sout[147] + sout[163] + sout[171]
   - sout[187] + xorout[154] - xorout[155] >= -4
 R6129: sout[131] - sout[139] - sout[147] - sout[163] - sout[171]
   + sout[187] + xorout[154] - xorout[155] >= -4
 R6130: - sout[131] + sout[139] - sout[147] - sout[163] - sout[171]
   + sout[187] + xorout[154] - xorout[155] >= -4
 R6131: - sout[131] - sout[139] + sout[147] - sout[163] - sout[171]
   + sout[187] + xorout[154] - xorout[155] >= -4
 R6132: - sout[131] - sout[139] - sout[147] + sout[163] - sout[171]
   + sout[187] + xorout[154] - xorout[155] >= -4
 R6133: - sout[131] - sout[139] - sout[147] - sout[163] + sout[171]
   + sout[187] + xorout[154] - xorout[155] >= -4
 R6134: sout[131] + sout[139] + sout[147] + sout[163] - sout[171]
   - sout[187] + xorout[154] - xorout[155] >= -2
 R6135: sout[131] + sout[139] + sout[147] - sout[163] + sout[171]
   - sout[187] + xorout[154] - xorout[155] >= -2
 R6136: sout[131] + sout[139] - sout[147] + sout[163] + sout[171]
   - sout[187] + xorout[154] - xorout[155] >= -2
 R6137: sout[131] - sout[139] + sout[147] + sout[163] + sout[171]
   - sout[187] + xorout[154] - xorout[155] >= -2
 R6138: - sout[131] + sout[139] + sout[147] + sout[163] + sout[171]
   - sout[187] + xorout[154] - xorout[155] >= -2
 R6139: sout[131] + sout[139] + sout[147] - sout[163] - sout[171]
   + sout[187] + xorout[154] - xorout[155] >= -2
 R6140: sout[131] + sout[139] - sout[147] + sout[163] - sout[171]
   + sout[187] + xorout[154] - xorout[155] >= -2
 R6141: sout[131] - sout[139] + sout[147] + sout[163] - sout[171]
   + sout[187] + xorout[154] - xorout[155] >= -2
 R6142: - sout[131] + sout[139] + sout[147] + sout[163] - sout[171]
   + sout[187] + xorout[154] - xorout[155] >= -2
 R6143: sout[131] + sout[139] - sout[147] - sout[163] + sout[171]
   + sout[187] + xorout[154] - xorout[155] >= -2
 R6144: sout[131] - sout[139] + sout[147] - sout[163] + sout[171]
   + sout[187] + xorout[154] - xorout[155] >= -2
 R6145: - sout[131] + sout[139] + sout[147] - sout[163] + sout[171]
   + sout[187] + xorout[154] - xorout[155] >= -2
 R6146: sout[131] - sout[139] - sout[147] + sout[163] + sout[171]
   + sout[187] + xorout[154] - xorout[155] >= -2
 R6147: - sout[131] + sout[139] - sout[147] + sout[163] + sout[171]
   + sout[187] + xorout[154] - xorout[155] >= -2
 R6148: - sout[131] - sout[139] + sout[147] + sout[163] + sout[171]
   + sout[187] + xorout[154] - xorout[155] >= -2
 R6149: sout[133] + sout[141] + sout[149] + sout[165] + sout[173]
   - sout[189] + xorout[157] >= 0
 R6150: sout[133] + sout[141] + sout[149] + sout[165] - sout[173]
   + sout[189] + xorout[157] >= 0
 R6151: sout[133] + sout[141] + sout[149] - sout[165] + sout[173]
   + sout[189] + xorout[157] >= 0
 R6152: sout[133] + sout[141] - sout[149] + sout[165] + sout[173]
   + sout[189] + xorout[157] >= 0
 R6153: sout[133] - sout[141] + sout[149] + sout[165] + sout[173]
   + sout[189] + xorout[157] >= 0
 R6154: - sout[133] + sout[141] + sout[149] + sout[165] + sout[173]
   + sout[189] + xorout[157] >= 0
 R6155: sout[133] + sout[141] + sout[149] + sout[165] + sout[173]
   + sout[189] - xorout[157] >= 0
 R6156: sout[133] + sout[141] + sout[149] - sout[165] - sout[173]
   - sout[189] + xorout[157] >= -2
 R6157: sout[133] + sout[141] - sout[149] + sout[165] - sout[173]
   - sout[189] + xorout[157] >= -2
 R6158: sout[133] - sout[141] + sout[149] + sout[165] - sout[173]
   - sout[189] + xorout[157] >= -2
 R6159: - sout[133] + sout[141] + sout[149] + sout[165] - sout[173]
   - sout[189] + xorout[157] >= -2
 R6160: sout[133] + sout[141] - sout[149] - sout[165] + sout[173]
   - sout[189] + xorout[157] >= -2
 R6161: sout[133] - sout[141] + sout[149] - sout[165] + sout[173]
   - sout[189] + xorout[157] >= -2
 R6162: - sout[133] + sout[141] + sout[149] - sout[165] + sout[173]
   - sout[189] + xorout[157] >= -2
 R6163: sout[133] - sout[141] - sout[149] + sout[165] + sout[173]
   - sout[189] + xorout[157] >= -2
 R6164: - sout[133] + sout[141] - sout[149] + sout[165] + sout[173]
   - sout[189] + xorout[157] >= -2
 R6165: - sout[133] - sout[141] + sout[149] + sout[165] + sout[173]
   - sout[189] + xorout[157] >= -2
 R6166: sout[133] + sout[141] - sout[149] - sout[165] - sout[173]
   + sout[189] + xorout[157] >= -2
 R6167: sout[133] - sout[141] + sout[149] - sout[165] - sout[173]
   + sout[189] + xorout[157] >= -2
 R6168: - sout[133] + sout[141] + sout[149] - sout[165] - sout[173]
   + sout[189] + xorout[157] >= -2
 R6169: sout[133] - sout[141] - sout[149] + sout[165] - sout[173]
   + sout[189] + xorout[157] >= -2
 R6170: - sout[133] + sout[141] - sout[149] + sout[165] - sout[173]
   + sout[189] + xorout[157] >= -2
 R6171: - sout[133] - sout[141] + sout[149] + sout[165] - sout[173]
   + sout[189] + xorout[157] >= -2
 R6172: sout[133] - sout[141] - sout[149] - sout[165] + sout[173]
   + sout[189] + xorout[157] >= -2
 R6173: - sout[133] + sout[141] - sout[149] - sout[165] + sout[173]
   + sout[189] + xorout[157] >= -2
 R6174: - sout[133] - sout[141] + sout[149] - sout[165] + sout[173]
   + sout[189] + xorout[157] >= -2
 R6175: - sout[133] - sout[141] - sout[149] + sout[165] + sout[173]
   + sout[189] + xorout[157] >= -2
 R6176: sout[133] - sout[141] - sout[149] - sout[165] - sout[173]
   - sout[189] + xorout[157] >= -4
 R6177: - sout[133] + sout[141] - sout[149] - sout[165] - sout[173]
   - sout[189] + xorout[157] >= -4
 R6178: - sout[133] - sout[141] + sout[149] - sout[165] - sout[173]
   - sout[189] + xorout[157] >= -4
 R6179: - sout[133] - sout[141] - sout[149] + sout[165] - sout[173]
   - sout[189] + xorout[157] >= -4
 R6180: - sout[133] - sout[141] - sout[149] - sout[165] + sout[173]
   - sout[189] + xorout[157] >= -4
 R6181: - sout[133] - sout[141] - sout[149] - sout[165] - sout[173]
   + sout[189] + xorout[157] >= -4
 R6182: - sout[188] + xorout[156] >= 0
 R6183: - sout[172] + xorout[156] >= 0
 R6184: - sout[164] + xorout[156] >= 0
 R6185: - sout[148] + xorout[156] >= 0
 R6186: - sout[140] + xorout[156] >= 0
 R6187: - sout[132] + xorout[156] >= 0
 R6188: - sout[132] + sout[133] >= 0
 R6189: - sout[140] + sout[141] >= 0
 R6190: - sout[148] + sout[149] >= 0
 R6191: - sout[164] + sout[165] >= 0
 R6192: - sout[172] + sout[173] >= 0
 R6193: - sout[188] + sout[189] >= 0
 R6194: sout[132] + sout[140] + sout[148] + sout[164] + sout[172]
   + sout[188] - xorout[156] >= 0
 R6195: - xorout[156] + xorout[157] >= 0
 R6196: - sout[133] - sout[141] - sout[149] - sout[165] - sout[173]
   - sout[189] + xorout[156] - xorout[157] >= -6
 R6197: sout[133] + sout[141] - sout[149] - sout[165] - sout[173]
   - sout[189] + xorout[156] - xorout[157] >= -4
 R6198: sout[133] - sout[141] + sout[149] - sout[165] - sout[173]
   - sout[189] + xorout[156] - xorout[157] >= -4
 R6199: - sout[133] + sout[141] + sout[149] - sout[165] - sout[173]
   - sout[189] + xorout[156] - xorout[157] >= -4
 R6200: sout[133] - sout[141] - sout[149] + sout[165] - sout[173]
   - sout[189] + xorout[156] - xorout[157] >= -4
 R6201: - sout[133] + sout[141] - sout[149] + sout[165] - sout[173]
   - sout[189] + xorout[156] - xorout[157] >= -4
 R6202: - sout[133] - sout[141] + sout[149] + sout[165] - sout[173]
   - sout[189] + xorout[156] - xorout[157] >= -4
 R6203: sout[133] - sout[141] - sout[149] - sout[165] + sout[173]
   - sout[189] + xorout[156] - xorout[157] >= -4
 R6204: - sout[133] + sout[141] - sout[149] - sout[165] + sout[173]
   - sout[189] + xorout[156] - xorout[157] >= -4
 R6205: - sout[133] - sout[141] + sout[149] - sout[165] + sout[173]
   - sout[189] + xorout[156] - xorout[157] >= -4
 R6206: - sout[133] - sout[141] - sout[149] + sout[165] + sout[173]
   - sout[189] + xorout[156] - xorout[157] >= -4
 R6207: sout[133] - sout[141] - sout[149] - sout[165] - sout[173]
   + sout[189] + xorout[156] - xorout[157] >= -4
 R6208: - sout[133] + sout[141] - sout[149] - sout[165] - sout[173]
   + sout[189] + xorout[156] - xorout[157] >= -4
 R6209: - sout[133] - sout[141] + sout[149] - sout[165] - sout[173]
   + sout[189] + xorout[156] - xorout[157] >= -4
 R6210: - sout[133] - sout[141] - sout[149] + sout[165] - sout[173]
   + sout[189] + xorout[156] - xorout[157] >= -4
 R6211: - sout[133] - sout[141] - sout[149] - sout[165] + sout[173]
   + sout[189] + xorout[156] - xorout[157] >= -4
 R6212: sout[133] + sout[141] + sout[149] + sout[165] - sout[173]
   - sout[189] + xorout[156] - xorout[157] >= -2
 R6213: sout[133] + sout[141] + sout[149] - sout[165] + sout[173]
   - sout[189] + xorout[156] - xorout[157] >= -2
 R6214: sout[133] + sout[141] - sout[149] + sout[165] + sout[173]
   - sout[189] + xorout[156] - xorout[157] >= -2
 R6215: sout[133] - sout[141] + sout[149] + sout[165] + sout[173]
   - sout[189] + xorout[156] - xorout[157] >= -2
 R6216: - sout[133] + sout[141] + sout[149] + sout[165] + sout[173]
   - sout[189] + xorout[156] - xorout[157] >= -2
 R6217: sout[133] + sout[141] + sout[149] - sout[165] - sout[173]
   + sout[189] + xorout[156] - xorout[157] >= -2
 R6218: sout[133] + sout[141] - sout[149] + sout[165] - sout[173]
   + sout[189] + xorout[156] - xorout[157] >= -2
 R6219: sout[133] - sout[141] + sout[149] + sout[165] - sout[173]
   + sout[189] + xorout[156] - xorout[157] >= -2
 R6220: - sout[133] + sout[141] + sout[149] + sout[165] - sout[173]
   + sout[189] + xorout[156] - xorout[157] >= -2
 R6221: sout[133] + sout[141] - sout[149] - sout[165] + sout[173]
   + sout[189] + xorout[156] - xorout[157] >= -2
 R6222: sout[133] - sout[141] + sout[149] - sout[165] + sout[173]
   + sout[189] + xorout[156] - xorout[157] >= -2
 R6223: - sout[133] + sout[141] + sout[149] - sout[165] + sout[173]
   + sout[189] + xorout[156] - xorout[157] >= -2
 R6224: sout[133] - sout[141] - sout[149] + sout[165] + sout[173]
   + sout[189] + xorout[156] - xorout[157] >= -2
 R6225: - sout[133] + sout[141] - sout[149] + sout[165] + sout[173]
   + sout[189] + xorout[156] - xorout[157] >= -2
 R6226: - sout[133] - sout[141] + sout[149] + sout[165] + sout[173]
   + sout[189] + xorout[156] - xorout[157] >= -2
 R6227: sout[135] + sout[143] + sout[151] + sout[167] + sout[175]
   - sout[191] + xorout[159] >= 0
 R6228: sout[135] + sout[143] + sout[151] + sout[167] - sout[175]
   + sout[191] + xorout[159] >= 0
 R6229: sout[135] + sout[143] + sout[151] - sout[167] + sout[175]
   + sout[191] + xorout[159] >= 0
 R6230: sout[135] + sout[143] - sout[151] + sout[167] + sout[175]
   + sout[191] + xorout[159] >= 0
 R6231: sout[135] - sout[143] + sout[151] + sout[167] + sout[175]
   + sout[191] + xorout[159] >= 0
 R6232: - sout[135] + sout[143] + sout[151] + sout[167] + sout[175]
   + sout[191] + xorout[159] >= 0
 R6233: sout[135] + sout[143] + sout[151] + sout[167] + sout[175]
   + sout[191] - xorout[159] >= 0
 R6234: sout[135] + sout[143] + sout[151] - sout[167] - sout[175]
   - sout[191] + xorout[159] >= -2
 R6235: sout[135] + sout[143] - sout[151] + sout[167] - sout[175]
   - sout[191] + xorout[159] >= -2
 R6236: sout[135] - sout[143] + sout[151] + sout[167] - sout[175]
   - sout[191] + xorout[159] >= -2
 R6237: - sout[135] + sout[143] + sout[151] + sout[167] - sout[175]
   - sout[191] + xorout[159] >= -2
 R6238: sout[135] + sout[143] - sout[151] - sout[167] + sout[175]
   - sout[191] + xorout[159] >= -2
 R6239: sout[135] - sout[143] + sout[151] - sout[167] + sout[175]
   - sout[191] + xorout[159] >= -2
 R6240: - sout[135] + sout[143] + sout[151] - sout[167] + sout[175]
   - sout[191] + xorout[159] >= -2
 R6241: sout[135] - sout[143] - sout[151] + sout[167] + sout[175]
   - sout[191] + xorout[159] >= -2
 R6242: - sout[135] + sout[143] - sout[151] + sout[167] + sout[175]
   - sout[191] + xorout[159] >= -2
 R6243: - sout[135] - sout[143] + sout[151] + sout[167] + sout[175]
   - sout[191] + xorout[159] >= -2
 R6244: sout[135] + sout[143] - sout[151] - sout[167] - sout[175]
   + sout[191] + xorout[159] >= -2
 R6245: sout[135] - sout[143] + sout[151] - sout[167] - sout[175]
   + sout[191] + xorout[159] >= -2
 R6246: - sout[135] + sout[143] + sout[151] - sout[167] - sout[175]
   + sout[191] + xorout[159] >= -2
 R6247: sout[135] - sout[143] - sout[151] + sout[167] - sout[175]
   + sout[191] + xorout[159] >= -2
 R6248: - sout[135] + sout[143] - sout[151] + sout[167] - sout[175]
   + sout[191] + xorout[159] >= -2
 R6249: - sout[135] - sout[143] + sout[151] + sout[167] - sout[175]
   + sout[191] + xorout[159] >= -2
 R6250: sout[135] - sout[143] - sout[151] - sout[167] + sout[175]
   + sout[191] + xorout[159] >= -2
 R6251: - sout[135] + sout[143] - sout[151] - sout[167] + sout[175]
   + sout[191] + xorout[159] >= -2
 R6252: - sout[135] - sout[143] + sout[151] - sout[167] + sout[175]
   + sout[191] + xorout[159] >= -2
 R6253: - sout[135] - sout[143] - sout[151] + sout[167] + sout[175]
   + sout[191] + xorout[159] >= -2
 R6254: sout[135] - sout[143] - sout[151] - sout[167] - sout[175]
   - sout[191] + xorout[159] >= -4
 R6255: - sout[135] + sout[143] - sout[151] - sout[167] - sout[175]
   - sout[191] + xorout[159] >= -4
 R6256: - sout[135] - sout[143] + sout[151] - sout[167] - sout[175]
   - sout[191] + xorout[159] >= -4
 R6257: - sout[135] - sout[143] - sout[151] + sout[167] - sout[175]
   - sout[191] + xorout[159] >= -4
 R6258: - sout[135] - sout[143] - sout[151] - sout[167] + sout[175]
   - sout[191] + xorout[159] >= -4
 R6259: - sout[135] - sout[143] - sout[151] - sout[167] - sout[175]
   + sout[191] + xorout[159] >= -4
 R6260: - sout[190] + xorout[158] >= 0
 R6261: - sout[174] + xorout[158] >= 0
 R6262: - sout[166] + xorout[158] >= 0
 R6263: - sout[150] + xorout[158] >= 0
 R6264: - sout[142] + xorout[158] >= 0
 R6265: - sout[134] + xorout[158] >= 0
 R6266: - sout[134] + sout[135] >= 0
 R6267: - sout[142] + sout[143] >= 0
 R6268: - sout[150] + sout[151] >= 0
 R6269: - sout[166] + sout[167] >= 0
 R6270: - sout[174] + sout[175] >= 0
 R6271: - sout[190] + sout[191] >= 0
 R6272: sout[134] + sout[142] + sout[150] + sout[166] + sout[174]
   + sout[190] - xorout[158] >= 0
 R6273: - xorout[158] + xorout[159] >= 0
 R6274: - sout[135] - sout[143] - sout[151] - sout[167] - sout[175]
   - sout[191] + xorout[158] - xorout[159] >= -6
 R6275: sout[135] + sout[143] - sout[151] - sout[167] - sout[175]
   - sout[191] + xorout[158] - xorout[159] >= -4
 R6276: sout[135] - sout[143] + sout[151] - sout[167] - sout[175]
   - sout[191] + xorout[158] - xorout[159] >= -4
 R6277: - sout[135] + sout[143] + sout[151] - sout[167] - sout[175]
   - sout[191] + xorout[158] - xorout[159] >= -4
 R6278: sout[135] - sout[143] - sout[151] + sout[167] - sout[175]
   - sout[191] + xorout[158] - xorout[159] >= -4
 R6279: - sout[135] + sout[143] - sout[151] + sout[167] - sout[175]
   - sout[191] + xorout[158] - xorout[159] >= -4
 R6280: - sout[135] - sout[143] + sout[151] + sout[167] - sout[175]
   - sout[191] + xorout[158] - xorout[159] >= -4
 R6281: sout[135] - sout[143] - sout[151] - sout[167] + sout[175]
   - sout[191] + xorout[158] - xorout[159] >= -4
 R6282: - sout[135] + sout[143] - sout[151] - sout[167] + sout[175]
   - sout[191] + xorout[158] - xorout[159] >= -4
 R6283: - sout[135] - sout[143] + sout[151] - sout[167] + sout[175]
   - sout[191] + xorout[158] - xorout[159] >= -4
 R6284: - sout[135] - sout[143] - sout[151] + sout[167] + sout[175]
   - sout[191] + xorout[158] - xorout[159] >= -4
 R6285: sout[135] - sout[143] - sout[151] - sout[167] - sout[175]
   + sout[191] + xorout[158] - xorout[159] >= -4
 R6286: - sout[135] + sout[143] - sout[151] - sout[167] - sout[175]
   + sout[191] + xorout[158] - xorout[159] >= -4
 R6287: - sout[135] - sout[143] + sout[151] - sout[167] - sout[175]
   + sout[191] + xorout[158] - xorout[159] >= -4
 R6288: - sout[135] - sout[143] - sout[151] + sout[167] - sout[175]
   + sout[191] + xorout[158] - xorout[159] >= -4
 R6289: - sout[135] - sout[143] - sout[151] - sout[167] + sout[175]
   + sout[191] + xorout[158] - xorout[159] >= -4
 R6290: sout[135] + sout[143] + sout[151] + sout[167] - sout[175]
   - sout[191] + xorout[158] - xorout[159] >= -2
 R6291: sout[135] + sout[143] + sout[151] - sout[167] + sout[175]
   - sout[191] + xorout[158] - xorout[159] >= -2
 R6292: sout[135] + sout[143] - sout[151] + sout[167] + sout[175]
   - sout[191] + xorout[158] - xorout[159] >= -2
 R6293: sout[135] - sout[143] + sout[151] + sout[167] + sout[175]
   - sout[191] + xorout[158] - xorout[159] >= -2
 R6294: - sout[135] + sout[143] + sout[151] + sout[167] + sout[175]
   - sout[191] + xorout[158] - xorout[159] >= -2
 R6295: sout[135] + sout[143] + sout[151] - sout[167] - sout[175]
   + sout[191] + xorout[158] - xorout[159] >= -2
 R6296: sout[135] + sout[143] - sout[151] + sout[167] - sout[175]
   + sout[191] + xorout[158] - xorout[159] >= -2
 R6297: sout[135] - sout[143] + sout[151] + sout[167] - sout[175]
   + sout[191] + xorout[158] - xorout[159] >= -2
 R6298: - sout[135] + sout[143] + sout[151] + sout[167] - sout[175]
   + sout[191] + xorout[158] - xorout[159] >= -2
 R6299: sout[135] + sout[143] - sout[151] - sout[167] + sout[175]
   + sout[191] + xorout[158] - xorout[159] >= -2
 R6300: sout[135] - sout[143] + sout[151] - sout[167] + sout[175]
   + sout[191] + xorout[158] - xorout[159] >= -2
 R6301: - sout[135] + sout[143] + sout[151] - sout[167] + sout[175]
   + sout[191] + xorout[158] - xorout[159] >= -2
 R6302: sout[135] - sout[143] - sout[151] + sout[167] + sout[175]
   + sout[191] + xorout[158] - xorout[159] >= -2
 R6303: - sout[135] + sout[143] - sout[151] + sout[167] + sout[175]
   + sout[191] + xorout[158] - xorout[159] >= -2
 R6304: - sout[135] - sout[143] + sout[151] + sout[167] + sout[175]
   + sout[191] + xorout[158] - xorout[159] >= -2
 R6305: sout[129] + sout[137] + sout[153] + sout[161] - sout[169]
   + xorout[161] >= 0
 R6306: sout[129] + sout[137] + sout[153] - sout[161] + sout[169]
   + xorout[161] >= 0
 R6307: sout[129] + sout[137] - sout[153] + sout[161] + sout[169]
   + xorout[161] >= 0
 R6308: sout[129] - sout[137] + sout[153] + sout[161] + sout[169]
   + xorout[161] >= 0
 R6309: - sout[129] + sout[137] + sout[153] + sout[161] + sout[169]
   + xorout[161] >= 0
 R6310: sout[129] + sout[137] + sout[153] + sout[161] + sout[169]
   - xorout[161] >= 0
 R6311: sout[129] + sout[137] - sout[153] - sout[161] - sout[169]
   + xorout[161] >= -2
 R6312: sout[129] - sout[137] + sout[153] - sout[161] - sout[169]
   + xorout[161] >= -2
 R6313: - sout[129] + sout[137] + sout[153] - sout[161] - sout[169]
   + xorout[161] >= -2
 R6314: sout[129] - sout[137] - sout[153] + sout[161] - sout[169]
   + xorout[161] >= -2
 R6315: - sout[129] + sout[137] - sout[153] + sout[161] - sout[169]
   + xorout[161] >= -2
 R6316: - sout[129] - sout[137] + sout[153] + sout[161] - sout[169]
   + xorout[161] >= -2
 R6317: sout[129] - sout[137] - sout[153] - sout[161] + sout[169]
   + xorout[161] >= -2
 R6318: - sout[129] + sout[137] - sout[153] - sout[161] + sout[169]
   + xorout[161] >= -2
 R6319: - sout[129] - sout[137] + sout[153] - sout[161] + sout[169]
   + xorout[161] >= -2
 R6320: - sout[129] - sout[137] - sout[153] + sout[161] + sout[169]
   + xorout[161] >= -2
 R6321: - sout[129] - sout[137] - sout[153] - sout[161] - sout[169]
   + xorout[161] >= -4
 R6322: - sout[168] + xorout[160] >= 0
 R6323: - sout[160] + xorout[160] >= 0
 R6324: - sout[152] + xorout[160] >= 0
 R6325: - sout[136] + xorout[160] >= 0
 R6326: - sout[128] + xorout[160] >= 0
 R6327: - sout[128] + sout[129] >= 0
 R6328: - sout[136] + sout[137] >= 0
 R6329: - sout[152] + sout[153] >= 0
 R6330: - sout[160] + sout[161] >= 0
 R6331: - sout[168] + sout[169] >= 0
 R6332: sout[128] + sout[136] + sout[152] + sout[160] + sout[168]
   - xorout[160] >= 0
 R6333: - xorout[160] + xorout[161] >= 0
 R6334: sout[129] - sout[137] - sout[153] - sout[161] - sout[169]
   + xorout[160] - xorout[161] >= -4
 R6335: - sout[129] + sout[137] - sout[153] - sout[161] - sout[169]
   + xorout[160] - xorout[161] >= -4
 R6336: - sout[129] - sout[137] + sout[153] - sout[161] - sout[169]
   + xorout[160] - xorout[161] >= -4
 R6337: - sout[129] - sout[137] - sout[153] + sout[161] - sout[169]
   + xorout[160] - xorout[161] >= -4
 R6338: - sout[129] - sout[137] - sout[153] - sout[161] + sout[169]
   + xorout[160] - xorout[161] >= -4
 R6339: sout[129] + sout[137] + sout[153] - sout[161] - sout[169]
   + xorout[160] - xorout[161] >= -2
 R6340: sout[129] + sout[137] - sout[153] + sout[161] - sout[169]
   + xorout[160] - xorout[161] >= -2
 R6341: sout[129] - sout[137] + sout[153] + sout[161] - sout[169]
   + xorout[160] - xorout[161] >= -2
 R6342: - sout[129] + sout[137] + sout[153] + sout[161] - sout[169]
   + xorout[160] - xorout[161] >= -2
 R6343: sout[129] + sout[137] - sout[153] - sout[161] + sout[169]
   + xorout[160] - xorout[161] >= -2
 R6344: sout[129] - sout[137] + sout[153] - sout[161] + sout[169]
   + xorout[160] - xorout[161] >= -2
 R6345: - sout[129] + sout[137] + sout[153] - sout[161] + sout[169]
   + xorout[160] - xorout[161] >= -2
 R6346: sout[129] - sout[137] - sout[153] + sout[161] + sout[169]
   + xorout[160] - xorout[161] >= -2
 R6347: - sout[129] + sout[137] - sout[153] + sout[161] + sout[169]
   + xorout[160] - xorout[161] >= -2
 R6348: - sout[129] - sout[137] + sout[153] + sout[161] + sout[169]
   + xorout[160] - xorout[161] >= -2
 R6349: sout[131] + sout[139] + sout[155] + sout[163] - sout[171]
   + xorout[163] >= 0
 R6350: sout[131] + sout[139] + sout[155] - sout[163] + sout[171]
   + xorout[163] >= 0
 R6351: sout[131] + sout[139] - sout[155] + sout[163] + sout[171]
   + xorout[163] >= 0
 R6352: sout[131] - sout[139] + sout[155] + sout[163] + sout[171]
   + xorout[163] >= 0
 R6353: - sout[131] + sout[139] + sout[155] + sout[163] + sout[171]
   + xorout[163] >= 0
 R6354: sout[131] + sout[139] + sout[155] + sout[163] + sout[171]
   - xorout[163] >= 0
 R6355: sout[131] + sout[139] - sout[155] - sout[163] - sout[171]
   + xorout[163] >= -2
 R6356: sout[131] - sout[139] + sout[155] - sout[163] - sout[171]
   + xorout[163] >= -2
 R6357: - sout[131] + sout[139] + sout[155] - sout[163] - sout[171]
   + xorout[163] >= -2
 R6358: sout[131] - sout[139] - sout[155] + sout[163] - sout[171]
   + xorout[163] >= -2
 R6359: - sout[131] + sout[139] - sout[155] + sout[163] - sout[171]
   + xorout[163] >= -2
 R6360: - sout[131] - sout[139] + sout[155] + sout[163] - sout[171]
   + xorout[163] >= -2
 R6361: sout[131] - sout[139] - sout[155] - sout[163] + sout[171]
   + xorout[163] >= -2
 R6362: - sout[131] + sout[139] - sout[155] - sout[163] + sout[171]
   + xorout[163] >= -2
 R6363: - sout[131] - sout[139] + sout[155] - sout[163] + sout[171]
   + xorout[163] >= -2
 R6364: - sout[131] - sout[139] - sout[155] + sout[163] + sout[171]
   + xorout[163] >= -2
 R6365: - sout[131] - sout[139] - sout[155] - sout[163] - sout[171]
   + xorout[163] >= -4
 R6366: - sout[170] + xorout[162] >= 0
 R6367: - sout[162] + xorout[162] >= 0
 R6368: - sout[154] + xorout[162] >= 0
 R6369: - sout[138] + xorout[162] >= 0
 R6370: - sout[130] + xorout[162] >= 0
 R6371: - sout[130] + sout[131] >= 0
 R6372: - sout[138] + sout[139] >= 0
 R6373: - sout[154] + sout[155] >= 0
 R6374: - sout[162] + sout[163] >= 0
 R6375: - sout[170] + sout[171] >= 0
 R6376: sout[130] + sout[138] + sout[154] + sout[162] + sout[170]
   - xorout[162] >= 0
 R6377: - xorout[162] + xorout[163] >= 0
 R6378: sout[131] - sout[139] - sout[155] - sout[163] - sout[171]
   + xorout[162] - xorout[163] >= -4
 R6379: - sout[131] + sout[139] - sout[155] - sout[163] - sout[171]
   + xorout[162] - xorout[163] >= -4
 R6380: - sout[131] - sout[139] + sout[155] - sout[163] - sout[171]
   + xorout[162] - xorout[163] >= -4
 R6381: - sout[131] - sout[139] - sout[155] + sout[163] - sout[171]
   + xorout[162] - xorout[163] >= -4
 R6382: - sout[131] - sout[139] - sout[155] - sout[163] + sout[171]
   + xorout[162] - xorout[163] >= -4
 R6383: sout[131] + sout[139] + sout[155] - sout[163] - sout[171]
   + xorout[162] - xorout[163] >= -2
 R6384: sout[131] + sout[139] - sout[155] + sout[163] - sout[171]
   + xorout[162] - xorout[163] >= -2
 R6385: sout[131] - sout[139] + sout[155] + sout[163] - sout[171]
   + xorout[162] - xorout[163] >= -2
 R6386: - sout[131] + sout[139] + sout[155] + sout[163] - sout[171]
   + xorout[162] - xorout[163] >= -2
 R6387: sout[131] + sout[139] - sout[155] - sout[163] + sout[171]
   + xorout[162] - xorout[163] >= -2
 R6388: sout[131] - sout[139] + sout[155] - sout[163] + sout[171]
   + xorout[162] - xorout[163] >= -2
 R6389: - sout[131] + sout[139] + sout[155] - sout[163] + sout[171]
   + xorout[162] - xorout[163] >= -2
 R6390: sout[131] - sout[139] - sout[155] + sout[163] + sout[171]
   + xorout[162] - xorout[163] >= -2
 R6391: - sout[131] + sout[139] - sout[155] + sout[163] + sout[171]
   + xorout[162] - xorout[163] >= -2
 R6392: - sout[131] - sout[139] + sout[155] + sout[163] + sout[171]
   + xorout[162] - xorout[163] >= -2
 R6393: sout[133] + sout[141] + sout[157] + sout[165] - sout[173]
   + xorout[165] >= 0
 R6394: sout[133] + sout[141] + sout[157] - sout[165] + sout[173]
   + xorout[165] >= 0
 R6395: sout[133] + sout[141] - sout[157] + sout[165] + sout[173]
   + xorout[165] >= 0
 R6396: sout[133] - sout[141] + sout[157] + sout[165] + sout[173]
   + xorout[165] >= 0
 R6397: - sout[133] + sout[141] + sout[157] + sout[165] + sout[173]
   + xorout[165] >= 0
 R6398: sout[133] + sout[141] + sout[157] + sout[165] + sout[173]
   - xorout[165] >= 0
 R6399: sout[133] + sout[141] - sout[157] - sout[165] - sout[173]
   + xorout[165] >= -2
 R6400: sout[133] - sout[141] + sout[157] - sout[165] - sout[173]
   + xorout[165] >= -2
 R6401: - sout[133] + sout[141] + sout[157] - sout[165] - sout[173]
   + xorout[165] >= -2
 R6402: sout[133] - sout[141] - sout[157] + sout[165] - sout[173]
   + xorout[165] >= -2
 R6403: - sout[133] + sout[141] - sout[157] + sout[165] - sout[173]
   + xorout[165] >= -2
 R6404: - sout[133] - sout[141] + sout[157] + sout[165] - sout[173]
   + xorout[165] >= -2
 R6405: sout[133] - sout[141] - sout[157] - sout[165] + sout[173]
   + xorout[165] >= -2
 R6406: - sout[133] + sout[141] - sout[157] - sout[165] + sout[173]
   + xorout[165] >= -2
 R6407: - sout[133] - sout[141] + sout[157] - sout[165] + sout[173]
   + xorout[165] >= -2
 R6408: - sout[133] - sout[141] - sout[157] + sout[165] + sout[173]
   + xorout[165] >= -2
 R6409: - sout[133] - sout[141] - sout[157] - sout[165] - sout[173]
   + xorout[165] >= -4
 R6410: - sout[172] + xorout[164] >= 0
 R6411: - sout[164] + xorout[164] >= 0
 R6412: - sout[156] + xorout[164] >= 0
 R6413: - sout[140] + xorout[164] >= 0
 R6414: - sout[132] + xorout[164] >= 0
 R6415: - sout[132] + sout[133] >= 0
 R6416: - sout[140] + sout[141] >= 0
 R6417: - sout[156] + sout[157] >= 0
 R6418: - sout[164] + sout[165] >= 0
 R6419: - sout[172] + sout[173] >= 0
 R6420: sout[132] + sout[140] + sout[156] + sout[164] + sout[172]
   - xorout[164] >= 0
 R6421: - xorout[164] + xorout[165] >= 0
 R6422: sout[133] - sout[141] - sout[157] - sout[165] - sout[173]
   + xorout[164] - xorout[165] >= -4
 R6423: - sout[133] + sout[141] - sout[157] - sout[165] - sout[173]
   + xorout[164] - xorout[165] >= -4
 R6424: - sout[133] - sout[141] + sout[157] - sout[165] - sout[173]
   + xorout[164] - xorout[165] >= -4
 R6425: - sout[133] - sout[141] - sout[157] + sout[165] - sout[173]
   + xorout[164] - xorout[165] >= -4
 R6426: - sout[133] - sout[141] - sout[157] - sout[165] + sout[173]
   + xorout[164] - xorout[165] >= -4
 R6427: sout[133] + sout[141] + sout[157] - sout[165] - sout[173]
   + xorout[164] - xorout[165] >= -2
 R6428: sout[133] + sout[141] - sout[157] + sout[165] - sout[173]
   + xorout[164] - xorout[165] >= -2
 R6429: sout[133] - sout[141] + sout[157] + sout[165] - sout[173]
   + xorout[164] - xorout[165] >= -2
 R6430: - sout[133] + sout[141] + sout[157] + sout[165] - sout[173]
   + xorout[164] - xorout[165] >= -2
 R6431: sout[133] + sout[141] - sout[157] - sout[165] + sout[173]
   + xorout[164] - xorout[165] >= -2
 R6432: sout[133] - sout[141] + sout[157] - sout[165] + sout[173]
   + xorout[164] - xorout[165] >= -2
 R6433: - sout[133] + sout[141] + sout[157] - sout[165] + sout[173]
   + xorout[164] - xorout[165] >= -2
 R6434: sout[133] - sout[141] - sout[157] + sout[165] + sout[173]
   + xorout[164] - xorout[165] >= -2
 R6435: - sout[133] + sout[141] - sout[157] + sout[165] + sout[173]
   + xorout[164] - xorout[165] >= -2
 R6436: - sout[133] - sout[141] + sout[157] + sout[165] + sout[173]
   + xorout[164] - xorout[165] >= -2
 R6437: sout[135] + sout[143] + sout[159] + sout[167] - sout[175]
   + xorout[167] >= 0
 R6438: sout[135] + sout[143] + sout[159] - sout[167] + sout[175]
   + xorout[167] >= 0
 R6439: sout[135] + sout[143] - sout[159] + sout[167] + sout[175]
   + xorout[167] >= 0
 R6440: sout[135] - sout[143] + sout[159] + sout[167] + sout[175]
   + xorout[167] >= 0
 R6441: - sout[135] + sout[143] + sout[159] + sout[167] + sout[175]
   + xorout[167] >= 0
 R6442: sout[135] + sout[143] + sout[159] + sout[167] + sout[175]
   - xorout[167] >= 0
 R6443: sout[135] + sout[143] - sout[159] - sout[167] - sout[175]
   + xorout[167] >= -2
 R6444: sout[135] - sout[143] + sout[159] - sout[167] - sout[175]
   + xorout[167] >= -2
 R6445: - sout[135] + sout[143] + sout[159] - sout[167] - sout[175]
   + xorout[167] >= -2
 R6446: sout[135] - sout[143] - sout[159] + sout[167] - sout[175]
   + xorout[167] >= -2
 R6447: - sout[135] + sout[143] - sout[159] + sout[167] - sout[175]
   + xorout[167] >= -2
 R6448: - sout[135] - sout[143] + sout[159] + sout[167] - sout[175]
   + xorout[167] >= -2
 R6449: sout[135] - sout[143] - sout[159] - sout[167] + sout[175]
   + xorout[167] >= -2
 R6450: - sout[135] + sout[143] - sout[159] - sout[167] + sout[175]
   + xorout[167] >= -2
 R6451: - sout[135] - sout[143] + sout[159] - sout[167] + sout[175]
   + xorout[167] >= -2
 R6452: - sout[135] - sout[143] - sout[159] + sout[167] + sout[175]
   + xorout[167] >= -2
 R6453: - sout[135] - sout[143] - sout[159] - sout[167] - sout[175]
   + xorout[167] >= -4
 R6454: - sout[174] + xorout[166] >= 0
 R6455: - sout[166] + xorout[166] >= 0
 R6456: - sout[158] + xorout[166] >= 0
 R6457: - sout[142] + xorout[166] >= 0
 R6458: - sout[134] + xorout[166] >= 0
 R6459: - sout[134] + sout[135] >= 0
 R6460: - sout[142] + sout[143] >= 0
 R6461: - sout[158] + sout[159] >= 0
 R6462: - sout[166] + sout[167] >= 0
 R6463: - sout[174] + sout[175] >= 0
 R6464: sout[134] + sout[142] + sout[158] + sout[166] + sout[174]
   - xorout[166] >= 0
 R6465: - xorout[166] + xorout[167] >= 0
 R6466: sout[135] - sout[143] - sout[159] - sout[167] - sout[175]
   + xorout[166] - xorout[167] >= -4
 R6467: - sout[135] + sout[143] - sout[159] - sout[167] - sout[175]
   + xorout[166] - xorout[167] >= -4
 R6468: - sout[135] - sout[143] + sout[159] - sout[167] - sout[175]
   + xorout[166] - xorout[167] >= -4
 R6469: - sout[135] - sout[143] - sout[159] + sout[167] - sout[175]
   + xorout[166] - xorout[167] >= -4
 R6470: - sout[135] - sout[143] - sout[159] - sout[167] + sout[175]
   + xorout[166] - xorout[167] >= -4
 R6471: sout[135] + sout[143] + sout[159] - sout[167] - sout[175]
   + xorout[166] - xorout[167] >= -2
 R6472: sout[135] + sout[143] - sout[159] + sout[167] - sout[175]
   + xorout[166] - xorout[167] >= -2
 R6473: sout[135] - sout[143] + sout[159] + sout[167] - sout[175]
   + xorout[166] - xorout[167] >= -2
 R6474: - sout[135] + sout[143] + sout[159] + sout[167] - sout[175]
   + xorout[166] - xorout[167] >= -2
 R6475: sout[135] + sout[143] - sout[159] - sout[167] + sout[175]
   + xorout[166] - xorout[167] >= -2
 R6476: sout[135] - sout[143] + sout[159] - sout[167] + sout[175]
   + xorout[166] - xorout[167] >= -2
 R6477: - sout[135] + sout[143] + sout[159] - sout[167] + sout[175]
   + xorout[166] - xorout[167] >= -2
 R6478: sout[135] - sout[143] - sout[159] + sout[167] + sout[175]
   + xorout[166] - xorout[167] >= -2
 R6479: - sout[135] + sout[143] - sout[159] + sout[167] + sout[175]
   + xorout[166] - xorout[167] >= -2
 R6480: - sout[135] - sout[143] + sout[159] + sout[167] + sout[175]
   + xorout[166] - xorout[167] >= -2
 R6481: sout[129] + sout[137] + sout[145] + sout[169] - sout[177]
   + xorout[169] >= 0
 R6482: sout[129] + sout[137] + sout[145] - sout[169] + sout[177]
   + xorout[169] >= 0
 R6483: sout[129] + sout[137] - sout[145] + sout[169] + sout[177]
   + xorout[169] >= 0
 R6484: sout[129] - sout[137] + sout[145] + sout[169] + sout[177]
   + xorout[169] >= 0
 R6485: - sout[129] + sout[137] + sout[145] + sout[169] + sout[177]
   + xorout[169] >= 0
 R6486: sout[129] + sout[137] + sout[145] + sout[169] + sout[177]
   - xorout[169] >= 0
 R6487: sout[129] + sout[137] - sout[145] - sout[169] - sout[177]
   + xorout[169] >= -2
 R6488: sout[129] - sout[137] + sout[145] - sout[169] - sout[177]
   + xorout[169] >= -2
 R6489: - sout[129] + sout[137] + sout[145] - sout[169] - sout[177]
   + xorout[169] >= -2
 R6490: sout[129] - sout[137] - sout[145] + sout[169] - sout[177]
   + xorout[169] >= -2
 R6491: - sout[129] + sout[137] - sout[145] + sout[169] - sout[177]
   + xorout[169] >= -2
 R6492: - sout[129] - sout[137] + sout[145] + sout[169] - sout[177]
   + xorout[169] >= -2
 R6493: sout[129] - sout[137] - sout[145] - sout[169] + sout[177]
   + xorout[169] >= -2
 R6494: - sout[129] + sout[137] - sout[145] - sout[169] + sout[177]
   + xorout[169] >= -2
 R6495: - sout[129] - sout[137] + sout[145] - sout[169] + sout[177]
   + xorout[169] >= -2
 R6496: - sout[129] - sout[137] - sout[145] + sout[169] + sout[177]
   + xorout[169] >= -2
 R6497: - sout[129] - sout[137] - sout[145] - sout[169] - sout[177]
   + xorout[169] >= -4
 R6498: - sout[176] + xorout[168] >= 0
 R6499: - sout[168] + xorout[168] >= 0
 R6500: - sout[144] + xorout[168] >= 0
 R6501: - sout[136] + xorout[168] >= 0
 R6502: - sout[128] + xorout[168] >= 0
 R6503: - sout[128] + sout[129] >= 0
 R6504: - sout[136] + sout[137] >= 0
 R6505: - sout[144] + sout[145] >= 0
 R6506: - sout[168] + sout[169] >= 0
 R6507: - sout[176] + sout[177] >= 0
 R6508: sout[128] + sout[136] + sout[144] + sout[168] + sout[176]
   - xorout[168] >= 0
 R6509: - xorout[168] + xorout[169] >= 0
 R6510: sout[129] - sout[137] - sout[145] - sout[169] - sout[177]
   + xorout[168] - xorout[169] >= -4
 R6511: - sout[129] + sout[137] - sout[145] - sout[169] - sout[177]
   + xorout[168] - xorout[169] >= -4
 R6512: - sout[129] - sout[137] + sout[145] - sout[169] - sout[177]
   + xorout[168] - xorout[169] >= -4
 R6513: - sout[129] - sout[137] - sout[145] + sout[169] - sout[177]
   + xorout[168] - xorout[169] >= -4
 R6514: - sout[129] - sout[137] - sout[145] - sout[169] + sout[177]
   + xorout[168] - xorout[169] >= -4
 R6515: sout[129] + sout[137] + sout[145] - sout[169] - sout[177]
   + xorout[168] - xorout[169] >= -2
 R6516: sout[129] + sout[137] - sout[145] + sout[169] - sout[177]
   + xorout[168] - xorout[169] >= -2
 R6517: sout[129] - sout[137] + sout[145] + sout[169] - sout[177]
   + xorout[168] - xorout[169] >= -2
 R6518: - sout[129] + sout[137] + sout[145] + sout[169] - sout[177]
   + xorout[168] - xorout[169] >= -2
 R6519: sout[129] + sout[137] - sout[145] - sout[169] + sout[177]
   + xorout[168] - xorout[169] >= -2
 R6520: sout[129] - sout[137] + sout[145] - sout[169] + sout[177]
   + xorout[168] - xorout[169] >= -2
 R6521: - sout[129] + sout[137] + sout[145] - sout[169] + sout[177]
   + xorout[168] - xorout[169] >= -2
 R6522: sout[129] - sout[137] - sout[145] + sout[169] + sout[177]
   + xorout[168] - xorout[169] >= -2
 R6523: - sout[129] + sout[137] - sout[145] + sout[169] + sout[177]
   + xorout[168] - xorout[169] >= -2
 R6524: - sout[129] - sout[137] + sout[145] + sout[169] + sout[177]
   + xorout[168] - xorout[169] >= -2
 R6525: sout[131] + sout[139] + sout[147] + sout[171] - sout[179]
   + xorout[171] >= 0
 R6526: sout[131] + sout[139] + sout[147] - sout[171] + sout[179]
   + xorout[171] >= 0
 R6527: sout[131] + sout[139] - sout[147] + sout[171] + sout[179]
   + xorout[171] >= 0
 R6528: sout[131] - sout[139] + sout[147] + sout[171] + sout[179]
   + xorout[171] >= 0
 R6529: - sout[131] + sout[139] + sout[147] + sout[171] + sout[179]
   + xorout[171] >= 0
 R6530: sout[131] + sout[139] + sout[147] + sout[171] + sout[179]
   - xorout[171] >= 0
 R6531: sout[131] + sout[139] - sout[147] - sout[171] - sout[179]
   + xorout[171] >= -2
 R6532: sout[131] - sout[139] + sout[147] - sout[171] - sout[179]
   + xorout[171] >= -2
 R6533: - sout[131] + sout[139] + sout[147] - sout[171] - sout[179]
   + xorout[171] >= -2
 R6534: sout[131] - sout[139] - sout[147] + sout[171] - sout[179]
   + xorout[171] >= -2
 R6535: - sout[131] + sout[139] - sout[147] + sout[171] - sout[179]
   + xorout[171] >= -2
 R6536: - sout[131] - sout[139] + sout[147] + sout[171] - sout[179]
   + xorout[171] >= -2
 R6537: sout[131] - sout[139] - sout[147] - sout[171] + sout[179]
   + xorout[171] >= -2
 R6538: - sout[131] + sout[139] - sout[147] - sout[171] + sout[179]
   + xorout[171] >= -2
 R6539: - sout[131] - sout[139] + sout[147] - sout[171] + sout[179]
   + xorout[171] >= -2
 R6540: - sout[131] - sout[139] - sout[147] + sout[171] + sout[179]
   + xorout[171] >= -2
 R6541: - sout[131] - sout[139] - sout[147] - sout[171] - sout[179]
   + xorout[171] >= -4
 R6542: - sout[178] + xorout[170] >= 0
 R6543: - sout[170] + xorout[170] >= 0
 R6544: - sout[146] + xorout[170] >= 0
 R6545: - sout[138] + xorout[170] >= 0
 R6546: - sout[130] + xorout[170] >= 0
 R6547: - sout[130] + sout[131] >= 0
 R6548: - sout[138] + sout[139] >= 0
 R6549: - sout[146] + sout[147] >= 0
 R6550: - sout[170] + sout[171] >= 0
 R6551: - sout[178] + sout[179] >= 0
 R6552: sout[130] + sout[138] + sout[146] + sout[170] + sout[178]
   - xorout[170] >= 0
 R6553: - xorout[170] + xorout[171] >= 0
 R6554: sout[131] - sout[139] - sout[147] - sout[171] - sout[179]
   + xorout[170] - xorout[171] >= -4
 R6555: - sout[131] + sout[139] - sout[147] - sout[171] - sout[179]
   + xorout[170] - xorout[171] >= -4
 R6556: - sout[131] - sout[139] + sout[147] - sout[171] - sout[179]
   + xorout[170] - xorout[171] >= -4
 R6557: - sout[131] - sout[139] - sout[147] + sout[171] - sout[179]
   + xorout[170] - xorout[171] >= -4
 R6558: - sout[131] - sout[139] - sout[147] - sout[171] + sout[179]
   + xorout[170] - xorout[171] >= -4
 R6559: sout[131] + sout[139] + sout[147] - sout[171] - sout[179]
   + xorout[170] - xorout[171] >= -2
 R6560: sout[131] + sout[139] - sout[147] + sout[171] - sout[179]
   + xorout[170] - xorout[171] >= -2
 R6561: sout[131] - sout[139] + sout[147] + sout[171] - sout[179]
   + xorout[170] - xorout[171] >= -2
 R6562: - sout[131] + sout[139] + sout[147] + sout[171] - sout[179]
   + xorout[170] - xorout[171] >= -2
 R6563: sout[131] + sout[139] - sout[147] - sout[171] + sout[179]
   + xorout[170] - xorout[171] >= -2
 R6564: sout[131] - sout[139] + sout[147] - sout[171] + sout[179]
   + xorout[170] - xorout[171] >= -2
 R6565: - sout[131] + sout[139] + sout[147] - sout[171] + sout[179]
   + xorout[170] - xorout[171] >= -2
 R6566: sout[131] - sout[139] - sout[147] + sout[171] + sout[179]
   + xorout[170] - xorout[171] >= -2
 R6567: - sout[131] + sout[139] - sout[147] + sout[171] + sout[179]
   + xorout[170] - xorout[171] >= -2
 R6568: - sout[131] - sout[139] + sout[147] + sout[171] + sout[179]
   + xorout[170] - xorout[171] >= -2
 R6569: sout[133] + sout[141] + sout[149] + sout[173] - sout[181]
   + xorout[173] >= 0
 R6570: sout[133] + sout[141] + sout[149] - sout[173] + sout[181]
   + xorout[173] >= 0
 R6571: sout[133] + sout[141] - sout[149] + sout[173] + sout[181]
   + xorout[173] >= 0
 R6572: sout[133] - sout[141] + sout[149] + sout[173] + sout[181]
   + xorout[173] >= 0
 R6573: - sout[133] + sout[141] + sout[149] + sout[173] + sout[181]
   + xorout[173] >= 0
 R6574: sout[133] + sout[141] + sout[149] + sout[173] + sout[181]
   - xorout[173] >= 0
 R6575: sout[133] + sout[141] - sout[149] - sout[173] - sout[181]
   + xorout[173] >= -2
 R6576: sout[133] - sout[141] + sout[149] - sout[173] - sout[181]
   + xorout[173] >= -2
 R6577: - sout[133] + sout[141] + sout[149] - sout[173] - sout[181]
   + xorout[173] >= -2
 R6578: sout[133] - sout[141] - sout[149] + sout[173] - sout[181]
   + xorout[173] >= -2
 R6579: - sout[133] + sout[141] - sout[149] + sout[173] - sout[181]
   + xorout[173] >= -2
 R6580: - sout[133] - sout[141] + sout[149] + sout[173] - sout[181]
   + xorout[173] >= -2
 R6581: sout[133] - sout[141] - sout[149] - sout[173] + sout[181]
   + xorout[173] >= -2
 R6582: - sout[133] + sout[141] - sout[149] - sout[173] + sout[181]
   + xorout[173] >= -2
 R6583: - sout[133] - sout[141] + sout[149] - sout[173] + sout[181]
   + xorout[173] >= -2
 R6584: - sout[133] - sout[141] - sout[149] + sout[173] + sout[181]
   + xorout[173] >= -2
 R6585: - sout[133] - sout[141] - sout[149] - sout[173] - sout[181]
   + xorout[173] >= -4
 R6586: - sout[180] + xorout[172] >= 0
 R6587: - sout[172] + xorout[172] >= 0
 R6588: - sout[148] + xorout[172] >= 0
 R6589: - sout[140] + xorout[172] >= 0
 R6590: - sout[132] + xorout[172] >= 0
 R6591: - sout[132] + sout[133] >= 0
 R6592: - sout[140] + sout[141] >= 0
 R6593: - sout[148] + sout[149] >= 0
 R6594: - sout[172] + sout[173] >= 0
 R6595: - sout[180] + sout[181] >= 0
 R6596: sout[132] + sout[140] + sout[148] + sout[172] + sout[180]
   - xorout[172] >= 0
 R6597: - xorout[172] + xorout[173] >= 0
 R6598: sout[133] - sout[141] - sout[149] - sout[173] - sout[181]
   + xorout[172] - xorout[173] >= -4
 R6599: - sout[133] + sout[141] - sout[149] - sout[173] - sout[181]
   + xorout[172] - xorout[173] >= -4
 R6600: - sout[133] - sout[141] + sout[149] - sout[173] - sout[181]
   + xorout[172] - xorout[173] >= -4
 R6601: - sout[133] - sout[141] - sout[149] + sout[173] - sout[181]
   + xorout[172] - xorout[173] >= -4
 R6602: - sout[133] - sout[141] - sout[149] - sout[173] + sout[181]
   + xorout[172] - xorout[173] >= -4
 R6603: sout[133] + sout[141] + sout[149] - sout[173] - sout[181]
   + xorout[172] - xorout[173] >= -2
 R6604: sout[133] + sout[141] - sout[149] + sout[173] - sout[181]
   + xorout[172] - xorout[173] >= -2
 R6605: sout[133] - sout[141] + sout[149] + sout[173] - sout[181]
   + xorout[172] - xorout[173] >= -2
 R6606: - sout[133] + sout[141] + sout[149] + sout[173] - sout[181]
   + xorout[172] - xorout[173] >= -2
 R6607: sout[133] + sout[141] - sout[149] - sout[173] + sout[181]
   + xorout[172] - xorout[173] >= -2
 R6608: sout[133] - sout[141] + sout[149] - sout[173] + sout[181]
   + xorout[172] - xorout[173] >= -2
 R6609: - sout[133] + sout[141] + sout[149] - sout[173] + sout[181]
   + xorout[172] - xorout[173] >= -2
 R6610: sout[133] - sout[141] - sout[149] + sout[173] + sout[181]
   + xorout[172] - xorout[173] >= -2
 R6611: - sout[133] + sout[141] - sout[149] + sout[173] + sout[181]
   + xorout[172] - xorout[173] >= -2
 R6612: - sout[133] - sout[141] + sout[149] + sout[173] + sout[181]
   + xorout[172] - xorout[173] >= -2
 R6613: sout[135] + sout[143] + sout[151] + sout[175] - sout[183]
   + xorout[175] >= 0
 R6614: sout[135] + sout[143] + sout[151] - sout[175] + sout[183]
   + xorout[175] >= 0
 R6615: sout[135] + sout[143] - sout[151] + sout[175] + sout[183]
   + xorout[175] >= 0
 R6616: sout[135] - sout[143] + sout[151] + sout[175] + sout[183]
   + xorout[175] >= 0
 R6617: - sout[135] + sout[143] + sout[151] + sout[175] + sout[183]
   + xorout[175] >= 0
 R6618: sout[135] + sout[143] + sout[151] + sout[175] + sout[183]
   - xorout[175] >= 0
 R6619: sout[135] + sout[143] - sout[151] - sout[175] - sout[183]
   + xorout[175] >= -2
 R6620: sout[135] - sout[143] + sout[151] - sout[175] - sout[183]
   + xorout[175] >= -2
 R6621: - sout[135] + sout[143] + sout[151] - sout[175] - sout[183]
   + xorout[175] >= -2
 R6622: sout[135] - sout[143] - sout[151] + sout[175] - sout[183]
   + xorout[175] >= -2
 R6623: - sout[135] + sout[143] - sout[151] + sout[175] - sout[183]
   + xorout[175] >= -2
 R6624: - sout[135] - sout[143] + sout[151] + sout[175] - sout[183]
   + xorout[175] >= -2
 R6625: sout[135] - sout[143] - sout[151] - sout[175] + sout[183]
   + xorout[175] >= -2
 R6626: - sout[135] + sout[143] - sout[151] - sout[175] + sout[183]
   + xorout[175] >= -2
 R6627: - sout[135] - sout[143] + sout[151] - sout[175] + sout[183]
   + xorout[175] >= -2
 R6628: - sout[135] - sout[143] - sout[151] + sout[175] + sout[183]
   + xorout[175] >= -2
 R6629: - sout[135] - sout[143] - sout[151] - sout[175] - sout[183]
   + xorout[175] >= -4
 R6630: - sout[182] + xorout[174] >= 0
 R6631: - sout[174] + xorout[174] >= 0
 R6632: - sout[150] + xorout[174] >= 0
 R6633: - sout[142] + xorout[174] >= 0
 R6634: - sout[134] + xorout[174] >= 0
 R6635: - sout[134] + sout[135] >= 0
 R6636: - sout[142] + sout[143] >= 0
 R6637: - sout[150] + sout[151] >= 0
 R6638: - sout[174] + sout[175] >= 0
 R6639: - sout[182] + sout[183] >= 0
 R6640: sout[134] + sout[142] + sout[150] + sout[174] + sout[182]
   - xorout[174] >= 0
 R6641: - xorout[174] + xorout[175] >= 0
 R6642: sout[135] - sout[143] - sout[151] - sout[175] - sout[183]
   + xorout[174] - xorout[175] >= -4
 R6643: - sout[135] + sout[143] - sout[151] - sout[175] - sout[183]
   + xorout[174] - xorout[175] >= -4
 R6644: - sout[135] - sout[143] + sout[151] - sout[175] - sout[183]
   + xorout[174] - xorout[175] >= -4
 R6645: - sout[135] - sout[143] - sout[151] + sout[175] - sout[183]
   + xorout[174] - xorout[175] >= -4
 R6646: - sout[135] - sout[143] - sout[151] - sout[175] + sout[183]
   + xorout[174] - xorout[175] >= -4
 R6647: sout[135] + sout[143] + sout[151] - sout[175] - sout[183]
   + xorout[174] - xorout[175] >= -2
 R6648: sout[135] + sout[143] - sout[151] + sout[175] - sout[183]
   + xorout[174] - xorout[175] >= -2
 R6649: sout[135] - sout[143] + sout[151] + sout[175] - sout[183]
   + xorout[174] - xorout[175] >= -2
 R6650: - sout[135] + sout[143] + sout[151] + sout[175] - sout[183]
   + xorout[174] - xorout[175] >= -2
 R6651: sout[135] + sout[143] - sout[151] - sout[175] + sout[183]
   + xorout[174] - xorout[175] >= -2
 R6652: sout[135] - sout[143] + sout[151] - sout[175] + sout[183]
   + xorout[174] - xorout[175] >= -2
 R6653: - sout[135] + sout[143] + sout[151] - sout[175] + sout[183]
   + xorout[174] - xorout[175] >= -2
 R6654: sout[135] - sout[143] - sout[151] + sout[175] + sout[183]
   + xorout[174] - xorout[175] >= -2
 R6655: - sout[135] + sout[143] - sout[151] + sout[175] + sout[183]
   + xorout[174] - xorout[175] >= -2
 R6656: - sout[135] - sout[143] + sout[151] + sout[175] + sout[183]
   + xorout[174] - xorout[175] >= -2
 R6657: sout[137] + sout[145] + sout[153] + sout[177] - sout[185]
   + xorout[177] >= 0
 R6658: sout[137] + sout[145] + sout[153] - sout[177] + sout[185]
   + xorout[177] >= 0
 R6659: sout[137] + sout[145] - sout[153] + sout[177] + sout[185]
   + xorout[177] >= 0
 R6660: sout[137] - sout[145] + sout[153] + sout[177] + sout[185]
   + xorout[177] >= 0
 R6661: - sout[137] + sout[145] + sout[153] + sout[177] + sout[185]
   + xorout[177] >= 0
 R6662: sout[137] + sout[145] + sout[153] + sout[177] + sout[185]
   - xorout[177] >= 0
 R6663: sout[137] + sout[145] - sout[153] - sout[177] - sout[185]
   + xorout[177] >= -2
 R6664: sout[137] - sout[145] + sout[153] - sout[177] - sout[185]
   + xorout[177] >= -2
 R6665: - sout[137] + sout[145] + sout[153] - sout[177] - sout[185]
   + xorout[177] >= -2
 R6666: sout[137] - sout[145] - sout[153] + sout[177] - sout[185]
   + xorout[177] >= -2
 R6667: - sout[137] + sout[145] - sout[153] + sout[177] - sout[185]
   + xorout[177] >= -2
 R6668: - sout[137] - sout[145] + sout[153] + sout[177] - sout[185]
   + xorout[177] >= -2
 R6669: sout[137] - sout[145] - sout[153] - sout[177] + sout[185]
   + xorout[177] >= -2
 R6670: - sout[137] + sout[145] - sout[153] - sout[177] + sout[185]
   + xorout[177] >= -2
 R6671: - sout[137] - sout[145] + sout[153] - sout[177] + sout[185]
   + xorout[177] >= -2
 R6672: - sout[137] - sout[145] - sout[153] + sout[177] + sout[185]
   + xorout[177] >= -2
 R6673: - sout[137] - sout[145] - sout[153] - sout[177] - sout[185]
   + xorout[177] >= -4
 R6674: - sout[184] + xorout[176] >= 0
 R6675: - sout[176] + xorout[176] >= 0
 R6676: - sout[152] + xorout[176] >= 0
 R6677: - sout[144] + xorout[176] >= 0
 R6678: - sout[136] + xorout[176] >= 0
 R6679: - sout[136] + sout[137] >= 0
 R6680: - sout[144] + sout[145] >= 0
 R6681: - sout[152] + sout[153] >= 0
 R6682: - sout[176] + sout[177] >= 0
 R6683: - sout[184] + sout[185] >= 0
 R6684: sout[136] + sout[144] + sout[152] + sout[176] + sout[184]
   - xorout[176] >= 0
 R6685: - xorout[176] + xorout[177] >= 0
 R6686: sout[137] - sout[145] - sout[153] - sout[177] - sout[185]
   + xorout[176] - xorout[177] >= -4
 R6687: - sout[137] + sout[145] - sout[153] - sout[177] - sout[185]
   + xorout[176] - xorout[177] >= -4
 R6688: - sout[137] - sout[145] + sout[153] - sout[177] - sout[185]
   + xorout[176] - xorout[177] >= -4
 R6689: - sout[137] - sout[145] - sout[153] + sout[177] - sout[185]
   + xorout[176] - xorout[177] >= -4
 R6690: - sout[137] - sout[145] - sout[153] - sout[177] + sout[185]
   + xorout[176] - xorout[177] >= -4
 R6691: sout[137] + sout[145] + sout[153] - sout[177] - sout[185]
   + xorout[176] - xorout[177] >= -2
 R6692: sout[137] + sout[145] - sout[153] + sout[177] - sout[185]
   + xorout[176] - xorout[177] >= -2
 R6693: sout[137] - sout[145] + sout[153] + sout[177] - sout[185]
   + xorout[176] - xorout[177] >= -2
 R6694: - sout[137] + sout[145] + sout[153] + sout[177] - sout[185]
   + xorout[176] - xorout[177] >= -2
 R6695: sout[137] + sout[145] - sout[153] - sout[177] + sout[185]
   + xorout[176] - xorout[177] >= -2
 R6696: sout[137] - sout[145] + sout[153] - sout[177] + sout[185]
   + xorout[176] - xorout[177] >= -2
 R6697: - sout[137] + sout[145] + sout[153] - sout[177] + sout[185]
   + xorout[176] - xorout[177] >= -2
 R6698: sout[137] - sout[145] - sout[153] + sout[177] + sout[185]
   + xorout[176] - xorout[177] >= -2
 R6699: - sout[137] + sout[145] - sout[153] + sout[177] + sout[185]
   + xorout[176] - xorout[177] >= -2
 R6700: - sout[137] - sout[145] + sout[153] + sout[177] + sout[185]
   + xorout[176] - xorout[177] >= -2
 R6701: sout[139] + sout[147] + sout[155] + sout[179] - sout[187]
   + xorout[179] >= 0
 R6702: sout[139] + sout[147] + sout[155] - sout[179] + sout[187]
   + xorout[179] >= 0
 R6703: sout[139] + sout[147] - sout[155] + sout[179] + sout[187]
   + xorout[179] >= 0
 R6704: sout[139] - sout[147] + sout[155] + sout[179] + sout[187]
   + xorout[179] >= 0
 R6705: - sout[139] + sout[147] + sout[155] + sout[179] + sout[187]
   + xorout[179] >= 0
 R6706: sout[139] + sout[147] + sout[155] + sout[179] + sout[187]
   - xorout[179] >= 0
 R6707: sout[139] + sout[147] - sout[155] - sout[179] - sout[187]
   + xorout[179] >= -2
 R6708: sout[139] - sout[147] + sout[155] - sout[179] - sout[187]
   + xorout[179] >= -2
 R6709: - sout[139] + sout[147] + sout[155] - sout[179] - sout[187]
   + xorout[179] >= -2
 R6710: sout[139] - sout[147] - sout[155] + sout[179] - sout[187]
   + xorout[179] >= -2
 R6711: - sout[139] + sout[147] - sout[155] + sout[179] - sout[187]
   + xorout[179] >= -2
 R6712: - sout[139] - sout[147] + sout[155] + sout[179] - sout[187]
   + xorout[179] >= -2
 R6713: sout[139] - sout[147] - sout[155] - sout[179] + sout[187]
   + xorout[179] >= -2
 R6714: - sout[139] + sout[147] - sout[155] - sout[179] + sout[187]
   + xorout[179] >= -2
 R6715: - sout[139] - sout[147] + sout[155] - sout[179] + sout[187]
   + xorout[179] >= -2
 R6716: - sout[139] - sout[147] - sout[155] + sout[179] + sout[187]
   + xorout[179] >= -2
 R6717: - sout[139] - sout[147] - sout[155] - sout[179] - sout[187]
   + xorout[179] >= -4
 R6718: - sout[186] + xorout[178] >= 0
 R6719: - sout[178] + xorout[178] >= 0
 R6720: - sout[154] + xorout[178] >= 0
 R6721: - sout[146] + xorout[178] >= 0
 R6722: - sout[138] + xorout[178] >= 0
 R6723: - sout[138] + sout[139] >= 0
 R6724: - sout[146] + sout[147] >= 0
 R6725: - sout[154] + sout[155] >= 0
 R6726: - sout[178] + sout[179] >= 0
 R6727: - sout[186] + sout[187] >= 0
 R6728: sout[138] + sout[146] + sout[154] + sout[178] + sout[186]
   - xorout[178] >= 0
 R6729: - xorout[178] + xorout[179] >= 0
 R6730: sout[139] - sout[147] - sout[155] - sout[179] - sout[187]
   + xorout[178] - xorout[179] >= -4
 R6731: - sout[139] + sout[147] - sout[155] - sout[179] - sout[187]
   + xorout[178] - xorout[179] >= -4
 R6732: - sout[139] - sout[147] + sout[155] - sout[179] - sout[187]
   + xorout[178] - xorout[179] >= -4
 R6733: - sout[139] - sout[147] - sout[155] + sout[179] - sout[187]
   + xorout[178] - xorout[179] >= -4
 R6734: - sout[139] - sout[147] - sout[155] - sout[179] + sout[187]
   + xorout[178] - xorout[179] >= -4
 R6735: sout[139] + sout[147] + sout[155] - sout[179] - sout[187]
   + xorout[178] - xorout[179] >= -2
 R6736: sout[139] + sout[147] - sout[155] + sout[179] - sout[187]
   + xorout[178] - xorout[179] >= -2
 R6737: sout[139] - sout[147] + sout[155] + sout[179] - sout[187]
   + xorout[178] - xorout[179] >= -2
 R6738: - sout[139] + sout[147] + sout[155] + sout[179] - sout[187]
   + xorout[178] - xorout[179] >= -2
 R6739: sout[139] + sout[147] - sout[155] - sout[179] + sout[187]
   + xorout[178] - xorout[179] >= -2
 R6740: sout[139] - sout[147] + sout[155] - sout[179] + sout[187]
   + xorout[178] - xorout[179] >= -2
 R6741: - sout[139] + sout[147] + sout[155] - sout[179] + sout[187]
   + xorout[178] - xorout[179] >= -2
 R6742: sout[139] - sout[147] - sout[155] + sout[179] + sout[187]
   + xorout[178] - xorout[179] >= -2
 R6743: - sout[139] + sout[147] - sout[155] + sout[179] + sout[187]
   + xorout[178] - xorout[179] >= -2
 R6744: - sout[139] - sout[147] + sout[155] + sout[179] + sout[187]
   + xorout[178] - xorout[179] >= -2
 R6745: sout[141] + sout[149] + sout[157] + sout[181] - sout[189]
   + xorout[181] >= 0
 R6746: sout[141] + sout[149] + sout[157] - sout[181] + sout[189]
   + xorout[181] >= 0
 R6747: sout[141] + sout[149] - sout[157] + sout[181] + sout[189]
   + xorout[181] >= 0
 R6748: sout[141] - sout[149] + sout[157] + sout[181] + sout[189]
   + xorout[181] >= 0
 R6749: - sout[141] + sout[149] + sout[157] + sout[181] + sout[189]
   + xorout[181] >= 0
 R6750: sout[141] + sout[149] + sout[157] + sout[181] + sout[189]
   - xorout[181] >= 0
 R6751: sout[141] + sout[149] - sout[157] - sout[181] - sout[189]
   + xorout[181] >= -2
 R6752: sout[141] - sout[149] + sout[157] - sout[181] - sout[189]
   + xorout[181] >= -2
 R6753: - sout[141] + sout[149] + sout[157] - sout[181] - sout[189]
   + xorout[181] >= -2
 R6754: sout[141] - sout[149] - sout[157] + sout[181] - sout[189]
   + xorout[181] >= -2
 R6755: - sout[141] + sout[149] - sout[157] + sout[181] - sout[189]
   + xorout[181] >= -2
 R6756: - sout[141] - sout[149] + sout[157] + sout[181] - sout[189]
   + xorout[181] >= -2
 R6757: sout[141] - sout[149] - sout[157] - sout[181] + sout[189]
   + xorout[181] >= -2
 R6758: - sout[141] + sout[149] - sout[157] - sout[181] + sout[189]
   + xorout[181] >= -2
 R6759: - sout[141] - sout[149] + sout[157] - sout[181] + sout[189]
   + xorout[181] >= -2
 R6760: - sout[141] - sout[149] - sout[157] + sout[181] + sout[189]
   + xorout[181] >= -2
 R6761: - sout[141] - sout[149] - sout[157] - sout[181] - sout[189]
   + xorout[181] >= -4
 R6762: - sout[188] + xorout[180] >= 0
 R6763: - sout[180] + xorout[180] >= 0
 R6764: - sout[156] + xorout[180] >= 0
 R6765: - sout[148] + xorout[180] >= 0
 R6766: - sout[140] + xorout[180] >= 0
 R6767: - sout[140] + sout[141] >= 0
 R6768: - sout[148] + sout[149] >= 0
 R6769: - sout[156] + sout[157] >= 0
 R6770: - sout[180] + sout[181] >= 0
 R6771: - sout[188] + sout[189] >= 0
 R6772: sout[140] + sout[148] + sout[156] + sout[180] + sout[188]
   - xorout[180] >= 0
 R6773: - xorout[180] + xorout[181] >= 0
 R6774: sout[141] - sout[149] - sout[157] - sout[181] - sout[189]
   + xorout[180] - xorout[181] >= -4
 R6775: - sout[141] + sout[149] - sout[157] - sout[181] - sout[189]
   + xorout[180] - xorout[181] >= -4
 R6776: - sout[141] - sout[149] + sout[157] - sout[181] - sout[189]
   + xorout[180] - xorout[181] >= -4
 R6777: - sout[141] - sout[149] - sout[157] + sout[181] - sout[189]
   + xorout[180] - xorout[181] >= -4
 R6778: - sout[141] - sout[149] - sout[157] - sout[181] + sout[189]
   + xorout[180] - xorout[181] >= -4
 R6779: sout[141] + sout[149] + sout[157] - sout[181] - sout[189]
   + xorout[180] - xorout[181] >= -2
 R6780: sout[141] + sout[149] - sout[157] + sout[181] - sout[189]
   + xorout[180] - xorout[181] >= -2
 R6781: sout[141] - sout[149] + sout[157] + sout[181] - sout[189]
   + xorout[180] - xorout[181] >= -2
 R6782: - sout[141] + sout[149] + sout[157] + sout[181] - sout[189]
   + xorout[180] - xorout[181] >= -2
 R6783: sout[141] + sout[149] - sout[157] - sout[181] + sout[189]
   + xorout[180] - xorout[181] >= -2
 R6784: sout[141] - sout[149] + sout[157] - sout[181] + sout[189]
   + xorout[180] - xorout[181] >= -2
 R6785: - sout[141] + sout[149] + sout[157] - sout[181] + sout[189]
   + xorout[180] - xorout[181] >= -2
 R6786: sout[141] - sout[149] - sout[157] + sout[181] + sout[189]
   + xorout[180] - xorout[181] >= -2
 R6787: - sout[141] + sout[149] - sout[157] + sout[181] + sout[189]
   + xorout[180] - xorout[181] >= -2
 R6788: - sout[141] - sout[149] + sout[157] + sout[181] + sout[189]
   + xorout[180] - xorout[181] >= -2
 R6789: sout[143] + sout[151] + sout[159] + sout[183] - sout[191]
   + xorout[183] >= 0
 R6790: sout[143] + sout[151] + sout[159] - sout[183] + sout[191]
   + xorout[183] >= 0
 R6791: sout[143] + sout[151] - sout[159] + sout[183] + sout[191]
   + xorout[183] >= 0
 R6792: sout[143] - sout[151] + sout[159] + sout[183] + sout[191]
   + xorout[183] >= 0
 R6793: - sout[143] + sout[151] + sout[159] + sout[183] + sout[191]
   + xorout[183] >= 0
 R6794: sout[143] + sout[151] + sout[159] + sout[183] + sout[191]
   - xorout[183] >= 0
 R6795: sout[143] + sout[151] - sout[159] - sout[183] - sout[191]
   + xorout[183] >= -2
 R6796: sout[143] - sout[151] + sout[159] - sout[183] - sout[191]
   + xorout[183] >= -2
 R6797: - sout[143] + sout[151] + sout[159] - sout[183] - sout[191]
   + xorout[183] >= -2
 R6798: sout[143] - sout[151] - sout[159] + sout[183] - sout[191]
   + xorout[183] >= -2
 R6799: - sout[143] + sout[151] - sout[159] + sout[183] - sout[191]
   + xorout[183] >= -2
 R6800: - sout[143] - sout[151] + sout[159] + sout[183] - sout[191]
   + xorout[183] >= -2
 R6801: sout[143] - sout[151] - sout[159] - sout[183] + sout[191]
   + xorout[183] >= -2
 R6802: - sout[143] + sout[151] - sout[159] - sout[183] + sout[191]
   + xorout[183] >= -2
 R6803: - sout[143] - sout[151] + sout[159] - sout[183] + sout[191]
   + xorout[183] >= -2
 R6804: - sout[143] - sout[151] - sout[159] + sout[183] + sout[191]
   + xorout[183] >= -2
 R6805: - sout[143] - sout[151] - sout[159] - sout[183] - sout[191]
   + xorout[183] >= -4
 R6806: - sout[190] + xorout[182] >= 0
 R6807: - sout[182] + xorout[182] >= 0
 R6808: - sout[158] + xorout[182] >= 0
 R6809: - sout[150] + xorout[182] >= 0
 R6810: - sout[142] + xorout[182] >= 0
 R6811: - sout[142] + sout[143] >= 0
 R6812: - sout[150] + sout[151] >= 0
 R6813: - sout[158] + sout[159] >= 0
 R6814: - sout[182] + sout[183] >= 0
 R6815: - sout[190] + sout[191] >= 0
 R6816: sout[142] + sout[150] + sout[158] + sout[182] + sout[190]
   - xorout[182] >= 0
 R6817: - xorout[182] + xorout[183] >= 0
 R6818: sout[143] - sout[151] - sout[159] - sout[183] - sout[191]
   + xorout[182] - xorout[183] >= -4
 R6819: - sout[143] + sout[151] - sout[159] - sout[183] - sout[191]
   + xorout[182] - xorout[183] >= -4
 R6820: - sout[143] - sout[151] + sout[159] - sout[183] - sout[191]
   + xorout[182] - xorout[183] >= -4
 R6821: - sout[143] - sout[151] - sout[159] + sout[183] - sout[191]
   + xorout[182] - xorout[183] >= -4
 R6822: - sout[143] - sout[151] - sout[159] - sout[183] + sout[191]
   + xorout[182] - xorout[183] >= -4
 R6823: sout[143] + sout[151] + sout[159] - sout[183] - sout[191]
   + xorout[182] - xorout[183] >= -2
 R6824: sout[143] + sout[151] - sout[159] + sout[183] - sout[191]
   + xorout[182] - xorout[183] >= -2
 R6825: sout[143] - sout[151] + sout[159] + sout[183] - sout[191]
   + xorout[182] - xorout[183] >= -2
 R6826: - sout[143] + sout[151] + sout[159] + sout[183] - sout[191]
   + xorout[182] - xorout[183] >= -2
 R6827: sout[143] + sout[151] - sout[159] - sout[183] + sout[191]
   + xorout[182] - xorout[183] >= -2
 R6828: sout[143] - sout[151] + sout[159] - sout[183] + sout[191]
   + xorout[182] - xorout[183] >= -2
 R6829: - sout[143] + sout[151] + sout[159] - sout[183] + sout[191]
   + xorout[182] - xorout[183] >= -2
 R6830: sout[143] - sout[151] - sout[159] + sout[183] + sout[191]
   + xorout[182] - xorout[183] >= -2
 R6831: - sout[143] + sout[151] - sout[159] + sout[183] + sout[191]
   + xorout[182] - xorout[183] >= -2
 R6832: - sout[143] - sout[151] + sout[159] + sout[183] + sout[191]
   + xorout[182] - xorout[183] >= -2
 R6833: sout[129] + sout[145] + sout[153] + sout[161] - sout[185]
   + xorout[185] >= 0
 R6834: sout[129] + sout[145] + sout[153] - sout[161] + sout[185]
   + xorout[185] >= 0
 R6835: sout[129] + sout[145] - sout[153] + sout[161] + sout[185]
   + xorout[185] >= 0
 R6836: sout[129] - sout[145] + sout[153] + sout[161] + sout[185]
   + xorout[185] >= 0
 R6837: - sout[129] + sout[145] + sout[153] + sout[161] + sout[185]
   + xorout[185] >= 0
 R6838: sout[129] + sout[145] + sout[153] + sout[161] + sout[185]
   - xorout[185] >= 0
 R6839: sout[129] + sout[145] - sout[153] - sout[161] - sout[185]
   + xorout[185] >= -2
 R6840: sout[129] - sout[145] + sout[153] - sout[161] - sout[185]
   + xorout[185] >= -2
 R6841: - sout[129] + sout[145] + sout[153] - sout[161] - sout[185]
   + xorout[185] >= -2
 R6842: sout[129] - sout[145] - sout[153] + sout[161] - sout[185]
   + xorout[185] >= -2
 R6843: - sout[129] + sout[145] - sout[153] + sout[161] - sout[185]
   + xorout[185] >= -2
 R6844: - sout[129] - sout[145] + sout[153] + sout[161] - sout[185]
   + xorout[185] >= -2
 R6845: sout[129] - sout[145] - sout[153] - sout[161] + sout[185]
   + xorout[185] >= -2
 R6846: - sout[129] + sout[145] - sout[153] - sout[161] + sout[185]
   + xorout[185] >= -2
 R6847: - sout[129] - sout[145] + sout[153] - sout[161] + sout[185]
   + xorout[185] >= -2
 R6848: - sout[129] - sout[145] - sout[153] + sout[161] + sout[185]
   + xorout[185] >= -2
 R6849: - sout[129] - sout[145] - sout[153] - sout[161] - sout[185]
   + xorout[185] >= -4
 R6850: - sout[184] + xorout[184] >= 0
 R6851: - sout[160] + xorout[184] >= 0
 R6852: - sout[152] + xorout[184] >= 0
 R6853: - sout[144] + xorout[184] >= 0
 R6854: - sout[128] + xorout[184] >= 0
 R6855: - sout[128] + sout[129] >= 0
 R6856: - sout[144] + sout[145] >= 0
 R6857: - sout[152] + sout[153] >= 0
 R6858: - sout[160] + sout[161] >= 0
 R6859: - sout[184] + sout[185] >= 0
 R6860: sout[128] + sout[144] + sout[152] + sout[160] + sout[184]
   - xorout[184] >= 0
 R6861: - xorout[184] + xorout[185] >= 0
 R6862: sout[129] - sout[145] - sout[153] - sout[161] - sout[185]
   + xorout[184] - xorout[185] >= -4
 R6863: - sout[129] + sout[145] - sout[153] - sout[161] - sout[185]
   + xorout[184] - xorout[185] >= -4
 R6864: - sout[129] - sout[145] + sout[153] - sout[161] - sout[185]
   + xorout[184] - xorout[185] >= -4
 R6865: - sout[129] - sout[145] - sout[153] + sout[161] - sout[185]
   + xorout[184] - xorout[185] >= -4
 R6866: - sout[129] - sout[145] - sout[153] - sout[161] + sout[185]
   + xorout[184] - xorout[185] >= -4
 R6867: sout[129] + sout[145] + sout[153] - sout[161] - sout[185]
   + xorout[184] - xorout[185] >= -2
 R6868: sout[129] + sout[145] - sout[153] + sout[161] - sout[185]
   + xorout[184] - xorout[185] >= -2
 R6869: sout[129] - sout[145] + sout[153] + sout[161] - sout[185]
   + xorout[184] - xorout[185] >= -2
 R6870: - sout[129] + sout[145] + sout[153] + sout[161] - sout[185]
   + xorout[184] - xorout[185] >= -2
 R6871: sout[129] + sout[145] - sout[153] - sout[161] + sout[185]
   + xorout[184] - xorout[185] >= -2
 R6872: sout[129] - sout[145] + sout[153] - sout[161] + sout[185]
   + xorout[184] - xorout[185] >= -2
 R6873: - sout[129] + sout[145] + sout[153] - sout[161] + sout[185]
   + xorout[184] - xorout[185] >= -2
 R6874: sout[129] - sout[145] - sout[153] + sout[161] + sout[185]
   + xorout[184] - xorout[185] >= -2
 R6875: - sout[129] + sout[145] - sout[153] + sout[161] + sout[185]
   + xorout[184] - xorout[185] >= -2
 R6876: - sout[129] - sout[145] + sout[153] + sout[161] + sout[185]
   + xorout[184] - xorout[185] >= -2
 R6877: sout[131] + sout[147] + sout[155] + sout[163] - sout[187]
   + xorout[187] >= 0
 R6878: sout[131] + sout[147] + sout[155] - sout[163] + sout[187]
   + xorout[187] >= 0
 R6879: sout[131] + sout[147] - sout[155] + sout[163] + sout[187]
   + xorout[187] >= 0
 R6880: sout[131] - sout[147] + sout[155] + sout[163] + sout[187]
   + xorout[187] >= 0
 R6881: - sout[131] + sout[147] + sout[155] + sout[163] + sout[187]
   + xorout[187] >= 0
 R6882: sout[131] + sout[147] + sout[155] + sout[163] + sout[187]
   - xorout[187] >= 0
 R6883: sout[131] + sout[147] - sout[155] - sout[163] - sout[187]
   + xorout[187] >= -2
 R6884: sout[131] - sout[147] + sout[155] - sout[163] - sout[187]
   + xorout[187] >= -2
 R6885: - sout[131] + sout[147] + sout[155] - sout[163] - sout[187]
   + xorout[187] >= -2
 R6886: sout[131] - sout[147] - sout[155] + sout[163] - sout[187]
   + xorout[187] >= -2
 R6887: - sout[131] + sout[147] - sout[155] + sout[163] - sout[187]
   + xorout[187] >= -2
 R6888: - sout[131] - sout[147] + sout[155] + sout[163] - sout[187]
   + xorout[187] >= -2
 R6889: sout[131] - sout[147] - sout[155] - sout[163] + sout[187]
   + xorout[187] >= -2
 R6890: - sout[131] + sout[147] - sout[155] - sout[163] + sout[187]
   + xorout[187] >= -2
 R6891: - sout[131] - sout[147] + sout[155] - sout[163] + sout[187]
   + xorout[187] >= -2
 R6892: - sout[131] - sout[147] - sout[155] + sout[163] + sout[187]
   + xorout[187] >= -2
 R6893: - sout[131] - sout[147] - sout[155] - sout[163] - sout[187]
   + xorout[187] >= -4
 R6894: - sout[186] + xorout[186] >= 0
 R6895: - sout[162] + xorout[186] >= 0
 R6896: - sout[154] + xorout[186] >= 0
 R6897: - sout[146] + xorout[186] >= 0
 R6898: - sout[130] + xorout[186] >= 0
 R6899: - sout[130] + sout[131] >= 0
 R6900: - sout[146] + sout[147] >= 0
 R6901: - sout[154] + sout[155] >= 0
 R6902: - sout[162] + sout[163] >= 0
 R6903: - sout[186] + sout[187] >= 0
 R6904: sout[130] + sout[146] + sout[154] + sout[162] + sout[186]
   - xorout[186] >= 0
 R6905: - xorout[186] + xorout[187] >= 0
 R6906: sout[131] - sout[147] - sout[155] - sout[163] - sout[187]
   + xorout[186] - xorout[187] >= -4
 R6907: - sout[131] + sout[147] - sout[155] - sout[163] - sout[187]
   + xorout[186] - xorout[187] >= -4
 R6908: - sout[131] - sout[147] + sout[155] - sout[163] - sout[187]
   + xorout[186] - xorout[187] >= -4
 R6909: - sout[131] - sout[147] - sout[155] + sout[163] - sout[187]
   + xorout[186] - xorout[187] >= -4
 R6910: - sout[131] - sout[147] - sout[155] - sout[163] + sout[187]
   + xorout[186] - xorout[187] >= -4
 R6911: sout[131] + sout[147] + sout[155] - sout[163] - sout[187]
   + xorout[186] - xorout[187] >= -2
 R6912: sout[131] + sout[147] - sout[155] + sout[163] - sout[187]
   + xorout[186] - xorout[187] >= -2
 R6913: sout[131] - sout[147] + sout[155] + sout[163] - sout[187]
   + xorout[186] - xorout[187] >= -2
 R6914: - sout[131] + sout[147] + sout[155] + sout[163] - sout[187]
   + xorout[186] - xorout[187] >= -2
 R6915: sout[131] + sout[147] - sout[155] - sout[163] + sout[187]
   + xorout[186] - xorout[187] >= -2
 R6916: sout[131] - sout[147] + sout[155] - sout[163] + sout[187]
   + xorout[186] - xorout[187] >= -2
 R6917: - sout[131] + sout[147] + sout[155] - sout[163] + sout[187]
   + xorout[186] - xorout[187] >= -2
 R6918: sout[131] - sout[147] - sout[155] + sout[163] + sout[187]
   + xorout[186] - xorout[187] >= -2
 R6919: - sout[131] + sout[147] - sout[155] + sout[163] + sout[187]
   + xorout[186] - xorout[187] >= -2
 R6920: - sout[131] - sout[147] + sout[155] + sout[163] + sout[187]
   + xorout[186] - xorout[187] >= -2
 R6921: sout[133] + sout[149] + sout[157] + sout[165] - sout[189]
   + xorout[189] >= 0
 R6922: sout[133] + sout[149] + sout[157] - sout[165] + sout[189]
   + xorout[189] >= 0
 R6923: sout[133] + sout[149] - sout[157] + sout[165] + sout[189]
   + xorout[189] >= 0
 R6924: sout[133] - sout[149] + sout[157] + sout[165] + sout[189]
   + xorout[189] >= 0
 R6925: - sout[133] + sout[149] + sout[157] + sout[165] + sout[189]
   + xorout[189] >= 0
 R6926: sout[133] + sout[149] + sout[157] + sout[165] + sout[189]
   - xorout[189] >= 0
 R6927: sout[133] + sout[149] - sout[157] - sout[165] - sout[189]
   + xorout[189] >= -2
 R6928: sout[133] - sout[149] + sout[157] - sout[165] - sout[189]
   + xorout[189] >= -2
 R6929: - sout[133] + sout[149] + sout[157] - sout[165] - sout[189]
   + xorout[189] >= -2
 R6930: sout[133] - sout[149] - sout[157] + sout[165] - sout[189]
   + xorout[189] >= -2
 R6931: - sout[133] + sout[149] - sout[157] + sout[165] - sout[189]
   + xorout[189] >= -2
 R6932: - sout[133] - sout[149] + sout[157] + sout[165] - sout[189]
   + xorout[189] >= -2
 R6933: sout[133] - sout[149] - sout[157] - sout[165] + sout[189]
   + xorout[189] >= -2
 R6934: - sout[133] + sout[149] - sout[157] - sout[165] + sout[189]
   + xorout[189] >= -2
 R6935: - sout[133] - sout[149] + sout[157] - sout[165] + sout[189]
   + xorout[189] >= -2
 R6936: - sout[133] - sout[149] - sout[157] + sout[165] + sout[189]
   + xorout[189] >= -2
 R6937: - sout[133] - sout[149] - sout[157] - sout[165] - sout[189]
   + xorout[189] >= -4
 R6938: - sout[188] + xorout[188] >= 0
 R6939: - sout[164] + xorout[188] >= 0
 R6940: - sout[156] + xorout[188] >= 0
 R6941: - sout[148] + xorout[188] >= 0
 R6942: - sout[132] + xorout[188] >= 0
 R6943: - sout[132] + sout[133] >= 0
 R6944: - sout[148] + sout[149] >= 0
 R6945: - sout[156] + sout[157] >= 0
 R6946: - sout[164] + sout[165] >= 0
 R6947: - sout[188] + sout[189] >= 0
 R6948: sout[132] + sout[148] + sout[156] + sout[164] + sout[188]
   - xorout[188] >= 0
 R6949: - xorout[188] + xorout[189] >= 0
 R6950: sout[133] - sout[149] - sout[157] - sout[165] - sout[189]
   + xorout[188] - xorout[189] >= -4
 R6951: - sout[133] + sout[149] - sout[157] - sout[165] - sout[189]
   + xorout[188] - xorout[189] >= -4
 R6952: - sout[133] - sout[149] + sout[157] - sout[165] - sout[189]
   + xorout[188] - xorout[189] >= -4
 R6953: - sout[133] - sout[149] - sout[157] + sout[165] - sout[189]
   + xorout[188] - xorout[189] >= -4
 R6954: - sout[133] - sout[149] - sout[157] - sout[165] + sout[189]
   + xorout[188] - xorout[189] >= -4
 R6955: sout[133] + sout[149] + sout[157] - sout[165] - sout[189]
   + xorout[188] - xorout[189] >= -2
 R6956: sout[133] + sout[149] - sout[157] + sout[165] - sout[189]
   + xorout[188] - xorout[189] >= -2
 R6957: sout[133] - sout[149] + sout[157] + sout[165] - sout[189]
   + xorout[188] - xorout[189] >= -2
 R6958: - sout[133] + sout[149] + sout[157] + sout[165] - sout[189]
   + xorout[188] - xorout[189] >= -2
 R6959: sout[133] + sout[149] - sout[157] - sout[165] + sout[189]
   + xorout[188] - xorout[189] >= -2
 R6960: sout[133] - sout[149] + sout[157] - sout[165] + sout[189]
   + xorout[188] - xorout[189] >= -2
 R6961: - sout[133] + sout[149] + sout[157] - sout[165] + sout[189]
   + xorout[188] - xorout[189] >= -2
 R6962: sout[133] - sout[149] - sout[157] + sout[165] + sout[189]
   + xorout[188] - xorout[189] >= -2
 R6963: - sout[133] + sout[149] - sout[157] + sout[165] + sout[189]
   + xorout[188] - xorout[189] >= -2
 R6964: - sout[133] - sout[149] + sout[157] + sout[165] + sout[189]
   + xorout[188] - xorout[189] >= -2
 R6965: sout[135] + sout[151] + sout[159] + sout[167] - sout[191]
   + xorout[191] >= 0
 R6966: sout[135] + sout[151] + sout[159] - sout[167] + sout[191]
   + xorout[191] >= 0
 R6967: sout[135] + sout[151] - sout[159] + sout[167] + sout[191]
   + xorout[191] >= 0
 R6968: sout[135] - sout[151] + sout[159] + sout[167] + sout[191]
   + xorout[191] >= 0
 R6969: - sout[135] + sout[151] + sout[159] + sout[167] + sout[191]
   + xorout[191] >= 0
 R6970: sout[135] + sout[151] + sout[159] + sout[167] + sout[191]
   - xorout[191] >= 0
 R6971: sout[135] + sout[151] - sout[159] - sout[167] - sout[191]
   + xorout[191] >= -2
 R6972: sout[135] - sout[151] + sout[159] - sout[167] - sout[191]
   + xorout[191] >= -2
 R6973: - sout[135] + sout[151] + sout[159] - sout[167] - sout[191]
   + xorout[191] >= -2
 R6974: sout[135] - sout[151] - sout[159] + sout[167] - sout[191]
   + xorout[191] >= -2
 R6975: - sout[135] + sout[151] - sout[159] + sout[167] - sout[191]
   + xorout[191] >= -2
 R6976: - sout[135] - sout[151] + sout[159] + sout[167] - sout[191]
   + xorout[191] >= -2
 R6977: sout[135] - sout[151] - sout[159] - sout[167] + sout[191]
   + xorout[191] >= -2
 R6978: - sout[135] + sout[151] - sout[159] - sout[167] + sout[191]
   + xorout[191] >= -2
 R6979: - sout[135] - sout[151] + sout[159] - sout[167] + sout[191]
   + xorout[191] >= -2
 R6980: - sout[135] - sout[151] - sout[159] + sout[167] + sout[191]
   + xorout[191] >= -2
 R6981: - sout[135] - sout[151] - sout[159] - sout[167] - sout[191]
   + xorout[191] >= -4
 R6982: - sout[190] + xorout[190] >= 0
 R6983: - sout[166] + xorout[190] >= 0
 R6984: - sout[158] + xorout[190] >= 0
 R6985: - sout[150] + xorout[190] >= 0
 R6986: - sout[134] + xorout[190] >= 0
 R6987: - sout[134] + sout[135] >= 0
 R6988: - sout[150] + sout[151] >= 0
 R6989: - sout[158] + sout[159] >= 0
 R6990: - sout[166] + sout[167] >= 0
 R6991: - sout[190] + sout[191] >= 0
 R6992: sout[134] + sout[150] + sout[158] + sout[166] + sout[190]
   - xorout[190] >= 0
 R6993: - xorout[190] + xorout[191] >= 0
 R6994: sout[135] - sout[151] - sout[159] - sout[167] - sout[191]
   + xorout[190] - xorout[191] >= -4
 R6995: - sout[135] + sout[151] - sout[159] - sout[167] - sout[191]
   + xorout[190] - xorout[191] >= -4
 R6996: - sout[135] - sout[151] + sout[159] - sout[167] - sout[191]
   + xorout[190] - xorout[191] >= -4
 R6997: - sout[135] - sout[151] - sout[159] + sout[167] - sout[191]
   + xorout[190] - xorout[191] >= -4
 R6998: - sout[135] - sout[151] - sout[159] - sout[167] + sout[191]
   + xorout[190] - xorout[191] >= -4
 R6999: sout[135] + sout[151] + sout[159] - sout[167] - sout[191]
   + xorout[190] - xorout[191] >= -2
 R7000: sout[135] + sout[151] - sout[159] + sout[167] - sout[191]
   + xorout[190] - xorout[191] >= -2
 R7001: sout[135] - sout[151] + sout[159] + sout[167] - sout[191]
   + xorout[190] - xorout[191] >= -2
 R7002: - sout[135] + sout[151] + sout[159] + sout[167] - sout[191]
   + xorout[190] - xorout[191] >= -2
 R7003: sout[135] + sout[151] - sout[159] - sout[167] + sout[191]
   + xorout[190] - xorout[191] >= -2
 R7004: sout[135] - sout[151] + sout[159] - sout[167] + sout[191]
   + xorout[190] - xorout[191] >= -2
 R7005: - sout[135] + sout[151] + sout[159] - sout[167] + sout[191]
   + xorout[190] - xorout[191] >= -2
 R7006: sout[135] - sout[151] - sout[159] + sout[167] + sout[191]
   + xorout[190] - xorout[191] >= -2
 R7007: - sout[135] + sout[151] - sout[159] + sout[167] + sout[191]
   + xorout[190] - xorout[191] >= -2
 R7008: - sout[135] - sout[151] + sout[159] + sout[167] + sout[191]
   + xorout[190] - xorout[191] >= -2
 R7009: - x[321] + x[385] + xorout[129] >= 0
 R7010: x[321] + x[385] - xorout[129] >= 0
 R7011: - x[320] + x[384] >= 0
 R7012: x[384] - xorout[128] >= 0
 R7013: x[321] - x[385] + xorout[129] >= 0
 R7014: - xorout[128] + xorout[129] >= 0
 R7015: - x[320] + x[321] >= 0
 R7016: x[320] - x[384] + xorout[128] >= 0
 R7017: - x[384] + x[385] >= 0
 R7018: - x[321] + x[384] - x[385] - xorout[129] >= -2
 R7019: - x[323] + x[387] + xorout[131] >= 0
 R7020: x[323] + x[387] - xorout[131] >= 0
 R7021: - x[322] + x[386] >= 0
 R7022: x[386] - xorout[130] >= 0
 R7023: x[323] - x[387] + xorout[131] >= 0
 R7024: - xorout[130] + xorout[131] >= 0
 R7025: - x[322] + x[323] >= 0
 R7026: x[322] - x[386] + xorout[130] >= 0
 R7027: - x[386] + x[387] >= 0
 R7028: - x[323] + x[386] - x[387] - xorout[131] >= -2
 R7029: - x[325] + x[389] + xorout[133] >= 0
 R7030: x[325] + x[389] - xorout[133] >= 0
 R7031: - x[324] + x[388] >= 0
 R7032: x[388] - xorout[132] >= 0
 R7033: x[325] - x[389] + xorout[133] >= 0
 R7034: - xorout[132] + xorout[133] >= 0
 R7035: - x[324] + x[325] >= 0
 R7036: x[324] - x[388] + xorout[132] >= 0
 R7037: - x[388] + x[389] >= 0
 R7038: - x[325] + x[388] - x[389] - xorout[133] >= -2
 R7039: - x[327] + x[391] + xorout[135] >= 0
 R7040: x[327] + x[391] - xorout[135] >= 0
 R7041: - x[326] + x[390] >= 0
 R7042: x[390] - xorout[134] >= 0
 R7043: x[327] - x[391] + xorout[135] >= 0
 R7044: - xorout[134] + xorout[135] >= 0
 R7045: - x[326] + x[327] >= 0
 R7046: x[326] - x[390] + xorout[134] >= 0
 R7047: - x[390] + x[391] >= 0
 R7048: - x[327] + x[390] - x[391] - xorout[135] >= -2
 R7049: - x[329] + x[393] + xorout[137] >= 0
 R7050: x[329] + x[393] - xorout[137] >= 0
 R7051: - x[328] + x[392] >= 0
 R7052: x[392] - xorout[136] >= 0
 R7053: x[329] - x[393] + xorout[137] >= 0
 R7054: - xorout[136] + xorout[137] >= 0
 R7055: - x[328] + x[329] >= 0
 R7056: x[328] - x[392] + xorout[136] >= 0
 R7057: - x[392] + x[393] >= 0
 R7058: - x[329] + x[392] - x[393] - xorout[137] >= -2
 R7059: - x[331] + x[395] + xorout[139] >= 0
 R7060: x[331] + x[395] - xorout[139] >= 0
 R7061: - x[330] + x[394] >= 0
 R7062: x[394] - xorout[138] >= 0
 R7063: x[331] - x[395] + xorout[139] >= 0
 R7064: - xorout[138] + xorout[139] >= 0
 R7065: - x[330] + x[331] >= 0
 R7066: x[330] - x[394] + xorout[138] >= 0
 R7067: - x[394] + x[395] >= 0
 R7068: - x[331] + x[394] - x[395] - xorout[139] >= -2
 R7069: - x[333] + x[397] + xorout[141] >= 0
 R7070: x[333] + x[397] - xorout[141] >= 0
 R7071: - x[332] + x[396] >= 0
 R7072: x[396] - xorout[140] >= 0
 R7073: x[333] - x[397] + xorout[141] >= 0
 R7074: - xorout[140] + xorout[141] >= 0
 R7075: - x[332] + x[333] >= 0
 R7076: x[332] - x[396] + xorout[140] >= 0
 R7077: - x[396] + x[397] >= 0
 R7078: - x[333] + x[396] - x[397] - xorout[141] >= -2
 R7079: - x[335] + x[399] + xorout[143] >= 0
 R7080: x[335] + x[399] - xorout[143] >= 0
 R7081: - x[334] + x[398] >= 0
 R7082: x[398] - xorout[142] >= 0
 R7083: x[335] - x[399] + xorout[143] >= 0
 R7084: - xorout[142] + xorout[143] >= 0
 R7085: - x[334] + x[335] >= 0
 R7086: x[334] - x[398] + xorout[142] >= 0
 R7087: - x[398] + x[399] >= 0
 R7088: - x[335] + x[398] - x[399] - xorout[143] >= -2
 R7089: - x[337] + x[401] + xorout[145] >= 0
 R7090: x[337] + x[401] - xorout[145] >= 0
 R7091: - x[336] + x[400] >= 0
 R7092: x[400] - xorout[144] >= 0
 R7093: x[337] - x[401] + xorout[145] >= 0
 R7094: - xorout[144] + xorout[145] >= 0
 R7095: - x[336] + x[337] >= 0
 R7096: x[336] - x[400] + xorout[144] >= 0
 R7097: - x[400] + x[401] >= 0
 R7098: - x[337] + x[400] - x[401] - xorout[145] >= -2
 R7099: - x[339] + x[403] + xorout[147] >= 0
 R7100: x[339] + x[403] - xorout[147] >= 0
 R7101: - x[338] + x[402] >= 0
 R7102: x[402] - xorout[146] >= 0
 R7103: x[339] - x[403] + xorout[147] >= 0
 R7104: - xorout[146] + xorout[147] >= 0
 R7105: - x[338] + x[339] >= 0
 R7106: x[338] - x[402] + xorout[146] >= 0
 R7107: - x[402] + x[403] >= 0
 R7108: - x[339] + x[402] - x[403] - xorout[147] >= -2
 R7109: - x[341] + x[405] + xorout[149] >= 0
 R7110: x[341] + x[405] - xorout[149] >= 0
 R7111: - x[340] + x[404] >= 0
 R7112: x[404] - xorout[148] >= 0
 R7113: x[341] - x[405] + xorout[149] >= 0
 R7114: - xorout[148] + xorout[149] >= 0
 R7115: - x[340] + x[341] >= 0
 R7116: x[340] - x[404] + xorout[148] >= 0
 R7117: - x[404] + x[405] >= 0
 R7118: - x[341] + x[404] - x[405] - xorout[149] >= -2
 R7119: - x[343] + x[407] + xorout[151] >= 0
 R7120: x[343] + x[407] - xorout[151] >= 0
 R7121: - x[342] + x[406] >= 0
 R7122: x[406] - xorout[150] >= 0
 R7123: x[343] - x[407] + xorout[151] >= 0
 R7124: - xorout[150] + xorout[151] >= 0
 R7125: - x[342] + x[343] >= 0
 R7126: x[342] - x[406] + xorout[150] >= 0
 R7127: - x[406] + x[407] >= 0
 R7128: - x[343] + x[406] - x[407] - xorout[151] >= -2
 R7129: - x[345] + x[409] + xorout[153] >= 0
 R7130: x[345] + x[409] - xorout[153] >= 0
 R7131: - x[344] + x[408] >= 0
 R7132: x[408] - xorout[152] >= 0
 R7133: x[345] - x[409] + xorout[153] >= 0
 R7134: - xorout[152] + xorout[153] >= 0
 R7135: - x[344] + x[345] >= 0
 R7136: x[344] - x[408] + xorout[152] >= 0
 R7137: - x[408] + x[409] >= 0
 R7138: - x[345] + x[408] - x[409] - xorout[153] >= -2
 R7139: - x[347] + x[411] + xorout[155] >= 0
 R7140: x[347] + x[411] - xorout[155] >= 0
 R7141: - x[346] + x[410] >= 0
 R7142: x[410] - xorout[154] >= 0
 R7143: x[347] - x[411] + xorout[155] >= 0
 R7144: - xorout[154] + xorout[155] >= 0
 R7145: - x[346] + x[347] >= 0
 R7146: x[346] - x[410] + xorout[154] >= 0
 R7147: - x[410] + x[411] >= 0
 R7148: - x[347] + x[410] - x[411] - xorout[155] >= -2
 R7149: - x[349] + x[413] + xorout[157] >= 0
 R7150: x[349] + x[413] - xorout[157] >= 0
 R7151: - x[348] + x[412] >= 0
 R7152: x[412] - xorout[156] >= 0
 R7153: x[349] - x[413] + xorout[157] >= 0
 R7154: - xorout[156] + xorout[157] >= 0
 R7155: - x[348] + x[349] >= 0
 R7156: x[348] - x[412] + xorout[156] >= 0
 R7157: - x[412] + x[413] >= 0
 R7158: - x[349] + x[412] - x[413] - xorout[157] >= -2
 R7159: - x[351] + x[415] + xorout[159] >= 0
 R7160: x[351] + x[415] - xorout[159] >= 0
 R7161: - x[350] + x[414] >= 0
 R7162: x[414] - xorout[158] >= 0
 R7163: x[351] - x[415] + xorout[159] >= 0
 R7164: - xorout[158] + xorout[159] >= 0
 R7165: - x[350] + x[351] >= 0
 R7166: x[350] - x[414] + xorout[158] >= 0
 R7167: - x[414] + x[415] >= 0
 R7168: - x[351] + x[414] - x[415] - xorout[159] >= -2
 R7169: - x[353] + x[417] + xorout[161] >= 0
 R7170: x[353] + x[417] - xorout[161] >= 0
 R7171: - x[352] + x[416] >= 0
 R7172: x[416] - xorout[160] >= 0
 R7173: x[353] - x[417] + xorout[161] >= 0
 R7174: - xorout[160] + xorout[161] >= 0
 R7175: - x[352] + x[353] >= 0
 R7176: x[352] - x[416] + xorout[160] >= 0
 R7177: - x[416] + x[417] >= 0
 R7178: - x[353] + x[416] - x[417] - xorout[161] >= -2
 R7179: - x[355] + x[419] + xorout[163] >= 0
 R7180: x[355] + x[419] - xorout[163] >= 0
 R7181: - x[354] + x[418] >= 0
 R7182: x[418] - xorout[162] >= 0
 R7183: x[355] - x[419] + xorout[163] >= 0
 R7184: - xorout[162] + xorout[163] >= 0
 R7185: - x[354] + x[355] >= 0
 R7186: x[354] - x[418] + xorout[162] >= 0
 R7187: - x[418] + x[419] >= 0
 R7188: - x[355] + x[418] - x[419] - xorout[163] >= -2
 R7189: - x[357] + x[421] + xorout[165] >= 0
 R7190: x[357] + x[421] - xorout[165] >= 0
 R7191: - x[356] + x[420] >= 0
 R7192: x[420] - xorout[164] >= 0
 R7193: x[357] - x[421] + xorout[165] >= 0
 R7194: - xorout[164] + xorout[165] >= 0
 R7195: - x[356] + x[357] >= 0
 R7196: x[356] - x[420] + xorout[164] >= 0
 R7197: - x[420] + x[421] >= 0
 R7198: - x[357] + x[420] - x[421] - xorout[165] >= -2
 R7199: - x[359] + x[423] + xorout[167] >= 0
 R7200: x[359] + x[423] - xorout[167] >= 0
 R7201: - x[358] + x[422] >= 0
 R7202: x[422] - xorout[166] >= 0
 R7203: x[359] - x[423] + xorout[167] >= 0
 R7204: - xorout[166] + xorout[167] >= 0
 R7205: - x[358] + x[359] >= 0
 R7206: x[358] - x[422] + xorout[166] >= 0
 R7207: - x[422] + x[423] >= 0
 R7208: - x[359] + x[422] - x[423] - xorout[167] >= -2
 R7209: - x[361] + x[425] + xorout[169] >= 0
 R7210: x[361] + x[425] - xorout[169] >= 0
 R7211: - x[360] + x[424] >= 0
 R7212: x[424] - xorout[168] >= 0
 R7213: x[361] - x[425] + xorout[169] >= 0
 R7214: - xorout[168] + xorout[169] >= 0
 R7215: - x[360] + x[361] >= 0
 R7216: x[360] - x[424] + xorout[168] >= 0
 R7217: - x[424] + x[425] >= 0
 R7218: - x[361] + x[424] - x[425] - xorout[169] >= -2
 R7219: - x[363] + x[427] + xorout[171] >= 0
 R7220: x[363] + x[427] - xorout[171] >= 0
 R7221: - x[362] + x[426] >= 0
 R7222: x[426] - xorout[170] >= 0
 R7223: x[363] - x[427] + xorout[171] >= 0
 R7224: - xorout[170] + xorout[171] >= 0
 R7225: - x[362] + x[363] >= 0
 R7226: x[362] - x[426] + xorout[170] >= 0
 R7227: - x[426] + x[427] >= 0
 R7228: - x[363] + x[426] - x[427] - xorout[171] >= -2
 R7229: - x[365] + x[429] + xorout[173] >= 0
 R7230: x[365] + x[429] - xorout[173] >= 0
 R7231: - x[364] + x[428] >= 0
 R7232: x[428] - xorout[172] >= 0
 R7233: x[365] - x[429] + xorout[173] >= 0
 R7234: - xorout[172] + xorout[173] >= 0
 R7235: - x[364] + x[365] >= 0
 R7236: x[364] - x[428] + xorout[172] >= 0
 R7237: - x[428] + x[429] >= 0
 R7238: - x[365] + x[428] - x[429] - xorout[173] >= -2
 R7239: - x[367] + x[431] + xorout[175] >= 0
 R7240: x[367] + x[431] - xorout[175] >= 0
 R7241: - x[366] + x[430] >= 0
 R7242: x[430] - xorout[174] >= 0
 R7243: x[367] - x[431] + xorout[175] >= 0
 R7244: - xorout[174] + xorout[175] >= 0
 R7245: - x[366] + x[367] >= 0
 R7246: x[366] - x[430] + xorout[174] >= 0
 R7247: - x[430] + x[431] >= 0
 R7248: - x[367] + x[430] - x[431] - xorout[175] >= -2
 R7249: - x[369] + x[433] + xorout[177] >= 0
 R7250: x[369] + x[433] - xorout[177] >= 0
 R7251: - x[368] + x[432] >= 0
 R7252: x[432] - xorout[176] >= 0
 R7253: x[369] - x[433] + xorout[177] >= 0
 R7254: - xorout[176] + xorout[177] >= 0
 R7255: - x[368] + x[369] >= 0
 R7256: x[368] - x[432] + xorout[176] >= 0
 R7257: - x[432] + x[433] >= 0
 R7258: - x[369] + x[432] - x[433] - xorout[177] >= -2
 R7259: - x[371] + x[435] + xorout[179] >= 0
 R7260: x[371] + x[435] - xorout[179] >= 0
 R7261: - x[370] + x[434] >= 0
 R7262: x[434] - xorout[178] >= 0
 R7263: x[371] - x[435] + xorout[179] >= 0
 R7264: - xorout[178] + xorout[179] >= 0
 R7265: - x[370] + x[371] >= 0
 R7266: x[370] - x[434] + xorout[178] >= 0
 R7267: - x[434] + x[435] >= 0
 R7268: - x[371] + x[434] - x[435] - xorout[179] >= -2
 R7269: - x[373] + x[437] + xorout[181] >= 0
 R7270: x[373] + x[437] - xorout[181] >= 0
 R7271: - x[372] + x[436] >= 0
 R7272: x[436] - xorout[180] >= 0
 R7273: x[373] - x[437] + xorout[181] >= 0
 R7274: - xorout[180] + xorout[181] >= 0
 R7275: - x[372] + x[373] >= 0
 R7276: x[372] - x[436] + xorout[180] >= 0
 R7277: - x[436] + x[437] >= 0
 R7278: - x[373] + x[436] - x[437] - xorout[181] >= -2
 R7279: - x[375] + x[439] + xorout[183] >= 0
 R7280: x[375] + x[439] - xorout[183] >= 0
 R7281: - x[374] + x[438] >= 0
 R7282: x[438] - xorout[182] >= 0
 R7283: x[375] - x[439] + xorout[183] >= 0
 R7284: - xorout[182] + xorout[183] >= 0
 R7285: - x[374] + x[375] >= 0
 R7286: x[374] - x[438] + xorout[182] >= 0
 R7287: - x[438] + x[439] >= 0
 R7288: - x[375] + x[438] - x[439] - xorout[183] >= -2
 R7289: - x[377] + x[441] + xorout[185] >= 0
 R7290: x[377] + x[441] - xorout[185] >= 0
 R7291: - x[376] + x[440] >= 0
 R7292: x[440] - xorout[184] >= 0
 R7293: x[377] - x[441] + xorout[185] >= 0
 R7294: - xorout[184] + xorout[185] >= 0
 R7295: - x[376] + x[377] >= 0
 R7296: x[376] - x[440] + xorout[184] >= 0
 R7297: - x[440] + x[441] >= 0
 R7298: - x[377] + x[440] - x[441] - xorout[185] >= -2
 R7299: - x[379] + x[443] + xorout[187] >= 0
 R7300: x[379] + x[443] - xorout[187] >= 0
 R7301: - x[378] + x[442] >= 0
 R7302: x[442] - xorout[186] >= 0
 R7303: x[379] - x[443] + xorout[187] >= 0
 R7304: - xorout[186] + xorout[187] >= 0
 R7305: - x[378] + x[379] >= 0
 R7306: x[378] - x[442] + xorout[186] >= 0
 R7307: - x[442] + x[443] >= 0
 R7308: - x[379] + x[442] - x[443] - xorout[187] >= -2
 R7309: - x[381] + x[445] + xorout[189] >= 0
 R7310: x[381] + x[445] - xorout[189] >= 0
 R7311: - x[380] + x[444] >= 0
 R7312: x[444] - xorout[188] >= 0
 R7313: x[381] - x[445] + xorout[189] >= 0
 R7314: - xorout[188] + xorout[189] >= 0
 R7315: - x[380] + x[381] >= 0
 R7316: x[380] - x[444] + xorout[188] >= 0
 R7317: - x[444] + x[445] >= 0
 R7318: - x[381] + x[444] - x[445] - xorout[189] >= -2
 R7319: - x[383] + x[447] + xorout[191] >= 0
 R7320: x[383] + x[447] - xorout[191] >= 0
 R7321: - x[382] + x[446] >= 0
 R7322: x[446] - xorout[190] >= 0
 R7323: x[383] - x[447] + xorout[191] >= 0
 R7324: - xorout[190] + xorout[191] >= 0
 R7325: - x[382] + x[383] >= 0
 R7326: x[382] - x[446] + xorout[190] >= 0
 R7327: - x[446] + x[447] >= 0
 R7328: - x[383] + x[446] - x[447] - xorout[191] >= -2
 R7329: x[256] - x[448] = 0
 R7330: x[257] - x[449] = 0
 R7331: x[258] - x[450] = 0
 R7332: x[259] - x[451] = 0
 R7333: x[260] - x[452] = 0
 R7334: x[261] - x[453] = 0
 R7335: x[262] - x[454] = 0
 R7336: x[263] - x[455] = 0
 R7337: x[264] - x[456] = 0
 R7338: x[265] - x[457] = 0
 R7339: x[266] - x[458] = 0
 R7340: x[267] - x[459] = 0
 R7341: x[268] - x[460] = 0
 R7342: x[269] - x[461] = 0
 R7343: x[270] - x[462] = 0
 R7344: x[271] - x[463] = 0
 R7345: x[272] - x[464] = 0
 R7346: x[273] - x[465] = 0
 R7347: x[274] - x[466] = 0
 R7348: x[275] - x[467] = 0
 R7349: x[276] - x[468] = 0
 R7350: x[277] - x[469] = 0
 R7351: x[278] - x[470] = 0
 R7352: x[279] - x[471] = 0
 R7353: x[280] - x[472] = 0
 R7354: x[281] - x[473] = 0
 R7355: x[282] - x[474] = 0
 R7356: x[283] - x[475] = 0
 R7357: x[284] - x[476] = 0
 R7358: x[285] - x[477] = 0
 R7359: x[286] - x[478] = 0
 R7360: x[287] - x[479] = 0
 R7361: x[288] - x[480] = 0
 R7362: x[289] - x[481] = 0
 R7363: x[290] - x[482] = 0
 R7364: x[291] - x[483] = 0
 R7365: x[292] - x[484] = 0
 R7366: x[293] - x[485] = 0
 R7367: x[294] - x[486] = 0
 R7368: x[295] - x[487] = 0
 R7369: x[296] - x[488] = 0
 R7370: x[297] - x[489] = 0
 R7371: x[298] - x[490] = 0
 R7372: x[299] - x[491] = 0
 R7373: x[300] - x[492] = 0
 R7374: x[301] - x[493] = 0
 R7375: x[302] - x[494] = 0
 R7376: x[303] - x[495] = 0
 R7377: x[304] - x[496] = 0
 R7378: x[305] - x[497] = 0
 R7379: x[306] - x[498] = 0
 R7380: x[307] - x[499] = 0
 R7381: x[308] - x[500] = 0
 R7382: x[309] - x[501] = 0
 R7383: x[310] - x[502] = 0
 R7384: x[311] - x[503] = 0
 R7385: x[312] - x[504] = 0
 R7386: x[313] - x[505] = 0
 R7387: x[314] - x[506] = 0
 R7388: x[315] - x[507] = 0
 R7389: x[316] - x[508] = 0
 R7390: x[317] - x[509] = 0
 R7391: x[318] - x[510] = 0
 R7392: x[319] - x[511] = 0
 R7393: x[510] = 0
 R7394: x[511] = 1
Bounds
Binaries
 x[0] x[1] x[2] x[3] x[4] x[5] x[6] x[7] x[8] x[9] x[10] x[11] x[12] x[13]
 x[14] x[15] x[16] x[17] x[18] x[19] x[20] x[21] x[22] x[23] x[24] x[25]
 x[26] x[27] x[28] x[29] x[30] x[31] x[32] x[33] x[34] x[35] x[36] x[37]
 x[38] x[39] x[40] x[41] x[42] x[43] x[44] x[45] x[46] x[47] x[48] x[49]
 x[50] x[51] x[52] x[53] x[54] x[55] x[56] x[57] x[58] x[59] x[60] x[61]
 x[62] x[63] x[64] x[65] x[66] x[67] x[68] x[69] x[70] x[71] x[72] x[73]
 x[74] x[75] x[76] x[77] x[78] x[79] x[80] x[81] x[82] x[83] x[84] x[85]
 x[86] x[87] x[88] x[89] x[90] x[91] x[92] x[93] x[94] x[95] x[96] x[97]
 x[98] x[99] x[100] x[101] x[102] x[103] x[104] x[105] x[106] x[107] x[108]
 x[109] x[110] x[111] x[112] x[113] x[114] x[115] x[116] x[117] x[118]
 x[119] x[120] x[121] x[122] x[123] x[124] x[125] x[126] x[127] x[128]
 x[129] x[130] x[131] x[132] x[133] x[134] x[135] x[136] x[137] x[138]
 x[139] x[140] x[141] x[142] x[143] x[144] x[145] x[146] x[147] x[148]
 x[149] x[150] x[151] x[152] x[153] x[154] x[155] x[156] x[157] x[158]
 x[159] x[160] x[161] x[162] x[163] x[164] x[165] x[166] x[167] x[168]
 x[169] x[170] x[171] x[172] x[173] x[174] x[175] x[176] x[177] x[178]
 x[179] x[180] x[181] x[182] x[183] x[184] x[185] x[186] x[187] x[188]
 x[189] x[190] x[191] x[192] x[193] x[194] x[195] x[196] x[197] x[198]
 x[199] x[200] x[201] x[202] x[203] x[204] x[205] x[206] x[207] x[208]
 x[209] x[210] x[211] x[212] x[213] x[214] x[215] x[216] x[217] x[218]
 x[219] x[220] x[221] x[222] x[223] x[224] x[225] x[226] x[227] x[228]
 x[229] x[230] x[231] x[232] x[233] x[234] x[235] x[236] x[237] x[238]
 x[239] x[240] x[241] x[242] x[243] x[244] x[245] x[246] x[247] x[248]
 x[249] x[250] x[251] x[252] x[253] x[254] x[255] x[256] x[257] x[258]
 x[259] x[260] x[261] x[262] x[263] x[264] x[265] x[266] x[267] x[268]
 x[269] x[270] x[271] x[272] x[273] x[274] x[275] x[276] x[277] x[278]
 x[279] x[280] x[281] x[282] x[283] x[284] x[285] x[286] x[287] x[288]
 x[289] x[290] x[291] x[292] x[293] x[294] x[295] x[296] x[297] x[298]
 x[299] x[300] x[301] x[302] x[303] x[304] x[305] x[306] x[307] x[308]
 x[309] x[310] x[311] x[312] x[313] x[314] x[315] x[316] x[317] x[318]
 x[319] x[320] x[321] x[322] x[323] x[324] x[325] x[326] x[327] x[328]
 x[329] x[330] x[331] x[332] x[333] x[334] x[335] x[336] x[337] x[338]
 x[339] x[340] x[341] x[342] x[343] x[344] x[345] x[346] x[347] x[348]
 x[349] x[350] x[351] x[352] x[353] x[354] x[355] x[356] x[357] x[358]
 x[359] x[360] x[361] x[362] x[363] x[364] x[365] x[366] x[367] x[368]
 x[369] x[370] x[371] x[372] x[373] x[374] x[375] x[376] x[377] x[378]
 x[379] x[380] x[381] x[382] x[383] x[384] x[385] x[386] x[387] x[388]
 x[389] x[390] x[391] x[392] x[393] x[394] x[395] x[396] x[397] x[398]
 x[399] x[400] x[401] x[402] x[403] x[404] x[405] x[406] x[407] x[408]
 x[409] x[410] x[411] x[412] x[413] x[414] x[415] x[416] x[417] x[418]
 x[419] x[420] x[421] x[422] x[423] x[424] x[425] x[426] x[427] x[428]
 x[429] x[430] x[431] x[432] x[433] x[434] x[435] x[436] x[437] x[438]
 x[439] x[440] x[441] x[442] x[443] x[444] x[445] x[446] x[447] x[448]
 x[449] x[450] x[451] x[452] x[453] x[454] x[455] x[456] x[457] x[458]
 x[459] x[460] x[461] x[462] x[463] x[464] x[465] x[466] x[467] x[468]
 x[469] x[470] x[471] x[472] x[473] x[474] x[475] x[476] x[477] x[478]
 x[479] x[480] x[481] x[482] x[483] x[484] x[485] x[486] x[487] x[488]
 x[489] x[490] x[491] x[492] x[493] x[494] x[495] x[496] x[497] x[498]
 x[499] x[500] x[501] x[502] x[503] x[504] x[505] x[506] x[507] x[508]
 x[509] x[510] x[511] sout[0] sout[1] sout[2] sout[3] sout[4] sout[5]
 sout[6] sout[7] sout[8] sout[9] sout[10] sout[11] sout[12] sout[13]
 sout[14] sout[15] sout[16] sout[17] sout[18] sout[19] sout[20] sout[21]
 sout[22] sout[23] sout[24] sout[25] sout[26] sout[27] sout[28] sout[29]
 sout[30] sout[31] sout[32] sout[33] sout[34] sout[35] sout[36] sout[37]
 sout[38] sout[39] sout[40] sout[41] sout[42] sout[43] sout[44] sout[45]
 sout[46] sout[47] sout[48] sout[49] sout[50] sout[51] sout[52] sout[53]
 sout[54] sout[55] sout[56] sout[57] sout[58] sout[59] sout[60] sout[61]
 sout[62] sout[63] sout[64] sout[65] sout[66] sout[67] sout[68] sout[69]
 sout[70] sout[71] sout[72] sout[73] sout[74] sout[75] sout[76] sout[77]
 sout[78] sout[79] sout[80] sout[81] sout[82] sout[83] sout[84] sout[85]
 sout[86] sout[87] sout[88] sout[89] sout[90] sout[91] sout[92] sout[93]
 sout[94] sout[95] sout[96] sout[97] sout[98] sout[99] sout[100] sout[101]
 sout[102] sout[103] sout[104] sout[105] sout[106] sout[107] sout[108]
 sout[109] sout[110] sout[111] sout[112] sout[113] sout[114] sout[115]
 sout[116] sout[117] sout[118] sout[119] sout[120] sout[121] sout[122]
 sout[123] sout[124] sout[125] sout[126] sout[127] sout[128] sout[129]
 sout[130] sout[131] sout[132] sout[133] sout[134] sout[135] sout[136]
 sout[137] sout[138] sout[139] sout[140] sout[141] sout[142] sout[143]
 sout[144] sout[145] sout[146] sout[147] sout[148] sout[149] sout[150]
 sout[151] sout[152] sout[153] sout[154] sout[155] sout[156] sout[157]
 sout[158] sout[159] sout[160] sout[161] sout[162] sout[163] sout[164]
 sout[165] sout[166] sout[167] sout[168] sout[169] sout[170] sout[171]
 sout[172] sout[173] sout[174] sout[175] sout[176] sout[177] sout[178]
 sout[179] sout[180] sout[181] sout[182] sout[183] sout[184] sout[185]
 sout[186] sout[187] sout[188] sout[189] sout[190] sout[191] xorout[0]
 xorout[1] xorout[2] xorout[3] xorout[4] xorout[5] xorout[6] xorout[7]
 xorout[8] xorout[9] xorout[10] xorout[11] xorout[12] xorout[13] xorout[14]
 xorout[15] xorout[16] xorout[17] xorout[18] xorout[19] xorout[20]
 xorout[21] xorout[22] xorout[23] xorout[24] xorout[25] xorout[26]
 xorout[27] xorout[28] xorout[29] xorout[30] xorout[31] xorout[32]
 xorout[33] xorout[34] xorout[35] xorout[36] xorout[37] xorout[38]
 xorout[39] xorout[40] xorout[41] xorout[42] xorout[43] xorout[44]
 xorout[45] xorout[46] xorout[47] xorout[48] xorout[49] xorout[50]
 xorout[51] xorout[52] xorout[53] xorout[54] xorout[55] xorout[56]
 xorout[57] xorout[58] xorout[59] xorout[60] xorout[61] xorout[62]
 xorout[63] xorout[64] xorout[65] xorout[66] xorout[67] xorout[68]
 xorout[69] xorout[70] xorout[71] xorout[72] xorout[73] xorout[74]
 xorout[75] xorout[76] xorout[77] xorout[78] xorout[79] xorout[80]
 xorout[81] xorout[82] xorout[83] xorout[84] xorout[85] xorout[86]
 xorout[87] xorout[88] xorout[89] xorout[90] xorout[91] xorout[92]
 xorout[93] xorout[94] xorout[95] xorout[96] xorout[97] xorout[98]
 xorout[99] xorout[100] xorout[101] xorout[102] xorout[103] xorout[104]
 xorout[105] xorout[106] xorout[107] xorout[108] xorout[109] xorout[110]
 xorout[111] xorout[112] xorout[113] xorout[114] xorout[115] xorout[116]
 xorout[117] xorout[118] xorout[119] xorout[120] xorout[121] xorout[122]
 xorout[123] xorout[124] xorout[125] xorout[126] xorout[127] xorout[128]
 xorout[129] xorout[130] xorout[131] xorout[132] xorout[133] xorout[134]
 xorout[135] xorout[136] xorout[137] xorout[138] xorout[139] xorout[140]
 xorout[141] xorout[142] xorout[143] xorout[144] xorout[145] xorout[146]
 xorout[147] xorout[148] xorout[149] xorout[150] xorout[151] xorout[152]
 xorout[153] xorout[154] xorout[155] xorout[156] xorout[157] xorout[158]
 xorout[159] xorout[160] xorout[161] xorout[162] xorout[163] xorout[164]
 xorout[165] xorout[166] xorout[167] xorout[168] xorout[169] xorout[170]
 xorout[171] xorout[172] xorout[173] xorout[174] xorout[175] xorout[176]
 xorout[177] xorout[178] xorout[179] xorout[180] xorout[181] xorout[182]
 xorout[183] xorout[184] xorout[185] xorout[186] xorout[187] xorout[188]
 xorout[189] xorout[190] xorout[191]
End
