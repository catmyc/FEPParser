from fepInterpretor import *
name = r'D:\work\2012_Carrier\FEP_CD_Anihilation\Alchem_C1_Unbound\forward.fepout'
f = fepData(name, False)
f.genHist(50)
f.printHist(False)
f.printHist(True)
