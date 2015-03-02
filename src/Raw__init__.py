from BAR import *
from fepInterpretor import *
fwdfilename = r'D:\work\2012_Carrier\FEP_CD_Anihilation\Alchem_C1_Unbound\forward.fepout'
rvsfilename = r'D:\work\2012_Carrier\FEP_CD_Anihilation\Alchem_C1_Unbound\backward.fepout'
# Read and compile fep data
Fwd = fepData(fwdfilename)
Rvs = fepData(rvsfilename)
#----------------------------
# BAR estimation
Wins = Fwd.WinList.copy()
outputWins = [0] + Wins
F_list = [0] # list for Free energy
DeltaF_list = [] # list for free-energy changes

FfwdRead = [0]
FrvsRead = [0]
'''
Test
-----------------------------------------------------------------
Win=Fwd.WinList[0]
bar = BARestimator(Fwd.DeltaU_All[Win],Rvs.DeltaU_All[Win], 300)
F = bar.BARSC(Fwd.F_read[Win])
print(F)
------------------------------------------------------------------
Test OK!!!
'''
temp = 300
print('Temperature is set to 300K!')
F_now = 0 #accumulator for free energy change
for win in Wins:
    print('#' * 40)
    print('Calculating free energy change of win [ %4.2f %4.2f ]:' % (win-0.02, win) )
    bar = BARestimator(Fwd.DeltaU_All[win], Rvs.DeltaU_All[win], temp)
    DF = bar.BARSC(Fwd.F_read[win], 1e-9)
    DeltaF_list.append(DF)
    F_now += DF
    F_list.append(F_now)
    FfwdRead.append(Fwd.F_read[win])
    FrvsRead.append(Rvs.F_read[win])


# test output
outBuff = '=' * 25 + '\n'
outBuff += 'Free Energy Data:\n'
outBuff += '%-8s %-8s %-8s %-8s\n' % ('Lambda','F_BAR', 'F_fwd', 'F_rvs' )
for i in range(len(F_list)):
    outBuff += '%4.2f  %8.4f  %8.4f %8.4f\n' % (outputWins[i], F_list[i], FfwdRead[i], FrvsRead[i])
    
open('ParsedFep.dat', 'w').write(outBuff)
print(outBuff)