FfwdRead = range(50)
FrvsRead = range(50)
F_list = range(50)
outputWins = range(50)
# test output
outBuff = '=' * 25 + '\n'
#print('=' * 25)
outBuff += 'Free Energy Data:\n'
#print('Free Energy Data:')
outBuff += '%-8s %-8s %-8s %-8s\n' % ('Lambda','F_BAR', 'F_fwd', 'F_rvs' )
#print('%-8s %-8s %-8s %-8s' % ('Lambda','F_BAR', 'F_fwd', 'F_rvs' ))
for i in range(len(F_list)):
    outBuff += '%4.2f  %8.4f  %8.4f %8.4f\n' % (outputWins[i], F_list[i], FfwdRead[i], FrvsRead[i])
    #print('%4.2f  %8.4f  %8.4f %8.4f' % (outputWins[i], F_list[i], FfwdRead[i], FrvsRead[i]))

open('ParsedFep.dat', 'w').write(outBuff)
print(outBuff)