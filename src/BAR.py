__Author__ = "Yuncheng Mao"
__Email__ = '''
                catmyc@gmail.com
                maoyuncheng@mail.nankai.edu.cn
            '''

import math
def _logmean(exp_args):
    '''
    exp_args: a list of values
    _logmean(exp_arg) = log(mean(1/(1+exp(exp_args))))
    '''
    '''
    # Old version of codes
    S = 0
    for arg in exp_args:
        S += 1 / (1 + math.exp(arg))
    return math.log(S / len(exp_args))
    '''
    S = sum([1/(1 + math.exp(f)) for f in exp_args])
    return math.log(S / len(exp_args))
    
class BARestimator:
    '''
    The BARestimator class for BAR analysis.
    Using W_F and W_R to record the forward and backward work.
    This is for estimating the free-energy change of one window!
    ''' 
    def __init__(self, Fwd, Rvs, Temperature):
        '''
        parameters:
        -----------
        Fwd: Forward work
        Rvs: Reverse work
        Temperature: As the name says.
        '''
        self.W_F = Fwd
        self.W_R = Rvs
        self.Temp = Temperature
        if len(self.W_F) == 0 or len(self.W_R) == 0:
            print("This is an EMPTY BARestimator!!!")
            self.isempty = True
        else:
            info = "Read in forward work list of length {0}, reverse work list of length {1}".format(len(self.W_F), len(self.W_R))
            print(info)
            self.isempty = False
    
    def readFromFile(self, FwdName, RvsName, Temperature):
        '''
        Read the work list from file.
        '''
        self.Temp = Temperature
        pass # Come to it later
    
    def BARzero(self, DeltaF = 0):
        '''
        Single BAR-estimation run.
        The return value of this function is used to determine convergence.
        Parameter DeltaF is an initial guess or input value. 
        The DeltaF value from FEP output should be taken as the initial guess. If not specified, 0 will be taken as input. 
        '''
        k = 0.001987200 # Boltzmann constant
        kT = k * self.Temp #In unit kcal/mol
        exp_args_F = []
        exp_args_R = []
        for wf in self.W_F:
            exp_args_F.append(wf - DeltaF)
        for wr in self.W_R:
            exp_args_R.append(wr + DeltaF)
        logF = _logmean(exp_args_F)
        logR = _logmean(exp_args_R)
        '''
        Test use:
        print("logF = %8.4f" % logF)
        print("logR = %8.4f" % logR)
        '''
        return kT * (logR - logF)
        
    
    def BARSC(self, DeltaF = 0, convergence = 1e-8, MAXITER = 1000): #self-consistent estimation
        '''
        BARzero < convergence is satisfied when achieving convergence.
        Maximum number of BAR estimations is MAXITER, if convergence is not reached. 
        Input parameter DeltaF is an initial guess.
        '''
        for iteration in range(MAXITER):
            token = self.BARzero(DeltaF)
            if math.fabs(token) < convergence:
                break
            else:
                DeltaF = DeltaF + token
        if iteration < MAXITER:
            print("Convergence achieved, after %d iterations!" % (iteration + 1))
        else:
            print("Maximum number of iteration reached!")
        print('Estimated Free energy change is: %7.4f' % DeltaF)
        return DeltaF
        