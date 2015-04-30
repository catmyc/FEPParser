__Author__ = "Yuncheng Mao"
__Email__ = '''
            catmyc@gmail.com
            maoyuncheng@mail.nankai.edu.cn
            '''
'''
To-do list:
-------------------------------------
1. User define the maximun iterations numbers for bar estimation.
2. A GUI.
3. Plotting functions for free energy and histogram.
4. Using bootstrap to estimate the error of BAR calculation
-------------------------------------
'''
from sys import exit
from vector import *
import BAR

# forward and backward mark
FWD = 'fwd'
BWD = 'bwd'

def genHist(Datalist, min_value, max_value, Nintervals=250):
    '''
    Generate a unified histogram for given data.
    ======================
    Datalist: an iterable that contains all the data
    '''
    total = 0
    delta = (max_value - min_value) / Nintervals
    hist = [0] * Nintervals
    for f in Datalist:
        if f < min or f > max:
            # discard data out of range
            continue 
        total += 1
        ind = int((f - min) / delta)
        hist[ind] += 1
    return [(f / total) for f in hist]

class fepWin:
    '''
    For storing the data of a fep window. Including the attributes:
    lambda1
    lambda2
    direction: FWD or BWD
    NSamples: number of sampling
    W_list: a list for all recorded work
    var: the variance of works recorded, useful for error estimation
    label: labeling the window using lambda values, e.g. a string '0.01-0.02'
    label_format: by default %.2f-%.2f 
    max_num_win: by default 100 wins.
    meanW: average work
    '''
    def __init__(self, max_nWin = 100):
        format_of_label_format = '%.{0}f-%.{0}f'
        # determine the label format first
        i = 1
        while True:
            if 10**i >= max_nWin:
                self.max_num_win = 10**i
                self.label_format = format_of_label_format.format(i)
                break
            else:
                i += 1
        self.lambda1 = 0.0
        self.lambda2 = 0.0
        self.direction = None
        self.NSamples = 0
        self.W_list = []
        self.var = 0.0
        self.label = None 
        self.meanW = 0.0
        self.temperature = 303.15 # default temperature, required in BAR estimation.
    
    def set_temperature(self, temp):
        '''
        Method called to set the temperature of the fep window.
        '''
        self.temperature = temp
    
    def copy(self):
        '''
        Make a copy if sometimes necessary
        '''
        newWin = fepWin()
        newWin.lambda1 = self.lambda1
        newWin.lambda2 = self.lambda2
        newWin.direction = self.direction
        newWin.label = self.label
        newWin.label_format = self.label_format
        newWin.max_num_win = self.max_num_win
        newWin.NSamples = self.NSamples
        newWin.W_list = self.W_list.copy()
        newWin.var = self.var
        newWin.meanW = self.meanW
    
    def set(self, l1, l2, works = []):
        self.lambda1 = l1
        self.lambda2 = l2
        if self.lambda1 < self.lambda2:
            self.direction = FWD
        elif self.lambda1 > self.lambda2:
            self.direction = BWD
        else:
            print('ERROR: Lambda1 is equal to Lambda2! Exit. Check you input. Lambda1 = Lambda2 = %.2f' % self.lambda1)
            exit()
        if len(works) == 0:
            print('ERROR: No recorded work given. Exit.')
            exit()
        self.NSamples = len(works)
        self.W_list = works
        self.var = vecvar(self.W_list)
        self.label = self.label_format % (min(l1, l2), max(l1, l2))
        self.meanW = vecmean(self.W_list)
        print('FEP window recorded. Lable: %s.' % self.label, self.direction)
    
    def set_from_file(self, filename):
        '''
        If the data of a single is extracted from the fepout file, then the user can read from such file the data of the FEP window.
        '''
        print('Reading from file. Note that the file must not contain more than one windows!')
        W = []
        for line in open(filename):
            buff = line.split()
            if line.startswith('#NEW'): 
                l1 = float(buff[6])
                l2 = float(buff[8])
            if line.startswith('FepEnergy:'):
                W.append(float(buff[6]))
        self.set(l1, l2, W)
    
    def clear_raw_data(self):
        '''
        Delete the work list to save memory if needed. E.g. such data is no longer needed after all calculations finished for a window pair.
        '''
        del self.W_list
        
        
    @classmethod
    def winYield(cls, filename, Temp, maxWin=100):
        '''
        Go through the fepout file and generate a series of fepWin objects.
        This is a fepout file reader.
        =====================
        Received parameters:
        filename: the fep output file
        Temp: Temperature assigned for BAR estimation, unit in K (Kelvin).
        '''
        l1 = 0
        l2 = 0
        W = []
        for line in open(filename):
            buff = line.split()
            if line.startswith('#Free energy change'): # End of a window
                fwin = fepWin(maxWin)
                fwin.set(l1, l2, W)
                fwin.set_temperature(Temp)
                # reset
                l1 = 0
                l2 = 0
                W = []
                yield fwin
            if line.startswith('#NEW'): # Beginning of a new window
                l1 = float(buff[6])
                l2 = float(buff[8])
            if line.startswith('FepEnergy:'):
                W.append(float(buff[6]))
    
class fepHistogram:
    def __init__(self, l, hist_f, hist_b):
        '''
        l: x coordinates
        hist_f: unified histogram of forward works
        hist_b: unified histogram of backward works
        '''
        self.fwd_hist = hist_f
        self.bwd_hist = hist_b

class winPair:
    '''
    Constructed by the corresponding forward and backward windows, and estimate the free energy change using BAR.
    The label of the forward and backward windows must match!
    -------------------------------------
    fwd_win: forward window
    bwd_win: backward window
    DF: DeltaF, the estimated free-energy change from BAR
    error_F: estimated error of free energy.
    
    calcDF(): calculate the free energy difference
    calcError(): estimate the error of BAR estimation
    clear_raw_data(): delete the work lists of the forward and backward windows.
    '''
    def __init__(self, win_f, win_b):
        if win_f.label != win_b.label:
            print('ERROR: Must be within the same windows! Exit.')
            exit()
        if win_f.direction == win_b.direction:
            print('ERROR: Must be opposite directions! Exit.')
            exit()
        if win_f.temperature != win_b.temperature:
            print('Windows are NOT assigned with the same temperature! Exit.')
            exit()
        self.fwd_win = win_f
        self.bwd_win = win_b 
        self.label = win_f.label
        self.hist_f = None
        self.hist_b = None
        self.DF = 0
        self.error_F = 0
    
    def calcDF(self):
        '''
        This is the core function of fep analysis.
        '''
        barMachine = BAR.BARestimator(self.fwd_win.W_list, self.bwd_win.W_list)
        self.DF = barMachine.BARSC(self.fwd_win.meanW) # The default parameters are good enough
        print("Free energy change for window [ %s ] calculated." % self.label)
    
    def calcError(self):
        '''
        Use bootstrap to estimate the error of BAR estimation.
        '''
        pass
    
    def calcHist(self, Nintervals = 250):
        '''
        Generate the histogram based on both work lists
        '''
        min_W = min(min(self.fwd_win.W_list), min(self.bwd_win.W_list))
        max_W = max(max(self.fwd_win.W_list), max(self.bwd_win.W_list))
        self.hist_f = genHist(self.fwd_win.W_list, min_W, max_W, Nintervals)
        self.hist_b = genHist(self.bwd_win.W_list, min_W, max_W, Nintervals)
        
    
    def clear_raw_data(self):
        self.fwd_win.clear_raw_data()
        self.bwd_win.clear_raw_data()



def run():
    max_num_wins = 100
    temperature = 298
    fwd_filename = r'D:\work\2012_Carrier\FEP_CD_Anihilation\Alchem_C1_Unbound\forward.fepout'
    bwd_filename = r'D:\work\2012_Carrier\FEP_CD_Anihilation\Alchem_C1_Unbound\backward.fepout'
    print("Reading forward fepout file...")
    fwd_wins = [win for win in fepWin.winYield(fwd_filename, temperature, max_num_wins)]
    print("Reading backward fepout file...")
    bwd_wins = [win for win in fepWin.winYield(bwd_filename, temperature, max_num_wins)]
    # Debug
    fw = fwd_wins[0]
    for w in bwd_wins:
        if w.label == fw.label:
            print("Found!!!!!!!")
    
    print(w.label)
    
    return
    

if __name__ == "__main__":
    run()

"""
class fepData:
    '''
    Read in fep data from fepout file.
    Fep data is compiled spontaneously while creating the object.
    
    '''
    def __init__(self, Filename, ifSplit = False):
        self.associatedFilename = Filename
        self.doSplit = ifSplit
        self.L1 = []
        self.L2 = []
        self.WinList = []
        self.histograms = dict() # map each win to a histogram
        self.NumWins = 0
        self.fororback = ''
        self.F_read = dict() # The Free energy change read from file.
        if self.doSplit:
            self.doDirectRead = False
            print('Processing data by splitting fepout file.')
            self.tmpFileList = []
            self.splitfepout()
        else: # Default using direct reading.
            self.doDirectRead = True
            self.DeltaU_All = dict() #Use dictionary to store all DeltaU data, using 0.02, 0.04, ..., 1.00 as the key.
            print('Processing data by direct read. Enough memory is required on your machine!!!')
            self.directRead()
        
    def splitfepout(self): # if memory is not enough, trade space for efficiency
        '''
        Split the .fepout file into small files, each containing data for one window.
        Only the data of work in a window is recorded. The vdW and Elec terms are discarded.
        filename: name of the .fepout file
        '''
        infile = open(self.associatedFilename, 'r')
        for line in infile:
            buff = line.split()
            if line.startswith('#NEW'):
                print('-' * 20)
                self.NumWins += 1
                l1 = float(buff[6])
                l2 = float(buff[8])
                self.L1.append(l1)
                self.L2.append(l2)
                if l1 < l2:
                    self.fororback = 'Fwd'
                else:
                    self.fororback = 'Rvs'
                tmpFileName = self.fororback + '_%.2f_%.2f.dat' % (l1, l2)
                self.tmpFileList.append(tmpFileName)
                print('Writing tmp file: %s' % tmpFileName)
                outfile = open(tmpFileName, 'w')
                print('#%4s, Lambda1=%5.2f, Lambda2=%5.2f' %(self.fororback, l1, l2))
                outfile.writelines('#%4s, Lambda1=%5.2f, Lambda2=%5.2f' %(self.fororback, l1, l2))
                DeltaU = [] # initialize the DeltaU list
            if line.startswith('FepEnergy:'):
                DeltaU.append(float(buff[6]))
            if line.startswith('#Free'):
                print('Length of DeltaU list is: %d.' % len(DeltaU))
                initialGuess = float(buff[11])
                if self.fororback == 'Fwd':
                    self.F_read[l1] = initialGuess
                else:
                    self.F_read[l2] = initialGuess
                # print('initialGuess is', initialGuess)
                print('Free energy change recorded in file is: %12.4f' % initialGuess)
                outfile.write('InitialGuess: %12.4f\n' % initialGuess)
                outfile.write('DeltaU: ' + ('%12.4f ' * len(DeltaU) % tuple(DeltaU)) + '\n')
                outfile.close()
                print('Finished reading window [ %4.2f, %4.2f]' % (l1, l2))
                print('Totally %d data entries read!!!' % len(DeltaU))
                print('Temporary file %s written.' % tmpFileName)
    
    def directRead(self):
        '''
        Directly read and compile data from .fepout files.
        All the data is stored in RAM. Enough memory is required for this method.
        This method is expected to have better performance.
        '''
        infile = open(self.associatedFilename, 'r')
        for line in infile:
            buff = line.split()
            if line.startswith('#NEW'):
                print('-' * 20)
                DeltaU = []
                self.NumWins += 1
                l1 = float(buff[6])
                l2 = float(buff[8])
                self.L1.append(l1)
                self.L2.append(l2)
                from math import fabs
                self.WinWidth = fabs(l1 - l2)
                if l1 < l2:
                    self.fororback = 'Fwd'
                else:
                    self.fororback = 'Rvs'
                print('#%4s, Lambda1=%5.2f, Lambda2=%5.2f' %(self.fororback, l1, l2))
            if line.startswith('FepEnergy:'):
                DeltaU.append(float(buff[6]))
            if line.startswith('#Free'):
                print('Length of DeltaU list is: %d.' % len(DeltaU))
                initialGuess = float(buff[11])
                if self.fororback == 'Fwd':
                    self.F_read[l2] = initialGuess
                else:
                    self.F_read[l1] = initialGuess
                if self.fororback == 'Fwd':
                    self.DeltaU_All[l2] = DeltaU
                else:
                    self.DeltaU_All[l1] = DeltaU
                print('Finished reading window [ %4.2f, %4.2f]' % (l1, l2))
                print('Free energy change recorded in file is: %12.4f' % initialGuess)
                print('Totally %d data entries read!!!' % len(DeltaU))
        if self.fororback == 'Fwd':
            self.WinList = self.L2
        else:
            self.WinList = self.L1
            
    def genHist(self, Nintervals):
        '''
        Generate the histogram for each window
        '''
        print('=' * 40)
        print('Generating histograms for each window...')
        for win in self.WinList:
            print('-' * 40)
            print('Processing win: [ %4.2f, %4.2f ]' % (win-self.WinWidth, win))
            self.histograms[win] = Histogram(self.DeltaU_All[win])
            self.histograms[win].stat(Nintervals)
        print('Histogram generated!')
        
    def printHist(self, ifWriteFile = True):
        '''
        Print histograms to stdout or file.
        Default write to file.
        '''
        if len(self.histograms) == 0:
            print('Histograms do not exist, need to be generated first!!!')
            return
        print('-' * 40)
        print('Output histogram data')
        if ifWriteFile:
            filename = self.fororback + '_Histogram.dat'
            print('Writing to file %s...' % filename)
            file = open(filename, 'w')
            for win in self.WinList:
                print('-' * 40)
                print('Writing Probability of Window [ %4.2f, %4.2f ]' % (win - self.WinWidth, win))
                file.write('# Probability of Window [ %4.2f, %4.2f ]\n' % (win - self.WinWidth, win))
                self.histograms[win].printhist(file)
            file.close()
            print('File %s written.' % filename)
        else:
            print('Writing histogram to stdout...')
            for win in self.WinList:
                print('# Probability of Window [ %4.2f, %4.2f ]' % (win - self.WinWidth, win))
                self.histograms[win].printhist()

"""                