__Author__ = "Yuncheng Mao"
__Email__ = '''
                catmyc@gmail.com
                maoyuncheng@mail.nankai.edu.cn
            '''
from Histogram import * 
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
        Split .fepout file into small files, each containing data for one window.
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