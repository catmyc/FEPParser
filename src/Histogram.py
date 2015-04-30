__Author__ = "Yuncheng Mao"
__Email__ = '''
                catmyc@gmail.com
                maoyuncheng@mail.nankai.edu.cn
            '''
            
class Histogram:
    '''
    data, stores raw data;
    A member hist, as a dictionary, is offered to store the calculated histogram.
    '''
    def __init__(self, DataArray): # Read in data
        '''
        Read data from a list or dictionary
        '''
        print('-' * 40)
        if isinstance(DataArray, list):
            self.data = DataArray # reference, not copy!
        if isinstance(DataArray, dict):
            self.data = [v for (k,v) in DataArray.items()]
            del k
        print('Raw data for histogram analysis read.')
        print('Sorting data...')
        self.data.sort()
        '''
        Print test:
        '''
        #print(self.data)
        self.size = len(self.data)
        self.min = self.data[0]
        self.max = self.data[-1]
        print('Input summary: Max value is %.4f, Min value is %.4f, number of data points is %d.' % (self.max, self.min, self.size))
    
    def stat(self, Nintervals = 200):
        '''
        Calculate histogram.
        By default, 100 statistical intervals will be applied.
        '''
        incr = (self.max - self.min) / Nintervals # increment
        current = self.min + incr
        count = [0] * Nintervals # hits in each interval
        index = 0
        print('Scanning through sorted data list...')
        for f in self.data: # scan through the list
            if f <= current:
                count[index] += 1
            elif f < self.max:
                index += 1
                count[index] += 1
                current += incr
            else: # to prevent the out-of-range error
                count[Nintervals -1] += 1
        x = [self.min + incr * (0.5 + i) for i in range(Nintervals)]
        P = [c / self.size for c in count]
        self.hist = dict(zip(x,P))
        print('Statistics completed.')
        
    def printhist(self, filehandle=None):
        outBuff = '# Value       Probability\n'
        from operator import itemgetter
        for (V, P) in sorted(self.hist.items(), key=itemgetter(0)):
            outBuff += '%-12.4f  %-8.4f\n' % (V, P)
        outBuff += '#' + '-' * 40 + '\n'
        if filehandle == None: # print to stdout
            print(outBuff)
        else:
            filehandle.write(outBuff) 