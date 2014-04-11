# CONSTANTS:
SAMPLES_PER_ROW = 
POS_CONTROLS_PER_ROW =
NEG_CONTROLS_PER_ROW =  

import sys
import math

def read_spectra_strings(filename):
    """Reads time layout text output file from SpectraMax"""
    with open(filename) as file:
        plates = file.read().split('~End') #separates each plate
    plates.pop() #removes extra '\n' character
    for i in range(len(plates)):
        plates[i] = plates[i].split('Temperature(\xa1C)') #separates 'labeling info' from 'well names' and 'values'
        plates[i] = plates[i][1].split('\t\n\t') #separates 'well names' and 'values'
        for j in range(len(plates[i])):
            plates[i][j] = plates[i][j].split('\t') #separates each 'well name' and 'value'
            plates[i][0][0] = 'Temperature' #Adds temperature label that was previously removed
        plates[i] = dict(zip(plates[i][0], plates[i][1])) #('well name', 'value')
    return plates

def load_assay_layout(filename):
    """Loads assay plate layout"""
    with open(filename) as file:
        layout = file.read().replace('\n',',')
    layout = layout.split(',,,,,,,,,,,,,,,,,,,,,,,,,') #separates well descriptor from well location
    for i in range(len(layout)):
        layout[i] = layout[i].split(',') #separates each well by comma
    layout[1].pop() # removes extra '\n' character from end        
    layout[0] = dict(zip(layout[1],layout[0])) 
    layout.pop()
    return layout

def read_info(filename, plates):
    """Loads concentration/compound information"""
    with open(filename) as file:
        file.readline()
        info = file.readlines()
        x = [{} for i in range(plates)]
        for i in range(len(info)):
            temp = info[i].split(',')
            j = int(temp[0]) - 1
            x[j][temp[1]] = [temp[2], temp[3]] #Well, Concentration, Sample ID
        info = x
    return info

def mean(alist, lognorm = False):
    if lognorm == False:
        m = math.fsum(alist) / len(alist)
    elif lognorm == True:
        xlist = list()
        for i in range(len(alist)):
            xlist.append(math.log10(alist[i]))
        m = mean(xlist)
    return m

def stdv(alist, lognorm = False): 
    sd = 0
    m = mean(alist, lognorm)
    if lognorm == False:
        for i in range(len(alist)):
            sd += math.pow(alist[i] - m, 2)
        sd = sd / (len(alist) - 1)
        sd = math.sqrt(sd)
    if lognorm == True:
        xlist = list()
        for i in range(len(alist)):
            xlist.append(math.log10(alist[i]))
        sd = stdv(xlist)
    return sd

def zscore(pmean, pstdev, nmean, nstdev):
    ans = 1 - (3*(pstdev + nstdev)/math.fabs(pmean - nmean))
    return ans

def median(alist):
    """Finds the median of a list"""
    alist.sort()
    center = len(alist) / 2.0
    if len(alist) % 2 == 0:
        median = (alist[int(center) - 1] + alist[int(center)])/2.0
    else:
        median = float(alist[int(center)])
    return {'median':median, 'index':center}

def quartile(alist):
    """Finds 1st and 3rd quartiles"""
    lower = alist[0:int(math.floor(median(alist)['index']))]
    upper = alist[int(math.ceil(median(alist)['index'])):]
    return {'Q1':median(lower)['median'], 'Q3':median(upper)['median']}

def bound(alist, Range = 1.5):
    quartiles = quartile(alist)
    IQR = abs(quartiles['Q3'] - quartiles['Q1'])
    lowerbound = quartiles['Q1'] - Range*IQR
    upperbound = quartiles['Q3'] + Range*IQR
    return {'lower':lowerbound, 'upper':upperbound}

def outlier(alist, Range = 1.5):
    bounds = bound(alist, Range)
    outliers = list()
    for i in range(len(alist)):
        if alist[i] < bounds['lower'] or alist[i] > bounds['upper']:
            outliers.append(alist[i])
    return outliers

def normalize(sample, positive, negative, lognorm = False):
    if lognorm == False:
        normalized = (sample - positive) / (negative - positive) * 100
    if lognorm == True:
        normalized = (math.log10(sample) - positive) / (negative - positive) * 100        
    return normalized

def in_list(alist):
    values = list()
    for i in range(len(alist)):
        if alist[i] not in values:
            values.append(alist[i])
    return sorted(values)

argc = len(sys.argv)

if argc is 1:
    print 'Command line error: No input file'
    sys.exit(1)
elif argc is 2: #Only SpectraMax input (uses standard plate layout)
    layout = load_assay_layout("default_assay_template.csv")
    info = read_info('drug_info.csv', 23)
elif argc is 3: #SpectraMax input and custom plate layout
    layout = load_assay_layout("%s" % sys.argv[3])
elif argc > 3:
    print 'Command line error: Wrong Number of Arguments'
    sys.exit(1)



#Set-up options for normalization
print "\nSelect Y/N for the following options:"
if raw_input("Use Default Options? ") == 'Y':
    graphpad_format_use = False
    row_normal_use = False
    lognorm_use = False
    median_use = False
    outlier_use = False
    info_use = False
else:
    if raw_input("Format for GraphPad? ") == 'Y':
        graphpad_format_use = True
    else:
        graphpad_format_use = False
    if raw_input("Use Row Normalization? ") == 'Y':
        row_normal_use = True
    else:
        row_normal_use = False

    if raw_input("Log-Normalize? ") == 'Y':
        lognorm_use = True
    else:
        lognorm_use = False

    if raw_input("Use Median (instead of mean)? ") == 'Y':
        median_use = True
    else:
        median_use = False

    if raw_input("Remove outliers? ") == 'Y':
        outlier_use = True
        outlier_threshold = 1.5
        outlier_threshold = float(raw_input("Enter outlier threshold (default = 1.5): "))
    else:
        outlier_use = False 

#Constants
dilutions_num = 12

#Open file and separate postive/negative controls and samples
filename = sys.argv[1]
plates = read_spectra_strings("%s" % filename)
plates_len = len(plates)

samples =  [{} for i in range(plates_len)]
control_pos =  [{} for i in range(plates_len)]
control_neg =  [{} for i in range(plates_len)]
samples_keys = [list() for i in range(plates_len)] #Ordered list of keys for samples

for i in range(plates_len):
    for j in range(ord('A'), ord('Q')):
        for k in range(1, 25):
            RFU = float(plates[i]['%c%d' % (j, k)])    
            if layout[0]['%c%d' % (j, k)] == 'S':
                samples[i]['%c%d' % (j, k)] = RFU
                samples_keys[i].append('%c%d' % (j, k))
            elif layout[0]['%c%d' % (j, k)] == 'PC':
                control_pos[i]['%c%d' % (j, k)] = RFU
            elif layout[0]['%c%d' % (j, k)] == 'NC':
                control_neg[i]['%c%d' % (j, k)] = RFU


#Removal of Outliers (If 'Y' in Options)
if outlier_use:
    control_pos_all = list()
    control_neg_all = list()
    for i in range(plates_len):
        control_pos_all += control_pos[i].values()
    for i in range(plates_len):
        control_neg_all += control_neg[i].values()

    outlier_control_pos = outlier(control_pos_all)
    outlier_control_neg = outlier(control_neg_all)

    for i in range(plates_len):
        pos_temp = dict()
        neg_temp = dict()
        for j in range(len(control_pos[i].keys())):
            if control_pos[i][control_pos[i].keys()[j]] not in outlier_control_pos:
                pos_temp[control_pos[i].keys()[j]] = control_pos[i][control_pos[i].keys()[j]]
            else:
                print 'deleted positive control'
        for j in range(len(control_neg[i].keys())):
            if control_neg[i][control_neg[i].keys()[j]] not in outlier_control_neg:
                neg_temp[control_neg[i].keys()[j]] = control_neg[i][control_neg[i].keys()[j]]
            else:   
                print 'deleted negative control'
        control_pos[i] = pos_temp
        control_neg[i] = neg_temp
    del pos_temp, neg_temp

#normalize sample values based on positive and negative controls
control_pos_average = list()
control_neg_average = list()
samples_normalized = [{} for i in range(plates_len)]

#Plate averages of positive / negative controls (Default Option)
if not row_normal_use:
    for i in range(plates_len):
        if median_use:
            control_pos_average.append(median(control_pos[i].values())['median'])
            control_neg_average.append(median(control_neg[i].values())['median'])
        else:
            control_pos_average.append(mean(control_pos[i].values(), lognorm_use))
            control_neg_average.append(mean(control_neg[i].values(), lognorm_use))
        for j in range(len(samples[i])):
            samples_normalized[i][samples[i].keys()[j]] = normalize(samples[i][samples[i].keys()[j]], control_pos_average[i], control_neg_average[i], lognorm_use)
#2-Row Averages (If 'Y' in Options)
else:
    for i in range(plates_len):
        control_pos_average.append([])
            # Sorted list (by well) of controls from plate i
        control_pos_sorted = list(zip(*sorted(zip(control_pos[i].keys(), control_pos[i].values())))[1])
        control_neg_sorted = list(zip(*sorted(zip(control_neg[i].keys(), control_neg[i].values())))[1])
        samples_sorted_keys, samples_sorted_values = (list(l) for l in zip(*sorted(zip(samples[i].keys(), sample[i].values()))))
        if median_use:
            while control_pos_sorted:
                control_pos_rows = []
                for j in range(2*POS_CONTROLS_PER_ROW):
                    control_pos_rows.append(control_pos_sorted.pop(0))
                control_pos_average.[i].append(median(control_pos_rows)['median'])
            while control_neg_sorted:
                control_neg_rows = []
                for j in range(2*NEG_CONTROLS_PER_ROW):
                    control_neg_rows.append(control_neg_sorted.pop(0))
                control_neg_average.[i].append(mean(control_neg_rows)['median'])
        else:
            while control_pos_sorted:
                control_pos_rows = []
                for j in range(2*POS_CONTROLS_PER_ROW):
                    control_pos_rows.append(control_pos_sorted.pop(0))
                control_pos_average.[i].append(mean(control_pos_rows, lognorm_use))
            while control_neg_sorted:
                control_neg_rows = []
                for j in range(2*NEG_CONTROLS_PER_ROW):
                    control_neg_rows.append(control_neg_sorted.pop(0))
                control_neg_average.[i].append(mean(control_neg_rows, lognorm_use))
        while samples_sorted_values:
            j = 0
            for k in range(2*SAMPLES_PER_ROW):
                samples_normalized[i][samples_sorted_keys.pop(0)] = normalize(samples_sorted_values.pop(0), control_pos_average[i][j], control_neg_average[i][j], lognorm_use)
                j += 1

concentration = list()
for i in range(len(info)):
    for j in range(len(info[i].values())):
        concentration.append(float(info[i].values()[j][0]))
ranges = [ [7,3.5,1.75,0.875,0.4375,0.21875,0.109375,0.0546875,0.02734375,0.013671875,0.006835938,0.003417969], [7,3.5,1.75,0.875,0.4375,0.21875,0.109375,0.0546875,0.02734375,0.013671875,0.006835938,0.003417969], [15,6.9659676,3.234980307,1.502317867,0.697673172,0.323997914,0.150463932,0.069875125,0.032449857,0.015069644,0.00699831,0.00325]]
concentration = in_list(concentration)
conc_range = {}
#ranges = range(dilutions_num)
for i in range(len(concentration)):
    # print 'Enter concentration range for %s' % concentration[i]
    # for j in range(dilutions_num):
    #     ranges[j] = float(raw_input('%d: ' % (j + 1)))
    conc_range["%s" % (concentration[i])] = sorted(ranges[i], reverse=True)
del ranges

#Format Normalized Samples for GraphPad (Default Option)
if graphpad_format_use:
    graphpad = ['' for i in range(plates_len)]

    for i in range(plates_len):
        count = -1
        for j in range(len(samples_keys[i])):
                if j % (dilutions_num*3)  == 0: #Dilutions*Replicates
                    count += 1
                if j % 3 == 0:
                    graphpad[i] += (',,,'*(count))
                    graphpad[i] += '%f,' % samples_normalized[i][samples_keys[i][j]]
                elif j % 3 == 1:
                    graphpad[i] += '%f,' % samples_normalized[i][samples_keys[i][j]]
                elif j % 3 == 2:
                    graphpad[i] += '%f\n' % samples_normalized[i][samples_keys[i][j]]
    #Write GraphPad-Formatted Data to File
    file = open('%s-data.csv' % (filename.split('.')[:1][0]), 'w+')
    for i in range(len(graphpad)):
        file.seek(0, 2)
        file.write('Plate %d\n' % (i+1))
        file.write(graphpad[i])
    file.close()
else:
    file = open('%s-data.csv' % (filename.split('.')[:1][0]), 'w+')
    file.write('Plate,Well,RFU,Normalied RFU')

    if info_use:
        file.write(',Concentration,Compound,Replicate,ID')
    for i in range(plates_len):
        l = -1
        for j in range(len(samples_normalized[i])):
            file.seek(0, 2)
            file.write('\n%d,%s,%f,%f' % ((i + 1), samples_keys[i][j],samples[i][samples_keys[i][j]] ,samples_normalized[i][samples_keys[i][j]]))
            if info_use:
                if samples_keys[i][j] in info[i].keys():
                    k = 0
                    l += 1
                #print conc_range[str(float(info[i].values()[k % 6][0]))][k % 12] 
                file.write(',%f,%s,%d,%s.%d' % (conc_range[str(float(info[i].values()[l][0]))][k % 12], info[i].values()[l][1].rstrip(), (j % 3) + 1, info[i].values()[l][1].rstrip(), (j % 3) + 1))
                if j % 3 == 2:
                    k += 1

    file.close()



