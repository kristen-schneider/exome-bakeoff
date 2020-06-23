import matplotlib.pyplot as plt
import numpy as np
from read_data_sheet import read_data

#FULL_data_by_LP ='/Users/krsc0813/exome-bakeoff/kristen_scripts/plotting/csv_files/59_ACMG/full_name_by_lp.csv'
#FULL_data_by_CAP ='/Users/krsc0813/exome-bakeoff/kristen_scripts/plotting/csv_files/oncogenes/full_name_by_capture.csv'
LP_data='/Users/kristen/PycharmProjects/exome-bakeoff/kristen_scripts/plotting/csv_files/59_ACMG/lp_name.csv'
CAP_data='/Users/kristen/PycharmProjects/exome-bakeoff/kristen_scripts/plotting/csv_files/59_ACMG/capture_name.csv'

def main():
    # plot by library prep, color by capture
    # sample_data_lp = read_data(FULL_data_by_LP)
    # plot_sample_by_library_prep(sample_data_lp)
    
    # plot by capture, color by library prep
    #sample_data_capture = read_data(FULL_data_by_CAP)
    #plot_sample_by_capture(sample_data_capture)

    # plot and color by library prep
    lp_data = read_data(LP_data)
    combined_lp_data = combine_data(lp_data)
    plot_lp(lp_data)
    plot_combined_lb(combined_lp_data)

    # plot and color by capture
    capture_data = read_data(CAP_data)
    combined_capture_data = combine_data(capture_data)
    plot_capture(capture_data)
    plot_combined_capture(combined_capture_data)

def combine_data(sample_data):
    labels = []
    for s in sample_data[1]:
        id = s.split('-')[0]
        if id not in labels:
            labels.append(id)
    #labels = sorted(list(set(sample_data[1].split('-')[0])))

    combined_label_data = []
    for l in labels:
        one_label_data = []
        for s in range(len(sample_data[1])):
            if l in sample_data[1][s]:
                for i in sample_data[2][s]:
                    one_label_data.append(i)
                    #one_label_data.append(i for i in sample_data[2][s])
        print(len(one_label_data))
        combined_label_data.append(one_label_data)

    return [labels, combined_label_data]



def plot_sample_by_library_prep(sample_data):
    plt.figure(figsize=(40, 15))
    box = plt.boxplot(sample_data[2], patch_artist=True)

    # colors
    colors = []
    Agilent = 'green'
    IDT = 'blue'
    Roche = 'black'
    Twist = 'yellow'
    for s in sample_data[1]:
        capture = s.split('-')[1]
        if capture == 'Agilent': colors.append(Agilent)
        elif capture == 'IDT': colors.append(IDT)
        elif capture == 'Roche': colors.append(Roche)
        elif capture == 'Twist': colors.append(Twist)
        else: print('no capture tech found')

    # label name
    x_label = []
    for s in sample_data[1]:
        x_label.append(s.split('-')[0])


    for b in range(len(box['boxes'])):
        plt.setp(box['boxes'][b], facecolor=colors[b])

    plt.plot([], c='green', label='Agilent')
    plt.plot([], c='blue', label='IDT')
    plt.plot([], c='black', label='Roche')
    plt.plot([], c='yellow', label='Twist')
    plt.legend()

    plt.title('Quality by Library Prep\n(colored by Capture Technology)')
    plt.xticks(np.arange(1, len(sample_data[1]) + 1), x_label, rotation=90)
    plt.xlim(-1, len(sample_data[1]) + 1)
    plt.savefig('sample_data_lp.png')

def plot_sample_by_capture(sample_data):
    plt.figure(figsize=(40, 15))
    box = plt.boxplot(sample_data[2], patch_artist=True)

    # colors
    colors = []
    AgilentQXT = 'lightpink'
    KAPACovaris = 'lightblue'
    KAPAEnzyme = 'lightyellow'
    Nextera = 'lightsalmon'
    NexteraCap = 'lightgreen'
    Twist = 'plum'
    for s in sample_data[1]:
        capture = s.split('-')[0]
        if capture == 'AgilentQXT': colors.append(AgilentQXT)
        elif capture == 'KAPACovaris': colors.append(KAPACovaris)
        elif capture == 'KAPAEnzyme': colors.append(KAPAEnzyme)
        elif capture == 'Nextera': colors.append(Nextera)
        elif capture == 'NexteraCap': colors.append(NexteraCap)
        elif capture == 'Twist': colors.append(Twist)
        else: print('no capture tech found')

    # label name
    x_label = []
    for s in sample_data[1]:
        x_label.append(s.split('-')[1])


    for b in range(len(box['boxes'])):
        plt.setp(box['boxes'][b], facecolor=colors[b])

    plt.plot([], c='lightpink', label='AgilentQXT')
    plt.plot([], c='lightblue', label='KAPACovaris')
    plt.plot([], c='lightyellow', label='KAPAEnzyme')
    plt.plot([], c='lightsalmon', label='Nextera')
    plt.plot([], c='lightgreen', label='NexteraCap')
    plt.plot([], c='plum', label='Twist')
    plt.legend()

    plt.title('Quality by Capture\n(colored by Library Prep Technology)')
    plt.xticks(np.arange(1, len(sample_data[1]) + 1), x_label, rotation=90)
    plt.xlim(-1, len(sample_data[1]) + 1)
    plt.savefig('sample_data_capture.png')

def plot_lp(lp_data):
    plt.figure(figsize=(40, 15))
    box = plt.boxplot(lp_data[2], patch_artist=True)

    colors = ['lightpink'] * 12
    for i in range(4): colors.append('lightblue')
    for i in range(12): colors.append('lightyellow')
    for i in range(16): colors.append('lightsalmon')
    for i in range(16): colors.append('lightgreen')
    for i in range(8): colors.append('plum')


    for b in range(len(box['boxes'])):
        plt.setp(box['boxes'][b], facecolor=colors[b])

    plt.title('Quality by Library Prep')
    plt.xticks(np.arange(1, len(lp_data[1])+1), lp_data[1], rotation=90)
    plt.xlim(-1, len(lp_data[1])+1)
    plt.savefig('lib_prep_onco.png')

def plot_combined_lb(combined_lp_data):
    plt.figure(figsize=(40, 15))
    box = plt.boxplot(combined_lp_data[1], patch_artist=True)

    colors = ['lightpink', 'lightblue', 'lightyellow', 'lightsalmon', 'lightgreen', 'plum']

    for b in range(len(box['boxes'])):
        plt.setp(box['boxes'][b], facecolor=colors[b])

    plt.title('Quality by Library Prep')
    plt.xticks(np.arange(1, len(combined_lp_data[0])+1), combined_lp_data[0], rotation=90)
    plt.xlim(0, len(combined_lp_data[0])+1)
    plt.savefig('lib_prep_combined_onco.png')

def plot_capture(capture_data):
    plt.figure(figsize=(40, 15))
    box = plt.boxplot(capture_data[2], patch_artist=True)

    colors = ['green'] * 8
    for i in range(20): colors.append('blue')
    for i in range(24): colors.append('black')
    for i in range(16): colors.append('yellow')


    for b in range(len(box['boxes'])):
        plt.setp(box['boxes'][b], facecolor=colors[b])



    plt.title('Quality by Capture')
    plt.xticks(np.arange(1, len(capture_data[1])+1), capture_data[1], rotation=90)
    plt.xlim(-1, len(capture_data[1])+1)
    plt.savefig('capture_onco.png')


def plot_combined_capture(combined_capture_data):
    plt.figure(figsize=(40, 15))
    box = plt.boxplot(combined_capture_data[1], patch_artist=True)

    colors = ['green', 'blue', 'black', 'yellow']

    for b in range(len(box['boxes'])):
        plt.setp(box['boxes'][b], facecolor=colors[b])


    plt.title('Quality by Capture')
    plt.xticks(np.arange(1, len(combined_capture_data[0])+1), combined_capture_data[0], rotation=90)
    plt.xlim(0, len(combined_capture_data[0])+1)
    plt.savefig('capture_combined_onco.png')


main()
