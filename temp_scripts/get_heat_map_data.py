import sys
import os

final_metrics = sys.argv[1]

def main():
    all_data
    for metric_file in os.listdir(final_metrics):
        if metric_file.endswith('quality_final.txt'):
            file_data = read_file(metric_file)
            

def read_file(metric_file):
    f = open(os.path.join(final_metrics_path, metrics_txt), 'r')
    data = []
    for line in f:
        line = line.strip().split()
        data.append(float(line[1].strip()))
    return data


# heatmap metric for quality
def get_avg_residual(data):
    avg = np.average(data)
    all_residuals = []
    for i in data:
        residual = np.absolute(avg - i)
        all_residuals.append(residual)
    #print (np.average(all_residuals))
    return np.average(all_residuals)

# heatmap metric for strand bias
def get_proportion_outside(data):
    # lower bound
    forward = get_forward_strand_mean_std(data)
    # upper bound 
    reverse = get_reverse_strand_mean_std(data)

    outside = 0
    for i in data:
        if i < forward or i > reverse: outside += 1
    #print("strandbias: ", outside/len(y))
    try: proportion_outside = outside/len(data)
    except ZeroDivisionError: proportion_outside = -1
    return proportion_outside

def get_forward_strand_mean_std(data):
    forwards = []
    for i in data:
        if i < 0.5: forwards.append(i)
    return np.average(forwards)-2*np.std(forwards)
def get_reverse_strand_mean_std(data):
    reverses = []
    for i in data:
        if i >= 0.5: reverses.append(i)
    return np.average(reverses)+2*np.std(reverses)

