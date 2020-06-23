def alphabetize_full(full_dict):
    sorted_names = sorted(full_dict)
    sorted_data = []
    
    for i in sorted_names:
        name_data = [i, full_dict[i]]
        sorted_data.append(name_data)
    
    return sorted_data   

def alphabetize_prep(prep_dict):
    sorted_names = sorted(prep_dict)
    sorted_data = []

    for i in sorted_names: 
        name_data = [i, prep_dict[i]]
        sorted_data.append(name_data)
    
    return sorted_data 
    
def alphabetize_capture(capture_dict):
    sorted_names = sorted(capture_dict)
    sorted_data = []
    
    for i in sorted_names:
        name_data = [i, capture_dict[i]]
        sorted_data.append(name_data)

    return sorted_data
