import numpy as np
from matplotlib import pyplot as plt



def calculate_squared_mean_deviation(array_1, array_2) -> float:
    # check if both arrays are of equal length
    if len(array_1) != len(array_2):
        print('The lengths of the arrays must be equal!')
        return 0
    mean_deviation = 0
    N = len(array_1)
    for i in range(N):
        mean_deviation += (array_1[i]**2-array_2[i]**2)/N
    return mean_deviation

def calculate_squared_deviation(array_1, array_2) -> list:
    # check if both arrays are of equal length
    if len(array_1) != len(array_2):
        print('The lengths of the arrays must be equal!')
        return 0
    mean_deviation = []
    N = len(array_1)
    for i in range(N):
        # ignore '0' values
        if array_1[i] == 0 or array_2[i] == 0:
            mean_deviation.append(0)
        else:
            mean_deviation.append(abs(array_1[i]-array_2[i]))
    return mean_deviation

def _shift_array_1d_by_index(array_1, array_2, index) -> tuple:
    """This function shifts the second array with respect to the first array. The first array is modified to keep the same length.

    Args:
        array_1 (np.array): first and reference array
        array_2 (np.array): second array
        index (int): how much should be shifted, to the left if negative
    """
    # print(f'in shift 1d, index: {index}')
    if index > 0:
        array_1_new = array_1
        array_2_new = array_2
        for i in range(abs(index)):
            # print(f'apply shift, i: {i}')
            array_1_new = np.append(array_1_new, 0)
            array_2_new = np.insert(array_2_new, 0, 0)
    elif index < 0:
        array_1_new = array_1
        array_2_new = array_2
        for i in range(abs(index)):
            array_1_new = np.insert(array_1_new, 0, 0)
            array_2_new = np.append(array_2_new, 0)
    else:
        array_1_new = array_1
        array_2_new = array_2
    return array_1_new, array_2_new

def shift_array_2d_by_index(array_1, array_2, index) -> tuple:
    """This function shifts the second 2d array relative to the first. Zeroes will be created on the outside.

    Args:
        array_1 (np.ndarray): array1
        array_2 (np.ndarray): array2
        index (int): how much to shift, negative means to the left

    Returns:
        tuple: shifted arrays
    """
    x_res = len(array_1[0])
    y_res = len(array_1)
    # print(f'xres={x_res} yres={y_res}, index={index}')
    x_res_new = x_res + abs(index)
    #create new lists with modified resolution
    array_1_new = np.zeros((y_res, x_res_new))
    array_2_new = np.zeros((y_res, x_res_new))
    for y in range(y_res):
        test_1, test_2 = _shift_array_1d_by_index(array_1[y], array_2[y], index)
        # print(f'len test1: {len(test_1)}, len test 2: {len(test_2)}')
        array_1_new[y] = test_1
        array_2_new[y] = test_2
        # array_1_new[y], array_2_new[y] = _shift_array_1d_by_index(array_1[y], array_2[y], index)
    return array_1_new, array_2_new

def create_mean_array(array_1, array_2):
    #check if dimensions are identical
    if len(array_1) != len(array_2) and len(array_1[0]) != len(array_2[0]):
        print('The length of the given arrays is not identical.')
        exit()
    else:
        x_res = len(array_1[0])
        y_res = len(array_1)
        new_array = np.zeros((y_res, x_res))
        for y in range(y_res):
            for x in range(x_res):
                new_array[y][x] = (array_1[y][x] + array_2[y][x])/2
        return new_array

def create_mean_array_V2(array_1, array_2, index):
    """This variant is meant to keep the size of the original array!

    Args:
        array_1 (_type_): _description_
        array_2 (_type_): _description_

    Returns:
        _type_: _description_
    """
    #check if dimensions are identical
    if len(array_1) != len(array_2) and len(array_1[0]) != len(array_2[0]):
        print('The length of the given arrays is not identical.')
        exit()
    else:
        x_res = len(array_1[0])
        y_res = len(array_1)
        new_array = np.zeros((y_res, x_res))
        for y in range(y_res):
            if index > 0:
                for x in range(x_res - index):
                    new_array[y][x] = (array_1[y][x+index] + array_2[y][x])/2
            if index < 0:
                for x in range(x_res - abs(index)):
                    new_array[y][x] = (array_1[y][x] + array_2[y][x+abs(index)])/2
        return new_array

def minimize_deviation_1d(array_1, array_2, n_tries=5, display=True) -> int:
    """This function tries to find the optimal shift between two arrays to find the lowest deviation. Best to use leveled data.

    Args:
        array_1 (np.array): first array, typically height data
        array_2 (np.array): second array, typically height data
        n_tries (int, optional): the maximum shift expected for optimal overlap, will be applied symmetrically to right and left shift. Defaults to 5.

    Returns:
        int: the optimal shift index, left shift if negative
    """
    # try shifting to the right
    mean_deviations_right = []
    for i in range(n_tries):
        array_1_new, array_2_new = _shift_array_1d_by_index(array_1, array_2, i+1)
        mean_deviation_array = calculate_squared_deviation(array_1_new, array_2_new)
        mean_deviation = np.mean(mean_deviation_array)
        mean_deviations_right.append(mean_deviation)
        x = range(len(array_1_new))
        # plt.plot(x, array_1_new, label='array_2')
        # plt.plot(x, array_2_new, label='array_1')
        # plt.plot(x, mean_deviation_array, label="Mean deviation_array")
        # plt.hlines(mean_deviation, label="Mean deviation", xmin=0, xmax=len(array_1))
        # plt.legend()
        # plt.show()
    min_dev_right = min(mean_deviations_right)

    # try shifting to the left
    mean_deviations_left = []
    for i in range(n_tries):
        array_1_new, array_2_new = _shift_array_1d_by_index(array_1, array_2, -(i+1))
        mean_deviation_array = calculate_squared_deviation(array_1_new, array_2_new)
        mean_deviation = np.mean(mean_deviation_array)
        mean_deviations_left.append(mean_deviation)
        x = range(len(array_1_new))
        # plt.plot(x, array_1_new, label='array_2')
        # plt.plot(x, array_2_new, label='array_1')
        # plt.plot(x, mean_deviation_array, label="Mean deviation_array")
        # plt.hlines(mean_deviation, label="Mean deviation", xmin=0, xmax=len(array_1))
        # plt.legend()
        # plt.show()
    min_dev_left = min(mean_deviations_left)

    # try not shifting at all
    mean_deviation_array = calculate_squared_deviation(array_1, array_2)
    min_dev_center = np.mean(mean_deviation_array)
    
    min_dev_index = 0
    if min_dev_left <= min_dev_right and min_dev_left < min_dev_center:
        # print(f'The minimal deviation of {min_dev_left} is reached for a shift of {mean_deviations_left.index(min_dev_left)+1} to the left')
        array_1_new, array_2_new = _shift_array_1d_by_index(array_1, array_2, -(mean_deviations_left.index(min_dev_left)+1))
        min_dev_index = -(mean_deviations_left.index(min_dev_left)+1)
    elif min_dev_left > min_dev_right and min_dev_right < min_dev_center:
        # print(f'The minimal deviation of {min_dev_right} is reached for a shift of {mean_deviations_right.index(min_dev_right)+1} to the right')
        array_1_new, array_2_new = _shift_array_1d_by_index(array_1, array_2, +(mean_deviations_right.index(min_dev_right)+1))
        min_dev_index = (mean_deviations_right.index(min_dev_right)+1)
    elif min_dev_center <= min_dev_left and min_dev_center <= min_dev_right:
        array_1_new = array_1
        array_2_new = array_2
        min_dev_index = 0
    else:
        print(f'min_dev_left: {min_dev_left}, min_dev_center: {min_dev_center}, min_dev_right: {min_dev_right}')
        print('Oops, something went wrong, could not do the optimization')
        
    
    if display == True:
        # redo the best deviation:
        x = range(len(array_1_new))
        plt.plot(x, array_1_new, label='array_1')
        plt.plot(x, array_2_new, label='array_2')
        mean_deviation_array = calculate_squared_deviation(array_1_new, array_2_new)
        mean_deviation = np.mean(mean_deviation_array)
        plt.plot(x, mean_deviation_array, label="Mean deviation_array")
        plt.hlines(mean_deviation, label="Mean deviation", xmin=0, xmax=len(array_1_new))
        plt.legend()
        plt.show()

    return min_dev_index

def minimize_deviation_2d(array_1, array_2, n_tries=5, display=True) -> int:
    """This function tries to find the optimal shift between two arrays to find the lowest deviation. Best to use leveled data.

    Args:
        array_1 (np.array): first array, typically height data
        array_2 (np.array): second array, typically height data
        n_tries (int, optional): the maximum shift expected for optimal overlap, will be applied symmetrically to right and left shift. Defaults to 5.

    Returns:
        int: the optimal shift index, left shift if negative
    """
    y_shifts = [] # contains the optimal shifts per y data line
    rows = len(array_1)
    for i in range(rows):
        y_shifts.append(minimize_deviation_1d(array_1[i], array_2[i], n_tries, display))
    # plt.plot(range(rows), y_shifts, label="optimal dev index per row")
    # plt.legend(loc="upper left")
    # plt.show()
    return int(round(np.mean(y_shifts)))

def _shift_data(data, y_shifts):
    min_shift = min(y_shifts)
    max_shift = max(y_shifts)
    # print('shifts: ', y_shifts)
    if min_shift < 0:
        min_shift = abs(min_shift)
    else:
        min_shift = 0
    if max_shift > 0:
        max_shift = max_shift
    else:
        max_shift = 0
    # min and max shift now correspond to the empty spaces which have to be added in the beginning(min) and end(max) of the first line
    total_shift = min_shift +  max_shift
    if total_shift != 0:
        y_res = len(data)
        x_res = len(data[0])
        x_res_new = x_res + min_shift + max_shift
        data_new = np.zeros((y_res, x_res_new))
        for i in range(len(y_shifts)):
            row_data = data[i]
            index = y_shifts[i]
            num_of_inserts = min_shift + index
            num_of_appends = max_shift - index
            if num_of_inserts > 0:
                for j in range(num_of_inserts):
                    row_data = np.insert(row_data, 0, 0)
            if num_of_appends > 0:
                for j in range(num_of_appends):
                    row_data = np.append(row_data, 0)
            data_new[i] = row_data
    else:
        return data
    '''
    max_shift = abs(max(y_shifts, key=abs))
    print('shifts: ', y_shifts)
    print('max shift: ', max_shift)
    y_res = len(data)
    x_res = len(data[0])
    x_res_new = x_res + 2*max_shift
    data_new = np.zeros((y_res, x_res_new))
    for i in range(len(y_shifts)):
        row_data = data[i]
        index = y_shifts[i]
        num_of_appends = max_shift - index
        num_of_inserts = max_shift + index
        if num_of_inserts > 0:
            for j in range(num_of_inserts):
                row_data = np.insert(row_data, 0, 0)
        if num_of_appends > 0:
            for j in range(num_of_appends):
                row_data = np.append(row_data, 0)
        data_new[i] = row_data
    '''    
    return data_new

def minimize_drift(data:np.array, n_tries=5, display=False) -> np.array:
    """This function tries to shift subsequent lines of a single measurement to minimize the deviation.
    The function will return the drift corrected data.

    Returns:
        np.array: the drift corrected data
    """
    rows = len(data)
    # cols = len(data[0])
    y_shifts = [0] # the first row will not be shifted, each subsequent shift is related to the first line
    for y in range(rows-1):
        y_shifts.append(minimize_deviation_1d(data[0], data[y+1], n_tries, display))
    return _shift_data(data, y_shifts)




