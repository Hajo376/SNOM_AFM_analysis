
def _Find_Header_Size(filepath) -> int:
    beginning_symb = '#'
    inside_header = True
    header = 0
    max_number = 100 # to avoid endless loop if beginning symbol changes
    with open(filepath, 'r') as file:
        while inside_header and header < max_number:
            line = file.readline()
            if line[0:1] == beginning_symb:
                header += 1
            else: inside_header = False
    return header

def Find_Index(header, filepath, channel):
    with open(filepath, 'r') as file:
        for i in range(header+1):
            line = file.readline()
    # print(line)
    split_line = line.split('\t')
    split_line.remove('\n')
    # print(split_line)
    return split_line.index(channel)

def _Remove_Empty_Spaces(line) -> str:
    # print('starting to replace empty spaces')
    try:
        line = line.replace(u'\xa0', '')
    except: pass
    try:
        line = line.replace(u'\t\t', '\t') # neaspec formatting sucks thus sometimes a simple \t is formatted as \t\xa0\t
    except:
        pass
    # seems like all lines have additional \t in front, so lets get rid of that
    try:
        line = line.replace(u'\t', '', 1)
    except: pass
    return line


def _Simplify_Line(line):
    # replace # in the beginning, might be different for different sofrware versions
    # print(line)
    line = line.replace('# ', '')
    # print(line)
    line = line.replace(':', 'split-here', 1) # only replace first occurence
    # print(line)
    line = line.split('split-here')
    # print(line)
    # only remove empty spaces only from second element, keep only the important empty spaces in second one as it is used as a separator for mutliple values
    line[1] = _Remove_Empty_Spaces(line[1])
    # print(line)
    line[1] = line[1].replace(u'\n', '')
    # print(line)
    # split second element into list if possible:
    try:
        line[1] = line[1].split(u'\t')
    except: pass
    # sometimes an empty element will be created, get rid of that too:
    if type(line[1]) is list and '' in line[1]:
        line[1].remove('')
    key = line[0]
    value = line[1]

    return key, value

def Convert_Header_To_Dict(filepath) -> dict:
    content = _Read_Parameters_Txt(filepath)
    parameters_dict = {}
    for i in range(len(content)):
        if i == 0:
            # first line is just the website adress
            pass
        else:
            key, value = _Simplify_Line(content[i])
            # sort dictornary entries correctly:
            # key: [value1, value2, ... , '[unit]'] but only if value contains a list with len > 1 and if value[0] is a string that cannot be converted to a number
            new_value = value.copy()
            if type(value) == list:
                # print('encountered a list', value)
                if len(value) > 1:
                    # print('len(list)', len(value))
                    try: float(value[0])
                    except: exception = True
                    else: exeption = False
                    # print('exception: ', exception)
                    if exception:
                        # print('exception found')
                        for i in range(len(value)-1):
                            new_value[i] = value[i+1]
                        new_value[-1] = value[0]
                        # print('old value: ', value, '\nnew value: ', new_value)


            
            # parameters_dict[key] = value
            parameters_dict[key] = new_value
    # for key in parameters_dict:
    #     print(f'{key}: {parameters_dict[key]}')
    return parameters_dict

def _Read_Parameters_Txt(filepath) -> list:
    content = []
    # get length of header, header lines should beginn with a '#'
    header = _Find_Header_Size(filepath)
    with open(filepath, 'r', encoding='UTF-8') as file:
        for i in range(header):
            line = file.readline()
            content.append(line)
    return content

def Get_Parameter_Values(parameters_dict, parameter) -> list:
    value = None
    if parameter in parameters_dict:
        value = parameters_dict[parameter]
    else:
        print('Parameter not in Parameter dict!')
    return value
