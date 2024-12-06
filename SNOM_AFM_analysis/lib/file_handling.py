
def _Find_Header_Size(filepath, header_indicator) -> int:
    # header_indicator = '#' # sadly header lines are not always starting with '#'
    inside_header = True
    header = 0
    max_number = 100 # to avoid endless loop if beginning symbol changes
    with open(filepath, 'r') as file:
        while inside_header and header < max_number:
            line = file.readline()
            if line[0:1] == header_indicator:
                header += 1
            else: inside_header = False
    if header == 0:
        # in case no line is starting with '#' just count the lines in the document
        header = sum(1 for line in open(filepath))
        if header > max_number:
            exit()
    return header

# def _Find_Header_Size_V2(filepath, header_indicator, tags) -> int:
#     pass

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

def _simplify_line_v2(line, separator, header_indicator, tags) -> tuple:
    # replace the header indicator if it is not ''
    # print('unsimplified line: ', line)
    if header_indicator != '':
        line = line.replace(header_indicator, '')
    # split line at separator
    if separator != '':
        line = line.split(separator)
        if len(line) <= 1:
            return None, None
        # print('simpli#fied line: ', line)
        # the first element is the key and should be identical to one of the tags
        # remove linebreak
        try: line[1] = line[1].replace(u'\n', '')
        except: pass
        # remove tabs and empty spaces from second element
        line[1] = _Remove_Empty_Spaces(line[1])
        # print('simplified line: ', line)
        # split second element into list if possible:
        try: line[1] = line[1].split(u'\t')
        except: pass
        # split second element into list if possible:
        # line[1] = line[1].split(u' ')
        # print('simplified line: ', line)
        # remove empty elements in second list
        try: line[1] = list(filter(('').__ne__, line[1])) # the date for example might contain a space inbetween date and time
        except: pass
        # print('simplified line: ', line)
        if len(line[1]) == 1:
            line[1] = line[1][0]
        # print('simplified line: ', line)
        key = line[0]
        value = line[1]

    else:
        try: line = line.replace(u'\xa0', '')
        except: pass
        try: line = line.replace(u'\t\t', '\t') # neaspec formatting sucks thus sometimes a simple \t is formatted as \t\xa0\t
        except: pass
        # print('simplified line: ', line)
        line = line.split('\t')
        # print('simplified line: ', line)
        # if len(line) <= 1:
        # remove empty elements
        line = list(filter(('').__ne__, line))
        # print('simplified line: ', line)

        # if line has more than 2 elements the first element is the key and the other elements as a list are the value
        if len(line) > 2:
            key = line[0]
            value = line[1:]
        else:
            key = line[0]
            value = line[1]
    # print('simplified line: ', line)
    # print('found key value pair: ', key, value)
    return key, value
    


def _get_version_number(content) -> str:
    # print('content: ', content)
    version = ''
    for i in range(len(content)):
        # print(content[i])
        if 'Version' in content[i]:
            # print(content[i])
            line = content[i].replace(u'\xa0', u' ') # replace all non-breaking spaces with normal spaces
            line = line.replace(u'\t', '') # remove all tabs
            line = line.replace(u'\n', '') # remove all linebreaks
            line = line.split(' ') # split line into list
            # print('line: ', line)
            version = line[-1]
            # print('version number: ' + version)
            break
    return version

def Convert_Header_To_Dict(filepath, separator=':', header_indicator='#', tags:list=None) -> dict:
    content = _Read_Parameters_Txt(filepath, header_indicator)
    # try to get version number
    print('trying to get version number')
    _get_version_number(content)

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

def convert_header_to_dict_v2(filepath, separator=':', header_indicator='#', tags_list:list=None) -> dict:
    # print('trying to convert header to dict')
    # print('separator: ', separator)
    # print('header_indicator: ', header_indicator)
    # tags_list = list(tags.values())
    try:
        header = _read_parameters_txt_v2(filepath, header_indicator, tags_list)
    except:
        print('could not read header')
        return None
    parameters_dict = {}
    for i in range(len(header)):
        key, value = _simplify_line_v2(header[i], separator, header_indicator, tags_list)
        # check if the key is in the tags list
        if key in tags_list:
            parameters_dict[key] = value
        else:
            print(f'key <{key}> not in tags list')
            return None
        if key is None:
            return None
    # check if every tag is in the dict
    for tag in tags_list:
        if tag not in parameters_dict:
            print(f'tag <{tag}> not in parameters_dict')
            return None
    return parameters_dict

def _Read_Parameters_Txt(filepath, header_indicator) -> list:
    content = []
    # get length of header, header lines should beginn with a '#'
    header = _Find_Header_Size(filepath, header_indicator)
    with open(filepath, 'r', encoding='UTF-8') as file:
        for i in range(header):
            line = file.readline()
            content.append(line)
    return content

def _read_parameters_txt_v2(filepath, header_indicator, tags) -> list:
    content = []
    with open(filepath, 'r', encoding='UTF-8') as file:
        all = file.read()
    # split content into lines
    all_lines = all.split('\n')
    # print('all_lines: ', all_lines)
    # convert tags dict to list of tag values
    # tags = list(tags.values())
    # print('tags: ', tags)
    # print('type(tags): ', type(tags))
    # check if the lines contain on of the tags
    for line in all_lines:
        # print('line: ', line)
        for tag in tags:
            # print('tag: ', tag)
            # print('tag.type(): ', type(tag))
            if type(tag) == str:
                if tag in line:
                    content.append(line)
            elif type(tag) == list:
                for subtag in tag:
                    if subtag in line:
                        content.append(line)
    # print('content: ', content)
    return content

def Get_Parameter_Values(parameters_dict, parameter) -> list:
    value = None
    if parameter in parameters_dict:
        value = parameters_dict[parameter]
    else:
        print('Parameter not in Parameter dict!')
    return value
