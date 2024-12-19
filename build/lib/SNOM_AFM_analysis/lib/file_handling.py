def Find_Index(header, filepath, channel):
    with open(filepath, 'r') as file:
        for i in range(header+1):
            line = file.readline()
    split_line = line.split('\t')
    split_line.remove('\n')
    return split_line.index(channel)

def Get_Parameter_Values(parameters_dict, parameter) -> list:
    value = None
    if parameter in parameters_dict:
        value = parameters_dict[parameter]
    else:
        print('Parameter not in Parameter dict!')
    return value

def convert_header_to_dict(filepath, separator=':', header_indicator='#', tags_list:list=None) -> dict:
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
            # print(f'key <{key}> not in tags list')
            return None
        if key is None:
            return None
    # check if every tag is in the dict
    for tag in tags_list:
        if tag not in parameters_dict:
            # print(f'tag <{tag}> not in parameters_dict')
            return None
    return parameters_dict

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

