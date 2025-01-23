def find_index(header, filepath, channel):
    with open(filepath, 'r') as file:
        for i in range(header+1):
            line = file.readline()
    split_line = line.split('\t')
    split_line.remove('\n')
    return split_line.index(channel)

def get_parameter_values(parameters_dict, parameter) -> list:
    value = None
    if parameter in parameters_dict:
        value = parameters_dict[parameter]
    else:
        print('Parameter not in Parameter dict!')
    return value

def convert_header_to_dict(filepath, separator=':', header_indicator='#', tags_list:list=None) -> dict:
    try:
        header = _read_parameters_txt(filepath, header_indicator, tags_list)
    except:
        print('could not read header')
        return None
    # print('header:', header)
    parameters_dict = {}
    for i in range(len(header)):
        key, value = _simplify_line(header[i], separator, header_indicator, tags_list)
        # print(f'key: <{key}>' , f'value: <{value}>')
        # check if the key is in the tags list
        if key in tags_list:
            parameters_dict[key] = value
        else:
            return None
        if key is None:
            return None
    # check if every tag is in the dict
    for tag in tags_list:
        if tag not in parameters_dict:
            return None
    return parameters_dict

def _read_parameters_txt(filepath, header_indicator, tags) -> list:
    content = []
    with open(filepath, 'r', encoding='UTF-8') as file:
        all = file.read()
    # split content into lines
    # print('all:', all)
    all_lines = all.split('\n')
    # print('all_lines:', all_lines)
    # print('tags:', tags)
    # check if the lines contain on of the tags
    for line in all_lines:
        for tag in tags:
            if type(tag) == str:
                if tag in line:
                    content.append(line)
            elif type(tag) == list:
                for subtag in tag:
                    if subtag in line:
                        content.append(line)
    # print('content:', content)
    return content

def _simplify_line(line, separator, header_indicator, tags) -> tuple:
    # replace the header indicator if it is not ''
    if header_indicator != '':
        line = line.replace(header_indicator, '')
    # split line at separator
    if separator != '':
        line = line.split(separator)
        if len(line) <= 1:
            return None, None
        # the first element is the key and should be identical to one of the tags
        # remove linebreak
        try: line[1] = line[1].replace(u'\n', '')
        except: pass
        # remove tabs and empty spaces from second element
        line[1] = _remove_empty_spaces(line[1])
        # split second element into list if possible:
        try: line[1] = line[1].split(u'\t')
        except: pass
        # split second element into list if possible:
        # remove empty elements in second list
        try: line[1] = list(filter(('').__ne__, line[1])) # the date for example might contain a space inbetween date and time
        except: pass
        if len(line[1]) == 1:
            line[1] = line[1][0]
        # try to simplyfy line[0] to only contain the tags, sometimes to much spaces might be in the line
        # however sometimes shorter tags might be part of longer tags, such as 'SCAN' and 'SCAN AREA'...
        for tag in tags:
            if tag in line[0]:
                # if there is a chracter around the tag it is not the tag we are looking for, optimally we want to make sure that the tag is sourrounded by 2 spaces or the end of the line
                # print('tag:', tag)
                # print('line[0]:', line[0])
                start_index = line[0].index(tag)
                end_index = start_index + len(tag)
                if start_index > 1 and line[0][start_index-2] != ' ':
                    continue
                elif start_index > 0 and line[0][start_index-1] != ' ':
                    continue
                if end_index < len(line[0])-1 and line[0][end_index+1] != ' ':
                    continue
                elif end_index < len(line[0]) and line[0][end_index] != ' ':
                    continue
                # if len(line[0]) > len(tag):
                #     if line[0][line[0].index(tag)-1] != ' ' or line[0][line[0].index(tag)+len(tag)] != ' ':
                #         if line[0][line[0].index(tag)-2] != '\t' or line[0][line[0].index(tag)+len(tag)] != '\t':
                #             continue
                line[0] = tag
                # print('line[0]:', line[0])

                break
        key = line[0]
        value = line[1]

    else:
        try: line = line.replace(u'\xa0', '')
        except: pass
        try: line = line.replace(u'\t\t', '\t') # neaspec formatting sucks thus sometimes a simple \t is formatted as \t\xa0\t
        except: pass
        line = line.split('\t')
        # remove empty elements
        line = list(filter(('').__ne__, line))
        # if line has more than 2 elements the first element is the key and the other elements as a list are the value
        if len(line) > 2:
            key = line[0]
            value = line[1:]
        else:
            key = line[0]
            value = line[1]
    return key, value

def _remove_empty_spaces(line) -> str:
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

