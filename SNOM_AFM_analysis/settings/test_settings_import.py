import json
from pathlib import Path, PurePath
this_files_path = Path(__file__).parent
from snom_colormaps import*

channel = 'O2A'


def replace_placeholder(dictionary):
    colormaps = {"<SNOM_amplitude>": SNOM_amplitude,
                 "<SNOM_height>": SNOM_height,
                 "<SNOM_phase>": SNOM_phase,
                 "<SNOM_realpart>": SNOM_realpart}
    placeholders = {'<channel>': channel}
    # first iterate through all placeholders and replace them in the dictionary
    for placeholder in placeholders:
        value = placeholders[placeholder]
        for key in dictionary:
            if placeholder in dictionary[key]:
                dictionary[key] = dictionary[key].replace(placeholder, value)
                # print('replaced channel!')
    # replace colormaps
    for key in dictionary:
        for colormap in colormaps:
            print(colormap, type(colormap))
            print(dictionary[key])
            if colormap in dictionary[key]:
                dictionary[key] = colormaps[colormap]
                break
    return dictionary

def test():
    parameters_path = this_files_path / Path('plotting_parameters.json')

    with open(parameters_path, 'r') as file:
        plotting_parameters = json.load(file)


    # print(plotting_parameters)
    # print(plotting_parameters['amplitude_title'])
    
    # try to replace placeholder in json with actual channel value

    plotting_parameters = replace_placeholder(plotting_parameters)

    # print(plotting_parameters['amplitude_title'])
    print(plotting_parameters)


    '''
    data = {
        'first': 'John',
        'last': 'Doe',
        'kids': {
            'first': 'number1',
            'second': 'number2',
            'third': 'number3'
        }
    }  

    print(data)
    print('{first} {last}, {kids[first]}, {kids[second]}, {kids[third]}'.format(**data))
    print(f'{data["first"]} {data["last"]}, {data["kids"]["first"]}, {data["kids"]["second"]}, {data["kids"]["third"]}')

    '''

def main():
    test()


if __name__ == '__main__':
    main()


