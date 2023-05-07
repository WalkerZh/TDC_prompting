import argparse

# from tdc.utils import retrieve_label_name_list
# label_list = retrieve_label_name_list('QM9')

# print(label_list)

from tdc.generation import RetroSyn
data = RetroSyn(name = 'USPTO-50K')
split = data.get_split()

print(split['valid'])

# print(split['test']['Reaction'][14])
# print(type(split['test']['Y'][1]))

# from tdc.utils import get_label_map
# print(get_label_map(name = 'TWOSIDES', task = 'DDI', name_column = 'Side Effect Name'))