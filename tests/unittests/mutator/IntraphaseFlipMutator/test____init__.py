import pytest

from pymatmc2.mutator import IntraphaseFlipMutator

def dev__implements_base_class():
    obj = IntraphaseFlipMutator()

    base_method_names = [
        'mutate_multicell',
        'accept_or_reject'
    ]

    for base_method_name in base_method_names:
        print(base_method_name)
    

if __name__ == "__main__":
    dev__implements_base_class()
    