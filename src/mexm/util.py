
def is_number(x):
    return isinstance(x, (int, float, complex)) and not isinstance(x, bool)
