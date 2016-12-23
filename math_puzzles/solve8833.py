import itertools as it
import pandas as pd


num_sets = {"": [8., 8., 3., 3.]}


def custom_math(pair, func):
    """
    Custom math function that deals with division by zero

    :param pair: tuple len = 2
        floats to combine by some function
    :param func: str
        str that is "__add__", "__sub__", "__div__", or "__mul__"
    :return: float
        the resulting value as a float
    """
    pair = [x for x in pair]
    if func == "__div__":
        func = "__mul__"
        try:
            pair[1] = 1. / pair[1]
        except ZeroDivisionError:
            pair[1] = 0
    new = getattr(pair[0], func)(pair[1])
    return new


def cycle_num_set(num_set):
    """
    Takes in a dictionary of lists which to combine two numbers from by some math

    :param num_set: dict of lists
        dict w/ key of previous operations and list of resulting numbers
    :return: dict of lists
        dict with all possible combinations of two combinbined numbers from prev list
    """
    methods = ["__add__", "__sub__", "__div__", "__mul__"]
    new_num_set = dict()
    for k in num_set.keys():
        unique_pairs = list(set(it.permutations(num_set[k], 2)))
        for func, pair in it.product(methods, unique_pairs):
            new_set = [x for x in num_set[k]]
            for i in pair:
                new_set.remove(i)
            new_key = "({0}{1}{2})".format(pair[0], func, pair[1])
            new_num_set[k + new_key] = [custom_math(pair, func)] + new_set
    return new_num_set


def cycle8833():
    """
    Get all possible combinations from combining 8, 8, and 3
    """
    num_set  = {"": [8., 8., 3., 3.]}
    for i in range(3):
        num_set = cycle_num_set(num_set)
    datur = {k: num_set[k][0] for k in num_set.keys()}
    return pd.DataFrame({"ops": list(datur.keys()), "val": list(datur.values())})


if __name__ == "__main__":
    df = cycle8833()
    print (df.query("val >= 23.99 & val <= 24.01")["ops"].values[0])


