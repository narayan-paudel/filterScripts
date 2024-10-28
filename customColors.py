import numpy as np

# def qualitative_colors(n):
#     if n < 1:
#         raise ValueError('Minimum number of qualitative colors is 1.')
#     elif n > 12:
#         raise ValueError('Maximum number of qualitative colors is 12.')
#     cols = ['#4477AA', '#332288', '#6699CC', '#88CCEE', '#44AA99', '#117733',
#             '#999933', '#DDCC77', '#661100', '#CC6677', '#AA4466', '#882255',
#             '#AA4499']
#     indices = [[0],
#                [0, 9],
#                [0, 7, 9],
#                [0, 5, 7, 9],
#                [1, 3, 5, 7, 9],
#                [1, 3, 5, 7, 9, 12],
#                [1, 3, 4, 5, 7, 9, 12],
#                [1, 3, 4, 5, 6, 7, 9, 12],
#                [1, 3, 4, 5, 6, 7, 9, 11, 12],
#                [1, 3, 4, 5, 6, 7, 8, 9, 11, 12],
#                [1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12],
#                [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]]
#     return [cols[ix] for ix in indices[n - 1]]


def qualitative_colors(n):
    if n < 1:
        raise ValueError('Minimum number of qualitative colors is 1.')
    elif n > 12:
        raise ValueError('Maximum number of qualitative colors is 12.')
    cols = ['#4477AA', '#332288', '#6699CC', '#88CCEE', '#44AA99', '#117733',
            '#999933', '#DDCC77',
            # '#e2514a',
            # '#fca55d',
            # '#fee999',
             # '#edf8a3',
              '#a2d9a4',
              # '#47a0b3', 
            '#661100','#882255', '#AA4466','#AA4499',  '#CC6677']

    indices = [[0],
               [0, 9],
               [0, 7, 9],
               [0, 5, 7, 9],
               [1, 3, 5, 7, 9],
               [1, 3, 5, 7, 9, 12],
               [1, 3, 4, 5, 7, 9, 12],
               [1, 3, 4, 5, 6, 7, 9, 12],
               [1, 3, 4, 5, 6, 7, 9, 11, 12],
               [1, 3, 4, 5, 6, 7, 8, 9, 11, 12],
               [1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12],
               [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]]
    return [cols[ix] for ix in indices[n - 1]]