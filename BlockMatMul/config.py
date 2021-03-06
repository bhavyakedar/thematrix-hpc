'''
problem_size = {'trapezoidal': [10**x for x in range(1, 9)],
                'vector': [10**x for x in range(1, 9)],
                'matrix_multiplication': [2**x for x in range(4, 11)],
                'pi_using_series': [10**x for x in range(1, 9)],
                'divide_and_conquer': [10**x for x in range(1, 9)] + [int(0.5*10**x) for x in range(5,9)],
                'image_warping': [2**x for x in range(3, 11)],
                'median_filtering': [2**x for x in range(3, 11)],
                'monte_carlo': [10**x for x in range(1, 9)] + [int(0.5*10**x) for x in range(5,9)] ,
                'prefix_sum': [10**x for x in range(1, 9)],
                'reduction': [10**x for x in range(1, 9)] + [int(0.5*10**x) for x in range(5,9)] ,
                'filter': [10**x for x in range(1, 9)]}
'''
prob = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096]
## warning, don't remove 0, main.py doesn't rectify
procs = [0]+[2, 4, 8, 16]
