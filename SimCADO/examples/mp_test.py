import multiprocessing as mp
import time
import numpy as np

def wait(i, j):
    print(i, j)
    time.sleep(2)
    return np.ones((2,2))*i

print(__name__)    
if __name__ == "__main__":    
    
    pool = mp.Pool(processes=7)
    ps = [pool.apply_async(wait, (i, 2)) for i in range(22)]
    results = [p.get() for p in ps]
    print(results)
        