import multiprocessing
import concurrent.futures

name = 'bounded_pool_executor'
class _BoundedPoolExecutor:

    semaphore = None
    
    def acquire(self):
        self.semaphore.acquire()
        
    def release(self, fn):
        self.semaphore.release()
        
    def submit(self, fn, *args, **kwargs):
        self.acquire()
        future = super().submit(fn, *args, **kwargs)
        future.add_done_callback(self.release)
        
        return future
    

class BoundedProcessPoolExecutor(_BoundedPoolExecutor, concurrent.futures.ProcessPoolExecutor):
    def __init__(self, max_workers=None,mp_context=None, initializer=None, initargs=()):
        super().__init__(max_workers,mp_context,initializer,initargs)
        self.semaphore = multiprocessing.BoundedSemaphore(max_workers)
                                                                                                        
