import multiprocessing as mpr


def dangerwrap(f, *args, **kwargs):
    event = mpr.Event()
    q = mpr.Queue()

    def signalling_f():
        q.put(f(*args, **kwargs))
        event.set()

    f_process = mpr.Process(target=signalling_f)
    f_process.start()
    try:
        event.wait()

    except KeyboardInterrupt:
        f_process.terminate()
        f_process.join()
        raise KeyboardInterrupt()

    return q.get()
