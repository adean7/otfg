

def abort(error_msg):
    print(error_msg)
    exit(1)

def run_time_out(currentParams):
    print('\nTime of run: {:>9.3f} s.'.format(currentParams.get_run_time()))
