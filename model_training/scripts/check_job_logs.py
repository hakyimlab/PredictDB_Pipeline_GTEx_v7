#! /usr/bin/env python

import os
import re
import sys

def relavent_logfiles(p):
    """Generator yielding relevant output log file name"""
    for file in os.listdir("../joblogs/"):
        if re.match(p, file):
            yield file


def check_job_logs(p):
    """Checks all log files in a jobs log directory to make sure they
    turn out ok.
    Makes sure the model script processed all genes it set out to"""
    
    nfiles = 0
    nprobs = 0
    for log_file in relavent_logfiles(p):
        nfiles += 1
        if os.stat('../joblogs/' + log_file).st_size == 0:
            # File is empty
            nprobs += 1
            print("{} is empty.".format(log_file))
        with open('../joblogs/' + log_file, 'r') as lf:
            # Go to last line.
            for line in lf:
                pass
            nums = line.strip().split(' / ')
            assert len(nums) == 2
            if nums[0] != nums[1] or nums[1] == '0':
                nprobs += 1
                print("Problem with {}".format(log_file))
                print(nums)
                #if nprobs > 20:
                #    print("Too many problems")
                #    print("{} files check so far".format(nfiles))
                #    return
    print("{} files found".format(nfiles))
    print("{} check out".format(nfiles - nprobs))

if __name__ == "__main__":
    p = 'gtex_training_.*\.out'
    check_job_logs(p)
