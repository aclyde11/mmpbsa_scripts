import os
import sys
from contextlib import contextmanager


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def parse_final(result_file):
    with open(result_file, 'r') as f:
        for line in f:
            if "DELTA TOTAL" in line:
                return float(line.split()[2]), float(line.split()[3])
        return None


def make_message_writer(verbose_, class_name_, enter_message_=False):
    class MessageWriter(object):
        class_name = class_name_

        def __init__(self, method_name=None, verbose=None, enter_message=None):
            self.verbose = verbose or verbose_
            self.method_name = method_name or "unknown_method"
            self.enter_message = enter_message or enter_message_

        def log(self, *args, **kwargs):
            if self.verbose:
                print(f"INFO [{self.class_name}:{self.method_name}]", *args, **kwargs)

        def error(self, *args, exit_all=False, **kwargs):
            print(f"{bcolors.WARNING}ERROR [{self.class_name}:{self.method_name}]{bcolors.ENDC}", *args, **kwargs,
                  file=sys.stderr)
            if exit_all:
                exit()

        def failure(self, *args, exit_all=False, **kwargs):
            print(f"{bcolors.WARNING}FAILURE [{self.class_name}:{self.method_name}]{bcolors.ENDC}", *args, **kwargs,
                  file=sys.stderr)
            if exit_all:
                exit()

        @classmethod
        def static_failure(cls, method_name, *args, exit_all=False, **kwargs):
            print(f"{bcolors.WARNING}ERROR [{cls.class_name}:{method_name}]{bcolors.ENDC}", *args, **kwargs,
                  file=sys.stderr)
            if exit_all:
                exit()

        @classmethod
        def static_log(cls, method_name, verbose, *args, **kwargs):
            if verbose:
                print(f"INFO [{cls.class_name}:{method_name}]", *args, **kwargs)

        def __enter__(self):
            if self.enter_message:
                self.log("Entering")
            return self

        def __exit__(self, *args, **kwargs):
            if self.enter_message:
                self.log("Exiting")

    return MessageWriter


@contextmanager
def working_directory(directory):
    owd = os.getcwd()
    try:
        os.chdir(directory)
        yield directory
    finally:
        os.chdir(owd)
