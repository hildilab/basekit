import datetime


# adapted from http://stackoverflow.com/questions/156330/get-timer-ticks-in-python

def now():
    return datetime.datetime.now()

# Prints one of the following formats*:
# 1.58 days
# 2.98 hours
# 9.28 minutes # Not actually added yet, oops.
# 5.60 seconds
# 790 milliseconds
# *Except I prefer abbreviated formats, so I print d,h,m,s, or ms. 
def format_delta(start,end):

    # Time in microseconds
    one_day = 86400000000
    one_hour = 3600000000
    one_second = 1000000
    one_millisecond = 1000

    delta = end - start

    build_time_us = delta.microseconds + delta.seconds * one_second + delta.days * one_day

    days = 0
    while build_time_us > one_day:
        build_time_us -= one_day
        days += 1

    if days > 0:
        time_str = "%.2fd" % ( days + build_time_us / float(one_day) )
    else:
        hours = 0
        while build_time_us > one_hour:
            build_time_us -= one_hour
            hours += 1
        if hours > 0:
            time_str = "%.2fh" % ( hours + build_time_us / float(one_hour) )
        else:
            seconds = 0
            while build_time_us > one_second:
                build_time_us -= one_second
                seconds += 1
            if seconds > 0:
                time_str = "%.2fs" % ( seconds + build_time_us / float(one_second) )
            else:
                ms = 0
                while build_time_us > one_millisecond:
                    build_time_us -= one_millisecond
                    ms += 1
                time_str = "%.2fms" % ( ms + build_time_us / float(one_millisecond) )
    return time_str




class Timer:
    def __init__(self, name="timer"):
        self.name = name
    def __enter__(self):
        self.begin = now()
        print "start '%s'" % self.name
    def __exit__(self, type, value, traceback):
        print "end '%s' in %s" % ( self.name, format_delta( self.begin, now() ) )

