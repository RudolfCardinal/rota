#!/usr/bin/python3.4
# -*- encoding: utf8 -*-

"""

Medical rota design/checking tool.

    See README.md

Copyright/licensing

    Copyright (C) 2015-2015 Rudolf Cardinal (rudolf@pobox.com).
    Department of Psychiatry, University of Cambridge.

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

        http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.

Sources

    See footnotes.

Single-file format

    Try to keep this software in a single file, to aid ad-hoc distribution.

Version history (see VERSION below)

    1.0 (2015-08-08 to 2015-08-12)
    - First version

    1.1 (2015-08-13)
    - CPFT Draft 2 rota spread out to be more realistic.
    - CPFT Draft 3 added.

    1.2 (2015-08-14)
    - Prototype display added.
    - Parallel processing for banding calculations, role coverage.
    - CPFT Draft 4 added.
    - Bugfix to get_work_interval(ignore_bank_holidays=True) -- wasn't
      ignoring bank holidays properly for prospective cover calculations.
    - CPFT North SHO components of drafts 2/3/4 changed to avoid Band 2 as a
      result (one extra day off).
    - SpR components of CPFT drafts 3/4 changed similarly. ***
"""

# =============================================================================
# Version
# =============================================================================

VERSION = 1.2

# =============================================================================
# Imports
# =============================================================================

import logging
LOG_FORMAT = '%(asctime)s.%(msecs)03d:%(levelname)s:%(name)s:%(message)s'
LOG_DATEFMT = '%Y-%m-%d %H:%M:%S'
logging.basicConfig(format=LOG_FORMAT, datefmt=LOG_DATEFMT)
logger = logging.getLogger("rota")
logger.setLevel(logging.DEBUG)

import argparse
import cgi
from collections import OrderedDict
import concurrent.futures
import datetime
import string
import sys


# =============================================================================
# AttrDict
# =============================================================================

class AttrDict(dict):
    # http://stackoverflow.com/questions/4984647
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self


# =============================================================================
# Constants
# =============================================================================

BANK_HOLIDAYS = [datetime.datetime.strptime(x, "%Y-%m-%d").date() for x in [
    # https://www.gov.uk/bank-holidays
    # All bank holiday dates vary, even the date-based ones; e.g. if Christmas
    # Day is a Sunday, then the Christmas Day substitute bank holiday is Tue 27
    # Dec, after the Boxing Day Monday bank holiday.

    # 2014
    "2014-01-01",  # New Year's Day
    "2014-04-18",  # Good Friday
    "2014-04-21",  # Easter Monday
    "2014-05-05",  # Early May Bank Holiday
    "2014-05-26",  # Spring Bank Holiday
    "2014-08-25",  # Summer Bank Holiday
    "2014-12-25",  # Christmas Day
    "2014-12-26",  # Boxing Day
    # 2015
    "2015-01-01",  # New Year's Day
    "2015-04-03",  # Good Friday
    "2015-04-06",  # Easter Monday
    "2015-05-04",  # Early May Bank Holiday
    "2015-05-25",  # Spring Bank Holiday
    "2015-08-31",  # Summer Bank Holiday
    "2015-12-25",  # Christmas Day
    "2015-12-28",  # Boxing Day (substitute)
    # 2016
    "2016-01-01",  # New Year's Day
    "2016-03-25",  # Good Friday
    "2016-03-28",  # Easter Monday
    "2016-05-02",  # Early May Bank Holiday
    "2016-05-30",  # Spring Bank Holiday
    "2016-08-29",  # Summer Bank Holiday
    "2016-12-26",  # Boxing Day
    "2016-12-27",  # Christmas Day (substitute)

    # Don't forget to add more in years to come.
]]

SUPPLEMENT = {
    # [BMA_1]
    "3": 1.0,
    "2A": 0.8,
    "2B": 0.5,
    "1A": 0.5,
    "1B": 0.4,
    "1C": 0.2,
    None: 0,
}

SHIFT_TYPES = AttrDict({
    "FULL": "Full shift",
    "PARTIAL": "Partial shift",
    "PARTIAL24": "24-hour partial shift",
    "ONCALL": "On-call",
    "NONE": None,
})

SECONDS_PER_MINUTE = 60
MINUTES_PER_HOUR = 60
HOURS_PER_DAY = 24
DAYS_PER_WEEK = 7

SECONDS_PER_HOUR = SECONDS_PER_MINUTE * MINUTES_PER_HOUR
SECONDS_PER_DAY = SECONDS_PER_HOUR * HOURS_PER_DAY
SECONDS_PER_WEEK = SECONDS_PER_DAY * DAYS_PER_WEEK

ARBITRARY_DATE = datetime.date(2015, 1, 1)
ARBITRARY_MONDAY_NEAR_BH = datetime.date(2015, 8, 17)
ARBITRARY_MONDAY_FAR_FROM_BH = datetime.date(2015, 6, 1)

DP = 2

NORMAL_DAY_START_H = 7
NORMAL_DAY_END_H = 19

COLOURS = AttrDict({
    "NWD": (100, 150, 100),
    "LATE_F": (255, 0, 0),
    "LATE_A": (255, 127, 0),
    "NIGHT_FA": (150, 150, 255),
    "LATE_C": (255, 255, 51),
    "LATE_P": (102, 255, 255),
    "NIGHT_CP": (153, 51, 255),
    "SPR_LATE": (251, 154, 153),
    "N_SPR_LATE": (251, 154, 153),
    "S_SPR_LATE": (0, 150, 150),
    "SPR_NIGHT": (255, 0, 255),
    "S_SPR_NIGHT": (255, 0, 255),
    "N_SPR_NIGHT": (0, 255, 0),
    "OFF": (255, 255, 255),
})


# =============================================================================
# Standalone functions
# =============================================================================

def is_bank_holiday(date):
    """Is the specified date (a datetime.date object) a UK bank holiday?
    Uses the BANK_HOLIDAYS list."""
    return date in BANK_HOLIDAYS


def is_weekend(date):
    """Is the specified date (a datetime.date object) a weekend?"""
    return date.weekday() in [5, 6]


def is_saturday(date):
    """Is the specified date (a datetime.date object) a Saturday?"""
    return date.weekday() == 5


def is_sunday(date):
    """Is the specified date (a datetime.date object) a Sunday?"""
    return date.weekday() == 6


def is_normal_working_day(date):
    """Is the specified date (a datetime.date object) a normal working day,
    i.e. not a weekend or a bank holiday?"""
    return not(is_weekend(date) or is_bank_holiday(date))


def css_compatible(name):
    """Is the name suitable for use as a CSS class name?
    This is rough and ready!"""
    for c in name:
        if not c.isalnum() and c != '_':
            return False
    return True


def yesno(x):
    """Maps truthy -> 'Yes' and falsy -> 'No'."""
    return "Yes" if x else "No"


def number_to_dp(number, dp, default="", en_dash_for_minus=True):
    """Format number to dp decimal places, optionally using a UTF-8 en dash
    for minus signs."""
    if number is None:
        return default
    s = "{:.{precision}f}".format(number, precision=dp)
    if en_dash_for_minus:
        s = s.replace("-", "–")  # hyphen becomes en dash for minus sign
    return s


def convert_duration(duration, units):
    """Convert a datetime.timedelta object -- a duration -- into other
    units. Possible units:
        s, sec, seconds
        m, min, minutes
        h, hr, hours
        d, days
        w, weeks
    """
    if duration is None:
        return None
    s = duration.total_seconds()
    if units in ['s', 'sec', 'seconds']:
        return s
    if units in ['m', 'min', 'minutes']:
        return s / SECONDS_PER_MINUTE
    if units in ['h', 'hr', 'hours']:
        return s / SECONDS_PER_HOUR
    if units in ['d', 'days']:
        return s / SECONDS_PER_DAY
    if units in ['w', 'weeks']:
        return s / SECONDS_PER_WEEK
    raise ValueError("Unknown units: {}".format(units))


def webify(v, preserve_newlines=True):
    """Converts a value into an HTML-safe str. Python 3 version.

    Converts value v to a string; escapes it to be safe in HTML
    format (escaping ampersands, replacing newlines with <br>, etc.).
    Returns str/unicode, depending on input. Returns "" for blank input.
    """
    nl = "<br>" if preserve_newlines else " "
    if v is None:
        return ""
    return cgi.escape(v).replace("\n", nl).replace("\\n", nl)


def formatdt(date):
    """Formats a date to ISO-8601 basic format, to minute accuracy with no
    timezone."""
    return date.strftime("%Y-%m-%dT%H:%M")


# =============================================================================
# Time intervals and lists of intervals
# =============================================================================

# -----------------------------------------------------------------------------
# Interval
# -----------------------------------------------------------------------------

class Interval(object):
    """Object representing a time interval, with start and end objects that are
    normally datetime.datetime objects (though with care, a subset of some
    methods are possible with datetime.date objects; caveat emptor, and some
    methods will crash).

    Does not handle open-ended intervals (−∞, +∞).

    There's probably an existing class for this...
    """
    def __init__(self, start, end):
        """Creates the interval."""
        if start is None or end is None:
            raise TypeError("Invalid interval creation")
        if start > end:
            (start, end) = (end, start)
        self.start = start
        self.end = end

    def __repr__(self):
        return "Interval(start={}, end={})".format(
            repr(self.start), repr(self.end))

    def __str__(self):
        return "{} − {}".format(formatdt(self.start), formatdt(self.end))

    def __add__(self, value):
        """Adds a constant (datetime.timedelta object) to the interval's start
        and end."""
        return Interval(self.start + value, self.end + value)

    def __lt__(self, other):
        """Allows sorting (on start time)."""
        return self.start < other.start

    def copy(self):
        """Makes a copy of the interval."""
        return Interval(self.start, self.end)

    def overlaps(self, other):
        """
        Does this interval overlap the other?
        Overlap:
                   S--------S     S---S            S---S
                     O---O          O---O        O---O
        Simpler method of testing is for non-overlap!
                   S---S              S---S
                       O---O      O---O
        """
        return not(self.end <= other.start or self.start >= other.end)

    def contiguous(self, other):
        """Does this interval overlap or touch the other?"""
        return not(self.end < other.start or self.start > other.end)

    def contains(self, time, inclusive=True):
        """Does the interval contain a momentary time?"""
        if inclusive:
            return (self.start <= time and time <= self.end)
        else:
            return (self.start < time and time < self.end)

    def union(self, other):
        """Returns an interval spanning the extent of this and the other."""
        return Interval(
            min(self.start, other.start),
            max(self.end, other.end)
        )

    def intersection(self, other):
        """Returns an interval representing the intersection of this and the
        other, or None if they don't overlap."""
        if not self.contiguous(other):
            return None
        return Interval(
            max(self.start, other.start),
            min(self.end, other.end)
        )

    def cut(self, times):
        """Returns a list of intervals produced by using times (a list of
        datetime.datetime objects, or a single such object) as a set of knives
        to slice this interval."""
        if not isinstance(times, list):
            # Single time
            time = times
            if not self.contains(time):
                return []
            return [
                Interval(self.start, time),
                Interval(time, self.end)
            ]
        else:
            # Multiple times
            times = [t for t in times if self.contains(t)]  # discard others
            times.sort()
            times = [self.start] + times + [self.end]
            intervals = []
            for i in range(len(times) - 1):
                intervals.append(Interval(times[i], times[i + 1]))
            return intervals

    def duration(self):
        """Returns a datetime.timedelta object representing the duration of
        this interval."""
        return self.end - self.start

    def duration_in(self, units):
        """Returns the duration of this interval in the specified units, as
        per convert_duration()."""
        return convert_duration(self.duration(), units)

    @staticmethod
    def wholeday(date):
        """Returns an Interval covering the date given (midnight at the start
        of that day to midnight at the start of the next day)."""
        start = datetime.datetime.combine(date, datetime.time())
        return Interval(
            start,
            start + datetime.timedelta(days=1)
        )

    @staticmethod
    def daytime(date, daybreak=datetime.time(NORMAL_DAY_START_H),
                nightfall=datetime.time(NORMAL_DAY_END_H)):
        """Returns an Interval representing daytime on the date given."""
        return Interval(
            datetime.datetime.combine(date, daybreak),
            datetime.datetime.combine(date, nightfall),
        )

    @staticmethod
    def dayspan(startdate, enddate, include_end=True):
        """Returns an Interval representing the date range given, from midnight
        at the start of the first day to midnight at the end of the last (i.e.
        at the start of the next day after the last), or if include_end is
        False, 24h before that."""
        if enddate < startdate:
            return None
        if enddate == startdate and include_end:
            return None
        start_dt = datetime.datetime.combine(startdate, datetime.time())
        end_dt = datetime.datetime.combine(enddate, datetime.time())
        if include_end:
            end_dt += datetime.timedelta(days=1)
        return Interval(start_dt, end_dt)

    def component_on_date(self, date):
        """Returns the part of this interval that falls on the date given,
        or None if the interval doesn't have any part during that date."""
        return self.intersection(Interval.wholeday(date))

    def day_night_duration(self, daybreak=datetime.time(NORMAL_DAY_START_H),
                           nightfall=datetime.time(NORMAL_DAY_END_H)):
        """Returns a (day, night) tuple of datetime.timedelta objects giving
        the duration of this interval that falls into day and night
        respectively."""
        daytotal = datetime.timedelta()
        nighttotal = datetime.timedelta()
        startdate = self.start.date()
        enddate = self.end.date()
        ndays = (enddate - startdate).days + 1
        for i in range(ndays):
            date = startdate + datetime.timedelta(days=i)
            component = self.component_on_date(date)
            # ... an interval on a single day
            day = Interval.daytime(date, daybreak, nightfall)
            daypart = component.intersection(day)
            if daypart is not None:
                daytotal += daypart.duration()
                nighttotal += component.duration() - daypart.duration()
            else:
                nighttotal += component.duration()
        return (daytotal, nighttotal)

    def duration_outside_nwh(self,
                             starttime=datetime.time(NORMAL_DAY_START_H),
                             endtime=datetime.time(NORMAL_DAY_END_H)):
        """Returns a duration (a datetime.timedelta object) representing the
        number of hours outside normal working hours."""
        ooh = datetime.timedelta()  # ooh = out of (normal) hours
        startdate = self.start.date()
        enddate = self.end.date()
        ndays = (enddate - startdate).days + 1
        for i in range(ndays):
            date = startdate + datetime.timedelta(days=i)
            component = self.component_on_date(date)
            # ... an interval on a single day
            if not is_normal_working_day(date):
                ooh += component.duration()  # all is out-of-normal-hours
            else:
                normalday = Interval.daytime(date, starttime, endtime)
                normalpart = component.intersection(normalday)
                if normalpart is not None:
                    ooh += component.duration() - normalpart.duration()
                else:
                    ooh += component.duration()
        return ooh

    def n_weekends(self):
        """Returns the number of weekends that this interval covers."""
        startdate = self.start.date()
        enddate = self.end.date()
        ndays = (enddate - startdate).days + 1
        in_weekend = False
        n_weekends = 0
        for i in range(ndays):
            date = startdate + datetime.timedelta(days=i)
            if not in_weekend and is_weekend(date):
                in_weekend = True
                n_weekends += 1
            elif in_weekend and not is_weekend(date):
                in_weekend = False
        return n_weekends

    def saturdays_of_weekends(self):
        """Returns the dates of all Saturdays that are part of weekends that
        this interval covers (representing a unique identifier for that
        weekend). The Saturday itself isn't necessarily the part of the weekend
        that the interval covers!"""
        startdate = self.start.date()
        enddate = self.end.date()
        ndays = (enddate - startdate).days + 1
        saturdays = set()
        for i in range(ndays):
            date = startdate + datetime.timedelta(days=i)
            if is_saturday(date):
                saturdays.add(date)
            elif is_sunday(date):
                saturdays.add(date - datetime.timedelta(days=1))
        return saturdays


# -----------------------------------------------------------------------------
# IntervalList
# -----------------------------------------------------------------------------

class IntervalList(object):
    """Object representing a list of Intervals.
    Maintains an internally sorted state (by interval start time)."""

    def __init__(self, intervals=None, no_overlap=True,
                 no_contiguous=True):
        """Creates the IntervalList."""
        # DO NOT USE intervals=[]; that's the route to a mutable default and
        # a huge amount of confusion as separate objects appear
        # non-independent.
        self.intervals = [] if intervals is None else list(intervals)
        self.no_overlap = no_overlap
        self.no_contiguous = no_contiguous
        for i in self.intervals:
            if not isinstance(i, Interval):
                raise TypeError(
                    "IntervalList creation failed: contents are not all "
                    "Interval: {}".format(repr(self.intervals)))
        self._tidy()

    def __repr__(self):
        return (
            "IntervalList(intervals={}, no_overlap={}, "
            "no_contiguous={})".format(
                repr(self.intervals),
                self.no_overlap,
                self.no_contiguous))

    def copy(self, no_overlap=None, no_contiguous=None):
        """Makes a copy of the IntervalList. The overlap/contiguous parameters
        can be changed."""
        if no_overlap is None:
            no_overlap = self.no_overlap
        if no_contiguous is None:
            no_contiguous = self.no_contiguous
        return IntervalList(self.intervals, no_overlap=no_overlap,
                            no_contiguous=no_contiguous)

    def add(self, interval):
        """Adds an interval to the list. If self.no_overlap is True, as is the
        default, it will merge any overlapping intervals thus created."""
        if interval is None:
            return
        if not isinstance(interval, Interval):
            raise TypeError(
                "Attempt to insert non-Interval into IntervalList")
        self.intervals.append(interval)
        self._tidy()

    def _tidy(self):
        """Removes overlaps, etc., and sorts."""
        if self.no_overlap:
            self.remove_overlap(self.no_contiguous)  # will sort
        else:
            self.sort()

    def sort(self):
        """Sorts (in place) by interval start."""
        self.intervals.sort()

    def list(self):
        """Returns the contained list."""
        return self.intervals

    def _remove_overlap_sub(self, also_remove_contiguous):
        # Returns True if overlap removed; False otherwise
        for i in range(len(self.intervals)):
            for j in range(i + 1, len(self.intervals)):
                first = self.intervals[i]
                second = self.intervals[j]
                if also_remove_contiguous:
                    test = first.contiguous(second)
                else:
                    test = first.overlaps(second)
                if test:
                    newint = first.union(second)
                    self.intervals.pop(j)
                    self.intervals.pop(i)
                    self.intervals.append(newint)
                    return True
        return False

    def remove_overlap(self, also_remove_contiguous=False):
        """Merges any overlapping intervals."""
        overlap = True
        while overlap:
            overlap = self._remove_overlap_sub(also_remove_contiguous)
        self.sort()

    def _any_overlap_or_contiguous(self, test_overlap):
        """Do any of the intervals overlap?"""
        for i in range(len(self.intervals)):
            for j in range(i + 1, len(self.intervals)):
                first = self.intervals[i]
                second = self.intervals[j]
                if test_overlap:
                    test = first.overlaps(second)
                else:
                    test = first.contiguous(second)
                if test:
                    return True
        return False

    def any_overlap(self):
        """Do any of the intervals overlap?"""
        return self._any_overlap_or_contiguous(test_overlap=True)

    def any_contiguous(self):
        """Are any of the intervals contiguous?"""
        return self._any_overlap_or_contiguous(test_overlap=False)

    def get_overlaps(self):
        """Returns an IntervalList containing intervals representing periods of
        overlap between intervals in this one."""
        overlaps = IntervalList()
        for i in range(len(self.intervals)):
            for j in range(i + 1, len(self.intervals)):
                first = self.intervals[i]
                second = self.intervals[j]
                ol = first.intersection(second)
                if ol is not None:
                    overlaps.add(ol)
        return overlaps

    def total_duration(self):
        """Returns a datetime.timedelta object with the total sum of durations.
        If there is overlap, time will be double-counted, so beware!"""
        total = datetime.timedelta()
        for interval in self.intervals:
            total += interval.duration()
        return total

    def n_weekends(self):
        """Returns the number of weekends that the intervals collectively
        touch."""
        saturdays = set()
        for interval in self.intervals:
            saturdays.update(interval.saturdays_of_weekends())
        return len(saturdays)

    def duration_outside_nwh(self,
                             starttime=datetime.time(NORMAL_DAY_START_H),
                             endtime=datetime.time(NORMAL_DAY_END_H)):
        """Returns the total duration outside normal working hours, i.e.
        evenings/nights, weekends (and Bank Holidays)."""
        total = datetime.timedelta()
        for interval in self.intervals:
            total += interval.duration_outside_nwh(starttime, endtime)
        return total

    def durations(self):
        """Returns a list of datetime.timedelta objects representing the
        duration of intervals."""
        return [x.duration() for x in self.intervals]

    def longest_duration(self):
        """Returns the duration of the longest interval, or None if none."""
        if not self.intervals:
            return None
        return max(self.durations())

    def longest_interval(self):
        """Returns the longest interval, or None if none."""
        longest_duration = self.longest_duration()
        for i in self.intervals:
            if i.duration() == longest_duration:
                return i
        return None

    def shortest_duration(self):
        """Returns the duration of the longest interval, or None if none."""
        if not self.intervals:
            return None
        return min(self.durations())

    def shortest_interval(self):
        """Returns the shortest interval, or None if none."""
        shortest_duration = self.shortest_duration()
        for i in self.intervals:
            if i.duration() == shortest_duration:
                return i
        return None

    def gaps(self):
        """Returns all the gaps between intervals, as an IntervalList."""
        if len(self.intervals) < 2:
            return IntervalList(None)
        gaps = []
        for i in range(len(self.intervals) - 1):
            gap = Interval(
                self.intervals[i].end,
                self.intervals[i + 1].start
            )
            gaps.append(gap)
        return IntervalList(gaps)

    def shortest_gap(self):
        """Find the shortest gap between intervals."""
        gaps = self.gaps()
        return gaps.shortest_interval()

    def shortest_gap_duration(self):
        """Find the duration of the shortest gap between intervals."""
        gaps = self.gaps()
        return gaps.shortest_duration()

    def start_date(self):
        if not self.intervals:
            return None
        return self.intervals[0].start.date()

    def end_date(self):
        if not self.intervals:
            return None
        return self.intervals[-1].end.date()

    def max_continuous_days(self):
        """
        The length of the longest sequence of days in which all days include
        an interval. Returns a tuple: (length, interval with start and end
        date of longest span).
        """
        if len(self.intervals) == 0:
            return None
        startdate = self.start_date()
        enddate = self.end_date()
        seq = ''
        ndays = (enddate - startdate).days + 1
        in_span = False
        for i in range(ndays):
            date = startdate + datetime.timedelta(days=i)
            wholeday = Interval.wholeday(date)
            if any([x.overlaps(wholeday) for x in self.intervals]):
                seq += '+'
            else:
                seq += ' '
        longest = max(seq.split(), key=len)
        longest_len = len(longest)
        longest_idx = seq.index(longest)
        longest_interval = Interval.dayspan(
            startdate + datetime.timedelta(days=longest_idx),
            startdate + datetime.timedelta(days=longest_idx + longest_len)
        )
        return (longest_len, longest_interval)

    def subset(self, interval, flexibility=2):
        """
        Returns an IntervalList that's a subset of this one, only containing
        intervals that meet the "interval" parameter criterion.

        flexibility == 0: permits only wholly contained intervals:

                        I----------------I
                N---N  N---N   Y---Y   N---N   N---N
                    N---N                N---N

        flexibility == 1: permits overlapping intervals as well:

                        I----------------I
                N---N  Y---Y   Y---Y   Y---Y   N---N
                    N---N                N---N

        flexibility == 2: permits adjoing intervals as well:

                        I----------------I
                N---N  Y---Y   Y---Y   Y---Y   N---N
                    Y---Y                Y---Y
        """
        if flexibility not in [0, 1, 2]:
            raise ValueError("subset: bad flexibility value")
        permitted = []
        for i in self.intervals:
            if flexibility == 0:
                ok = i.start > interval.start and i.end < interval.end
            elif flexibility == 1:
                ok = i.end > interval.start and i.start < interval.end
            else:
                ok = i.end >= interval.start and i.start <= interval.end
            if ok:
                permitted.append(i)
        return IntervalList(permitted)

    def gap_subset(self, interval, flexibility=2):
        """
        Returns an IntervalList that's a subset of this one, only containing
        *gaps* between intervals that meet the interval criterion.
        """
        return self.gaps().subset(interval, flexibility)

    def first_interval_starting(self, start):
        """Returns an interval starting with the start parameter, or None."""
        for i in self.intervals:
            if i.start == start:
                return i
        return None

    def first_interval_ending(self, end):
        """Returns an interval ending with the end parameter, or None."""
        for i in self.intervals:
            if i.end == end:
                return i
        return None

    def _sufficient_gaps(self, startdate, enddate, requiredgaps,
                         flexibility):
        """
        Are there sufficient gaps (specified by requiredgaps) in the date
        range specified? This is a worker function for sufficient_gaps.
        """
        requiredgaps = list(requiredgaps)  # make a copy
        interval = Interval.dayspan(startdate, enddate, include_end=True)
        # logger.debug(">>> _sufficient_gaps")
        gaps = self.gap_subset(interval, flexibility)
        gapdurations = gaps.durations()
        gaplist = gaps.list()
        gapdurations.sort(reverse=True)  # longest gap first
        requiredgaps.sort(reverse=True)  # longest gap first
        # logger.debug("... gaps = {}".format(gaps))
        # logger.debug("... gapdurations = {}".format(gapdurations))
        # logger.debug("... requiredgaps = {}".format(requiredgaps))
        while requiredgaps:
            # logger.debug("... processing gap")
            if not gapdurations:
                # logger.debug("<<< no gaps left")
                return False, None  # ***
            if gapdurations[0] < requiredgaps[0]:
                # logger.debug("<<< longest gap is too short")
                return False, self.first_interval_ending(gaplist[0].start)
            gapdurations.pop(0)
            requiredgaps.pop(0)
            gaplist.pop(0)
            # ... keeps gaplist and gapdurations mapped to each other
        # logger.debug("<<< success")
        return True, None

    def sufficient_gaps(self, every_n_days, requiredgaps, flexibility=2):
        """
        Are gaps present sufficiently often?
        For example:
            every_n_days=21
            requiredgaps=[
                datetime.timedelta(hours=62),
                datetime.timedelta(hours=48),
            ]
        ... means "is there at least one 62-hour gap and one (separate) 48-hour
        gap in every possible 21-day sequence within the IntervalList?

        If flexibility == 0:
            ... gaps must be WHOLLY WITHIN the interval.
        If flexibility == 1:
            ... gaps may OVERLAP the edges of the interval.
        If flexibility == 2:
            ... gaps may ABUT the edges of the interval.

        Returns (True, None) or (False, first_failure_interval).
        """
        if len(self.intervals) < 2:
            return False, None
        startdate = self.start_date()
        enddate = self.end_date()
        ndays = (enddate - startdate).days + 1
        if ndays <= every_n_days:
            # Our interval is too short, or just right
            return self._sufficient_gaps(startdate, enddate, requiredgaps,
                                         flexibility)
        for i in range(ndays - every_n_days):
            j = i + every_n_days
            a = startdate + datetime.timedelta(days=i)
            b = startdate + datetime.timedelta(days=j)
            (sufficient, ffi) = self._sufficient_gaps(a, b, requiredgaps,
                                                      flexibility)
            if not sufficient:
                return False, ffi
        return True, None


if False:
    TESTCODE = """
a = datetime.datetime(2015, 1, 1)
b = datetime.datetime(2015, 1, 6)
i = Interval(a, b)
j = i + datetime.timedelta(hours=3)
i.cut(datetime.datetime(2015, 1, 3))
    """


# =============================================================================
# Shift definition
# =============================================================================

class Shift(object):
    """Object representing a shift pattern."""

    def __init__(self, name, abbreviation, start_time=None, duration_h=None,
                 work=True, roles=None, nwd_only=False, resident=False,
                 shift_type=SHIFT_TYPES.NONE, rgb=(255, 255, 255)):
        """Initializes the Shift."""
        self.name = name
        self.abbreviation = abbreviation
        self.start_time = start_time
        self.duration_h = duration_h
        self.work = work
        self.roles = [] if roles is None else list(roles)
        self.nwd_only = nwd_only
        self.resident = resident
        self.shift_type = shift_type
        self.rgb = rgb
        # Validation
        if not css_compatible(self.name):
            raise ValueError(
                "Shift name illegal as a CSS name; change it: {}".format(
                    repr(self)))
        if not isinstance(start_time, datetime.time):
            raise TypeError(
                "Shift time not a datetime.time: {}".format(repr(self)))
        if shift_type not in SHIFT_TYPES.values():
            raise ValueError(
                "Shift type ({}) must be one of [{}]".format(
                    shift_type, ", ".join(SHIFT_TYPES.values())))

    def __repr__(self):
        return (
            "Shift(name={}, abbreviation={}, start_time={}, duration_h={}, "
            "work={}, roles={}, nwd_only={}, resident={}, shift_type={}, "
            "rgb={})".format(
                repr(self.name),
                repr(self.abbreviation),
                repr(self.start_time),
                self.duration_h,
                self.work,
                self.roles if self.roles else None,
                self.nwd_only,
                self.resident,
                repr(self.shift_type),
                self.rgb,
            )
        )

    def get_css_definition(self):
        """Gets a CSS fragment used for HTML representations of the shift."""
        return string.Template("""
            .$NAME {
                background-color: rgb($RED, $GREEN, $BLUE);
            }
        """).substitute(
            NAME=self.name,
            RED=self.rgb[0],
            GREEN=self.rgb[1],
            BLUE=self.rgb[2],
        )

    def get_day_td(self, date=None, daynum=None):
        """Gets a table cell (HTML td element) for this shift on a given date,
        used for the main rota depiction and the prototype."""
        if date is None and daynum is None:
            raise AssertionError("Specify either date or daynum")
        if date:
            bank_holiday = is_bank_holiday(date)
            weekend = is_weekend(date)
        else:
            bank_holiday = False
            weekend = daynum in [5, 6]
        if (bank_holiday or weekend) and self.nwd_only:
            classdetail = ' class="weekend"'
            content = "(off)"
        else:
            classdetail = ' class="{}"'.format(self.name)
            content = self.abbreviation
        return "<td{classdetail}>{content}</td>".format(
            classdetail=classdetail,
            content=content,
        )

    @staticmethod
    def get_shift_header_tr():
        """Returns a table header used for shift descriptions."""
        return """
            <tr>
                <th>Shift name</th>
                <th>Abbreviation</th>
                <th>Start time</th>
                <th>End time</th>
                <th>Duration (h)</th>
                <th>Shift type</th>
                <th>Work?</th>
                <th>Roles</th>
                <th>On normal working days only?</th>
                <th>Resident?</th>
            </tr>
        """

    def get_shift_tr(self):
        """Gets an HTML table row describing the shift."""
        return """
            <tr class="{name}">
                <td>{name}</td>
                <td>{abbreviation}</td>
                <td>{start_time}</td>
                <td>{end_time}</td>
                <td>{duration_h}</td>
                <td>{shift_type}</td>
                <td>{work}</td>
                <td>{roles}</td>
                <td>{nwd}</td>
                <td>{resident}</td>
            </tr>
        """.format(
            name=self.name,
            abbreviation=self.abbreviation,
            start_time=self.start_time.strftime("%H:%M"),
            end_time=(self.get_end_time().strftime("%H:%M")
                      if not self.invalid() else ""),
            duration_h=self.duration_h,
            shift_type=self.shift_type,
            work="Work" if self.work else "Off",
            roles=", ".join(self.roles),
            nwd="NWD only" if self.nwd_only else "",
            resident=yesno(self.resident),
        )

    def invalid(self):
        """Is the shift invalid, by virtue of not having a time/duration?"""
        return self.start_time is None or self.duration_h is None

    def get_end_time(self):
        """Get the shift's end time, as a datetime.time object, or None if it's
        invalid."""
        if self.invalid():
            return None
        interval = self.get_interval(ARBITRARY_DATE)
        return interval.end.time()

    def includes_7pm_to_7am(self):
        """Does the shift include any part between 7pm and 7am?"""
        if self.invalid():
            return False
        interval = self.get_interval(ARBITRARY_DATE)
        (day, night) = interval.day_night_duration(daybreak=datetime.time(7),
                                                   nightfall=datetime.time(19))
        return night.total_seconds() > 0

    def get_interval(self, startdate):
        """Returns an Interval representing the shift pattern, when the shift
        starts on the specified date. Returns None if it's invalid."""
        if self.invalid():
            return None
        start = datetime.datetime.combine(startdate, self.start_time)
        end = start + datetime.timedelta(hours=self.duration_h)
        return Interval(start, end)

    def get_work_interval(self, date, ignore_bank_holidays=True):
        """Returns an Interval representing the shift pattern, when the shift
        starts on the specified date, BUT ONLY if the shift involves actual
        work (i.e. not if it's a non-work shift, and not if it's a NWD shift
        on a weekend/Bank Holiday). Returns None if it's invalid."""
        if not self.work:
            return None
        if self.nwd_only and is_weekend(date):
            return None
        if (self.nwd_only and is_bank_holiday(date)
                and not ignore_bank_holidays):
            return None
        return self.get_interval(date)

    def roles_include(self, role):
        """Does the shift include the specified work role?"""
        return role in self.roles


# =============================================================================
# Doctor definition
# =============================================================================

class Doctor(object):
    """Object representing a doctor."""

    def __init__(self, name, daypattern, leave_weeks_per_year, ooh_roles=None):
        """
        Initializes the Doctor object.

        name: generic name for this doctor, e.g. "SHO 1"
        daypattern: list of Shift objects
        """
        self.name = name
        self.daypattern = list(daypattern)
        self.leave_weeks_per_year = leave_weeks_per_year
        self.ooh_roles = [] if ooh_roles is None else list(ooh_roles)
        # Check no shifts overlap
        intervals = []
        for i in range(len(daypattern)):
            startdate = ARBITRARY_DATE + datetime.timedelta(days=i)
            interval = daypattern[i].get_interval(startdate)
            if interval is not None:
                intervals.append(interval)
        il = IntervalList(intervals)
        if il.any_overlap():
            raise ValueError(
                "Doctor creation failed: overlapping shifts: {}".format(self))

    def __repr__(self):
        return (
            "Doctor(name={}, daypattern={}, leave_weeks_per_year={}, "
            "ooh_roles={})".format(
                repr(self.name),
                repr(self.daypattern),
                self.leave_weeks_per_year,
                self.ooh_roles if self.ooh_roles else None,
            )
        )

    def pad(self, target_length, shifts):
        """Pads the pattern of shifts to the target length (in days), using
        the sequence of shifts specified."""
        index = 0
        while len(self.daypattern) < target_length:
            self.daypattern.append(shifts[index])
            index += 1
            index %= len(shifts)

    def truncate(self, target_length):
        """Truncates the pattern of shift to the target length (in days)."""
        self.daypattern = self.daypattern[:target_length]

    def set_at(self, at, shifts, zero_based=False):
        n = len(shifts)
        zidx = at if zero_based else at - 1
        if zidx < 0 or zidx >= len(self.daypattern):
            raise ValueError("Bad 'at' parameter")
        for i in range(n):
            if zidx + i < len(self.daypattern):
                self.daypattern[zidx + i] = shifts[i]
            else:
                self.daypattern.append(shifts[i])

    def get_pattern_length(self):
        """Returns the length of the pattern, in days."""
        return len(self.daypattern)

    def ooh_roles_include(self, role):
        """Do the doctor's person-defined out-of-hours roles include the
        specified role?"""
        return role in self.ooh_roles

    def rotate(self, n=0):
        """Rotates the pattern forwards n days."""
        if n >= 0:
            self.daypattern = self.daypattern[-n:] + self.daypattern[:-n]
        else:
            # backwards
            n = -n
            self.daypattern = self.daypattern[n:] + self.daypattern[:n]

    def get_shift_for_day(self, dayindex):
        """Returns a shift for a given day (as a zero-based offset from the
        start of the rota)."""
        return self.daypattern[dayindex % self.get_pattern_length()]

    def get_shift_for_interval(self, interval, rota_start_date):
        """Returns the shift that the doctor will be working during the
        specified interval, or None if none can be found."""
        rota_start_dt = datetime.datetime.combine(rota_start_date,
                                                  datetime.time())
        if interval.start < rota_start_dt:
            raise AssertionError("Can't calculate shift for before rota start")
        # Start with the shift beginning on the day of the interval, and work
        # backwards
        shift_date = interval.start.date()
        days_into_rota = (interval.start.date() - rota_start_date).days
        l = self.get_pattern_length()
        pattern_idx = days_into_rota % l
        count = 0
        while count < l:
            s = self.daypattern[pattern_idx]
            s_int = s.get_interval(shift_date)
            if s_int is not None and s_int.overlaps(interval):
                return s
            # Otherwise, move on
            count += 1
            shift_date -= datetime.timedelta(days=1)
            pattern_idx -= 1
            pattern_idx %= l
        return None

    def copy(self, name=None, rotate=0):
        """Makes a copy of the Doctor object, optionally renaming it and
        rotating the shift pattern. Used to make lots of doctors with the same
        pattern as the starting point."""
        myclone = Doctor(name or self.name, self.daypattern,
                         self.leave_weeks_per_year, self.ooh_roles)
        l = self.get_pattern_length()
        rotate = rotate % l
        myclone.rotate(rotate)
        return myclone

    def get_all_times(self, rota_start_time):
        """Returns all shift start/end times, as a list of datetime.datetime
        objects."""
        times = set()
        for d, p in enumerate(self.daypattern):
            # Summer time ignored; see test()
            if p.invalid():
                continue
            interval = p.get_interval(rota_start_time
                                      + datetime.timedelta(days=d))
            times.add(interval.start)
            times.add(interval.end)
        return times

    def gen_daynum_idx_date(self, rota_start_date, rota_end_date):
        ndays = (rota_end_date - rota_start_date).days + 1
        l = len(self.daypattern)
        for daynum in range(ndays):
            idx = daynum % l
            date = rota_start_date + datetime.timedelta(days=daynum)
            yield (daynum, idx, date)

    def get_nwd_coverage(self, nwd_shifts, rota_start_date, rota_end_date):
        """On what proportion of normal working days (not weekends, not Bank
        Holidays) will this doctor be working a NWD pattern (one in the
        nwd_shifts list given)?"""
        yes = 0
        total = 0
        for (daynum, idx, date) in self.gen_daynum_idx_date(rota_start_date,
                                                            rota_end_date):
            shift = self.get_shift_for_day(daynum)
            if is_normal_working_day(date):
                total += 1
                if shift in nwd_shifts:
                    yes += 1
        return yes / total

    def get_work_pattern(self, start_date, end_date):
        """Returns an IntervalList with all working hours.
        Amalgamates adjacent periods of work.
        Treats Bank Holidays as work, because they are a statutory leave
        requirement in addition to annual leave, and leave isn't counted for
        the purposes of calculating working hours."""
        il = IntervalList(None)
        for (daynum, idx, date) in self.gen_daynum_idx_date(start_date,
                                                            end_date):
            shift = self.get_shift_for_day(daynum)
            il.add(shift.get_work_interval(date))
        return il

    def banding_info(self, rota):
        """Useful information for banding calculations. Returns the tuple:
        (
            hours_per_week,
            hours_per_week_inc_pc,
            fraction_weekends_worked,
            fraction_hours_ooh,
            max_continuous_hours,
            max_continuous_hours_interval,
            min_rest_between_duties_h,
            min_rest_between_duties_interval,
            max_continuous_duty_days,
            max_continuous_duty_interval,
            off_24h_every_7d,
            off_24h_every_7d_ffi,
            off_48h_every_14d,
            off_48h_every_14d_ffi,
            one_48h_and_one_62h_off_every_21d,
            one_48h_and_one_62h_off_every_21d_ffi,
            one_48h_and_one_62h_off_every_28d,
            one_48h_and_one_62h_off_every_28d_ffi,
        )
        This is the SLOWEST part.
        """
        logger.debug("  ... banding calcs part 1 for {}".format(self.name))
        rota_interval = Interval.dayspan(rota.start_date, rota.end_date,
                                         include_end=True)
        rota_weeks = convert_duration(rota_interval.duration(), 'w')
        workpattern = self.get_work_pattern(rota.start_date, rota.end_date)
        total_work_hours = convert_duration(workpattern.total_duration(), 'h')
        rota_weekends = rota_interval.n_weekends()
        weekends_worked = workpattern.n_weekends()
        ooh = convert_duration(
            workpattern.duration_outside_nwh(
                starttime=datetime.time(7),
                endtime=datetime.time(19)
            ), 'h')

        hours_per_week = total_work_hours / rota_weeks
        fraction_weekends_worked = weekends_worked / rota_weekends
        fraction_hours_ooh = ooh / total_work_hours

        if rota.prospective_cover:
            # the alternative Riddell formula uses the rota cycle rather than
            # the whole rota, but is equivalent; equivalents given
            leave_weeks = (self.leave_weeks_per_year / 52) * rota_weeks  # E
            non_leave_weeks = rota_weeks - leave_weeks  # B - E
            hours_if_no_leave_taken = total_work_hours  # D
            hours_per_leave_week_not_covered_by_others = (
                rota.hours_per_leave_week
            )  # C
            hours_removed_by_leave = (
                leave_weeks
                * hours_per_leave_week_not_covered_by_others
            )  # E * C
            hours_worked = (
                hours_if_no_leave_taken
                - hours_removed_by_leave
            )  # D - (E * C)
            hours_per_week_inc_pc = (
                hours_worked
                / non_leave_weeks
            )  # (D - (E * C)) / (B - E)
        else:
            hours_per_week_inc_pc = hours_per_week

        max_continuous_hours_interval = workpattern.longest_interval()
        max_continuous_hours = \
            max_continuous_hours_interval.duration_in('h')
        min_rest_between_duties_interval = workpattern.shortest_gap()
        min_rest_between_duties_h = \
            min_rest_between_duties_interval.duration_in('h')
        (max_continuous_duty_days,
         max_continuous_duty_interval) = workpattern.max_continuous_days()

        # _ffi: first failure interval
        (off_24h_every_7d, off_24h_every_7d_ffi) = workpattern.sufficient_gaps(
            every_n_days=7,
            requiredgaps=[
                datetime.timedelta(hours=24),
            ]
        )
        (off_48h_every_14d,
         off_48h_every_14d_ffi) = workpattern.sufficient_gaps(
            every_n_days=14,
            requiredgaps=[
                datetime.timedelta(hours=48),
            ]
        )
        (one_48h_and_one_62h_off_every_21d,
         one_48h_and_one_62h_off_every_21d_ffi) = workpattern.sufficient_gaps(
            every_n_days=21,
            requiredgaps=[
                datetime.timedelta(hours=62),
                datetime.timedelta(hours=48),
            ]
        )
        if one_48h_and_one_62h_off_every_21d:
            one_48h_and_one_62h_off_every_28d = True  # by definition
            one_48h_and_one_62h_off_every_28d_ffi = None
        else:
            (
                one_48h_and_one_62h_off_every_28d,
                one_48h_and_one_62h_off_every_28d_ffi
            ) = workpattern.sufficient_gaps(
                every_n_days=28,
                requiredgaps=[
                    datetime.timedelta(hours=62),
                    datetime.timedelta(hours=48),
                ]
            )

        return (
            hours_per_week,
            hours_per_week_inc_pc,
            fraction_weekends_worked,
            fraction_hours_ooh,
            max_continuous_hours,
            max_continuous_hours_interval,
            min_rest_between_duties_h,
            min_rest_between_duties_interval,
            max_continuous_duty_days,
            max_continuous_duty_interval,
            off_24h_every_7d,
            off_24h_every_7d_ffi,
            off_48h_every_14d,
            off_48h_every_14d_ffi,
            one_48h_and_one_62h_off_every_21d,
            one_48h_and_one_62h_off_every_21d_ffi,
            one_48h_and_one_62h_off_every_28d,
            one_48h_and_one_62h_off_every_28d_ffi,
        )

    def get_shift_types(self):
        shift_types = set()
        for shift in self.daypattern:
            shift_types.add(shift.shift_type)
        shift_types = [x for x in shift_types if x is not None]
        return shift_types

    def on_call_frequency(self):
        n = 0
        for shift in self.daypattern:
            if shift.shift_type == SHIFT_TYPES.ONCALL:
                n += 1
        return n / len(self.daypattern)

    def get_banding(self, rota, banding_info):
        """What banding should this doctor be paid?
        Takes as parameters:
            rota -- the rota
            banding_info -- the results of self.banding_info(rota)
                ... passed in manually for speed
        Returns a tuple: (banding description, decision list)"""
        decisions = []

        # Calculations used for the flowchart: preamble
        resident_after_7pm = any(
            s.resident and s.includes_7pm_to_7am()
            for s in self.daypattern
        )
        (
            hours_per_week,
            hours_per_week_inc_pc,
            fraction_weekends_worked,
            fraction_hours_ooh,
            max_continuous_hours,
            max_continuous_hours_interval,
            min_rest_between_duties_h,
            min_rest_between_duties_interval,
            max_continuous_duty_days,
            max_continuous_duty_interval,
            off_24h_every_7d,
            off_24h_every_7d_ffi,
            off_48h_every_14d,
            off_48h_every_14d_ffi,
            one_48h_and_one_62h_off_every_21d,
            one_48h_and_one_62h_off_every_21d_ffi,
            one_48h_and_one_62h_off_every_28d,
            one_48h_and_one_62h_off_every_28d_ffi,
        ) = banding_info
        on_call_freq = self.on_call_frequency()
        shift_types = self.get_shift_types()
        if shift_types:
            main_shift_type = shift_types[0]
            decisions.append("Rota type: {}".format(main_shift_type))
            if len(shift_types) > 1:
                decisions.append(
                    "WARNING: HYBRID ROTA; CONFUSED; ASSUMING {}".format(
                        main_shift_type))
        else:
            decisions.append("WARNING: NO SHIFT TYPES")
            main_shift_type = None

        # Rest compliance (New Deal/Working Time Regulations)
        rest_compliant = True
        # (a) Working hours
        if hours_per_week_inc_pc > 56:
            rest_compliant = False
            decisions.append("ND fail: >56 h/week on average (inc. PC)")
        if hours_per_week_inc_pc > 48:
            # rest_compliant = False
            decisions.append("WARNING: WTR fail: >48 h/week on average "
                             "(inc. PC); PROCEEDING ANYWAY.")
        if max_continuous_hours > 13:
            rest_compliant = False
            decisions.append(
                "WTR fail: >13 continuous hours (would need compensatory rest)"
                " (first failure: {})".format(max_continuous_hours_interval))
        if main_shift_type in [SHIFT_TYPES.FULL, SHIFT_TYPES.PARTIAL,
                               SHIFT_TYPES.PARTIAL24]:
            if min_rest_between_duties_h < 11:
                rest_compliant = False
                decisions.append(
                    "WTR fail: <11 h rest between duties (would need "
                    "compensatory rest) (first failure: {}".format(
                        min_rest_between_duties_interval))
        elif main_shift_type == SHIFT_TYPES.ONCALL:
            if min_rest_between_duties_h < 12:
                rest_compliant = False
                decisions.append(
                    "ND/WTR fail: <12 h rest between duties "
                    "(first failure: {})".format(
                        min_rest_between_duties_interval))
        # (b) Time off
        if main_shift_type in [SHIFT_TYPES.FULL, SHIFT_TYPES.PARTIAL,
                               SHIFT_TYPES.PARTIAL24]:
            if not one_48h_and_one_62h_off_every_28d:
                rest_compliant = False
                decisions.append(
                    "ND fail: does not provide one period of 48h off and one "
                    "period of 62h off every 28 days "
                    "(first failure: {})".format(
                        one_48h_and_one_62h_off_every_28d_ffi))
        elif main_shift_type == SHIFT_TYPES.ONCALL:
            if not one_48h_and_one_62h_off_every_21d:
                rest_compliant = False
                decisions.append(
                    "ND fail: does not provide one period of 48h off and one "
                    "period of 62h off every 21 days "
                    "(first failure: {})".format(
                        one_48h_and_one_62h_off_every_21d_ffi))
        if not(off_24h_every_7d or off_48h_every_14d):
            rest_compliant = False
            msg = "WTR fail: not 24h off every 7d, or 48h off every 14d"
            if not off_24h_every_7d:
                msg += " (24h q7d first failure: {})".format(
                    off_24h_every_7d_ffi)
            if not off_48h_every_14d:
                msg += " (48h q14d first failure: {})".format(
                    off_48h_every_14d_ffi)
            decisions.append(msg)
        if max_continuous_duty_days > 13:
            rest_compliant = False
            decisions.append(
                "WTR fail: >13 continuous duty days "
                "(first failure: {})".format(max_continuous_duty_interval))
        # (c) Rest at work
        decisions.append("Natural breaks assumed OK")
        decisions.append("Minimum rest assumed OK")
        decisions.append("Timing of continuous rest assumed OK")

        # Other calculations used for the flowchart
        more_than_48h = hours_per_week_inc_pc > 48
        on_call_rota = "On-call" in shift_types
        if on_call_rota:
            decisions.append("WARNING: DON’T UNDERSTAND PROSPECTIVE COVER")
        one_in_six_or_more_inc_pc = on_call_freq >= 1/6
        one_in_eight_or_more_inc_pc = on_call_freq >= 1/8
        one_in_eight_without_pc_or_less = on_call_freq < 1/8
        # ... these three depend on the prospective cover calculation somehow
        one_weekend_in_three_or_more = fraction_weekends_worked >= 1/3
        one_weekend_in_four_or_more = fraction_weekends_worked >= 1/4
        more_than_third_hours_outside_normal = fraction_hours_ooh > 1/3
        fulfil_criterion_r = (
            resident_after_7pm
            or rota.work4h_after7pm_halformore
        )
        criterion_r_reasoning = (
            " [CRITERION R WORKING: Resident after 7pm? "
            + yesno(resident_after_7pm)
            + " / 4h work after 7pm on >=50% occasions, from monitoring? "
            + yesno(rota.work4h_after7pm_halformore)
            + "]"
        )
        resident = any(s.resident for s in self.daypattern)

        # The flowchart from the BMA's guide [BMA_1]:
        band = None
        decisions.append("Complies with ND/WTR rest requirements? "
                         + yesno(rest_compliant))
        if not rest_compliant:
            return ("3", decisions)  # BAND 3
        decisions.append(">48h of actual work (inc. PC)? "
                         + yesno(more_than_48h))
        if more_than_48h:
            decisions.append("On call rota? " + yesno(on_call_rota))
            if on_call_rota:
                decisions.append(
                    "1:6 or more inc. PC? "
                    + yesno(one_in_six_or_more_inc_pc)
                    + " / 1:3 weekends or more? "
                    + yesno(one_weekend_in_three_or_more)
                )
                if one_in_six_or_more_inc_pc or one_weekend_in_three_or_more:
                    decisions.append("Fulfil criterion R? "
                                     + yesno(fulfil_criterion_r)
                                     + criterion_r_reasoning)
                    if fulfil_criterion_r:
                        return ("2A", decisions)  # BAND 2A
                    else:
                        return ("2B", decisions)  # BAND 2B
                else:
                    return ("2B", decisions)  # BAND 2B
            else:
                decisions.append(
                    ">1/3 hours outside 7am-7pm Mon-Fri? "
                    + yesno(more_than_third_hours_outside_normal)
                    + " / 1:3 weekends or more? "
                    + yesno(one_weekend_in_three_or_more)
                )
                if (more_than_third_hours_outside_normal
                        or one_weekend_in_three_or_more):
                    return ("2A", decisions)  # BAND 2A
                else:
                    return ("2B", decisions)  # BAND 2B
        else:
            decisions.append("On call rota? " + yesno(on_call_rota))
            if not on_call_rota:
                decisions.append(
                    ">1/3 hours outside 7am-7pm Mon-Fri? "
                    + yesno(more_than_third_hours_outside_normal)
                    + " / 1:4 weekends or more? "
                    + yesno(one_weekend_in_four_or_more)
                )
                if (more_than_third_hours_outside_normal
                        or one_weekend_in_four_or_more):
                    return ("1A", decisions)  # BAND 1A
                else:
                    return ("1B", decisions)  # BAND 1B
            else:
                decisions.append(
                    "1:6 or more inc. PC? "
                    + yesno(one_in_six_or_more_inc_pc)
                )
                if one_in_six_or_more_inc_pc:
                    return ("1A", decisions)  # BAND 1A
                else:
                    decisions.append(
                        "1:8 or more inc. PC? "
                        + yesno(one_in_eight_or_more_inc_pc)
                        + " / 1:4 weekends or more? "
                        + yesno(one_weekend_in_four_or_more)
                    )
                    if (one_in_eight_or_more_inc_pc
                            or one_weekend_in_four_or_more):
                        decisions.append("Fulfil criterion R? "
                                         + yesno(fulfil_criterion_r)
                                         + criterion_r_reasoning)
                        if fulfil_criterion_r:
                            return ("1A", decisions)  # BAND 1A
                        else:
                            return ("1B", decisions)  # BAND 1B
                    else:
                        decisions.append(
                            "1:8 without PC or less? "
                            + yesno(one_in_eight_without_pc_or_less)
                        )
                        if not one_in_eight_without_pc_or_less:
                            return ("1B", decisions)  # BAND 1B
                        else:
                            decisions.append(
                                "Resident (clinically/contractually)? "
                                + yesno(resident)
                            )
                            if resident:
                                return ("1B", decisions)  # BAND 1B
                            else:
                                return ("1C", decisions)  # BAND 1C
        return (band, decisions)


# =============================================================================
# Standalone functions for parallel processing. NOT YET WORKING.
# =============================================================================

def pfunc_doctor_banding(args):
    # logger.debug("args={}".format(args))
    rota = args[0]
    doc = args[1]
    # logger.debug("rota={}, doc={}".format(rota, doc))
    return doc.banding_info(rota)


def pfunc_role_coverage(args):
    rota = args[0]
    role = args[1]
    return rota.get_role_coverage(role)


# =============================================================================
# Rota definition
# =============================================================================

class Rota(object):
    """Object representing a rota."""

    def __init__(self, name, shifts, doctors, nwd_shifts,
                 prototypes=None,
                 start_date=datetime.date(2015, 8, 5),
                 end_date=datetime.date(2016, 2, 2),
                 prospective_cover=True,
                 prototype_rotation_to_monday_start=2,
                 work4h_after7pm_halformore=True,
                 hours_per_leave_week=40,
                 comments=None):
        """
        Initializes the rota.

        name: rota name
        shifts: list of Shift objects; all shifts used by the rota
        doctors: list of Doctor objects
        nwd_shifts: list of Shift objects; all normal-working-day shifts
        prototypes: list of Doctor objects with prototype day patterns
        start_date: date the rota starts
        end_date: last date of the rota
        prospective_cover: is prospective cover in use?
        prototype_rotation_to_monday_start: rotate() parameter required to get
            the prototype daypattern to start on a Monday (default 2, for
            patterns that start on a Wednesday).
        work4h_after7pm_halformore: does monitoring show that the doctors work
            4 hours after 7pm on half or more occasions?
        hours_per_leave_week: usually 40; the number of hours in a 'base' week,
            used for leave calculations
        comments: a list of HTML objects to be inserted in a <ul> of comments,
            or None.
        """
        self.name = name
        self.shifts = shifts
        self.doctors = doctors
        self.prototypes = [] if prototypes is None else prototypes
        self.nwd_shifts = nwd_shifts
        self.start_date = start_date
        self.end_date = end_date
        self.prospective_cover = prospective_cover
        self.prototype_rotation_to_monday_start = \
            prototype_rotation_to_monday_start
        self.work4h_after7pm_halformore = work4h_after7pm_halformore
        self.hours_per_leave_week = hours_per_leave_week
        self.comments = [] if comments is None else comments
        # Derived:
        self.start_time = datetime.datetime.combine(
            start_date, datetime.time())  # midnight on first day
        self.n_days = (end_date - start_date).days + 1

    def __repr__(self):
        return (
            "Rota(name={}, shifts={}, doctors={}, nwd_shifts={}, "
            "start_date={}, end_date={}, prospective_cover={}, "
            "work4h_after7pm_halformore={}, hours_per_leave_week={})".format(
                repr(self.name), self.shifts, self.doctors, self.nwd_shifts,
                self.start_date, self.end_date, self.prospective_cover,
                self.work4h_after7pm_halformore, self.hours_per_leave_week,
            )
        )

    def get_all_roles(self):
        """What on-call roles are defined in the rota's shifts, or on a person-
        specific basis for the doctors?"""
        roles = set()
        for s in self.shifts:
            roles.update(s.roles)
        for d in self.doctors:
            roles.update(d.ooh_roles)
        roles = list(roles)
        roles.sort()
        return roles

    def get_role_coverage(self, role):
        """For a given role, returns the tuple containing:
        (a) what proportion of time is the role fulfilled by at least one
            doctor? (Target is often 1.)
        (b) what proportion of time is the role fulfilled by more than one
            doctor? (Target is often close to 0, but note that some handover
            time is sensible, and will make this non-zero unless you build
            a separate handover shift, which is visually complex).
        """
        logger.debug("... role coverage for: {}".format(role))
        # 1. Build a list of intervals covering that role.
        coverage = IntervalList(no_overlap=False, no_contiguous=False)
        # ... an IntervalList that's happy to maintain overlapping intervals
        for d in self.doctors:
            for (daynum, idx, date) in d.gen_daynum_idx_date(self.start_date,
                                                             self.end_date):
                s = d.get_shift_for_day(daynum)
                if s is None:
                    continue
                if s.roles_include(role) or (s.work and not s.nwd_only and
                                             d.ooh_roles_include(role)):
                    # Shift or doctor provide the role
                    coverage.add(s.get_interval(date))
        # 2. Analyse.
        nonoverlapping_coverage = coverage.copy(no_overlap=True,
                                                no_contiguous=True)
        whole_rota_interval = Interval.dayspan(self.start_date, self.end_date)
        total_duration_s = whole_rota_interval.duration_in('s')
        coverage_duration = nonoverlapping_coverage.total_duration()
        duration_at_least_one_s = convert_duration(coverage_duration, 's')
        overlaps = coverage.get_overlaps()
        overlap_duration = overlaps.total_duration()
        duration_more_than_one_s = convert_duration(overlap_duration, 's')
        # 3. Done.
        prop_at_least_one = duration_at_least_one_s / total_duration_s
        prop_more_than_one = duration_more_than_one_s / total_duration_s
        return (prop_at_least_one, prop_more_than_one)

    def get_all_times(self):
        """Returns a list of datetime.datetime objects representing shift
        start/end times for all doctors in the rota."""
        times = set()
        for d in self.doctors:
            times.update(d.get_all_times(self.start_time))
        times = list(times)
        times.sort()
        return times

    def print_html(self, filename, daynums=False, skipanalytics=False,
                   noprototypes=False):
        """Write the rota summary/analysis to an HTML file."""
        logger.info("Writing HTML for rota {} to {}...".format(self.name,
                                                               filename))
        with open(filename, 'w') as f:
            print(
                """
                    <!DOCTYPE html>
                    <html>
                        <head>
                            <title>{title}</title>
                            <meta charset="utf-8">
                            <style type="text/css">
                                {css}
                            </style>
                        </head>
                        <body>
                            {body}
                        </body>
                    </html>
                """.format(title=webify(self.name),
                           css=self.get_css(),
                           body=self.get_html_body(daynums, skipanalytics,
                                                   noprototypes)),
                file=f
            )
        logger.info("HTML written.")

    def get_css(self):
        """Returns the CSS for the rota display."""
        css = """
            body {
                font-family: "Times New Roman", Georgia, Serif;
            }
            h1, h2, h3 {
                font-family: Arial, Helvetica, sans-serif;
            }
            h1 {
                font-size: 1.4em;
            }
            h2 {
                font-size: 1.2em;
            }
            h3 {
                font-size: 1.0em;
            }
            table {
                border-collapse: collapse;
                border: 1px solid black;
                padding: 0px;
            }
            tr, th, td {
                vertical-align: top;
                text-align: left;
                margin: 0px;
                padding: 1px;
                border: 1px solid black;
            }
            .endmatter {
                background-color: rgb(225, 225, 225);
            }
            .weekend {
                background-color: rgb(200, 200, 200);
            }
        """
        for s in self.shifts:
            css += s.get_css_definition()
        return css

    def get_html_body(self, daynums=False, skipanalytics=False,
                      noprototypes=False):
        """Returns the main part of the HTML display."""
        html = (
            "<h1>{}</h1>".format(self.name)
            + self.get_html_comments()
            + self.get_html_rota_settings()
            + self.get_html_shifts()
            + self.get_html_rota_pattern(daynums)
        )
        if not noprototypes:
            html += self.get_html_prototypes()
        if not skipanalytics:
            html += (
                self.get_html_role_coverage()
                + self.get_html_nwd_coverage()
                + self.get_html_bandings()
            )
        html += self.get_html_footnotes()
        return html

    def get_html_comments(self):
        """Adds in the comments, which should already be in valid HTML
        format."""
        if not self.comments:
            return ""
        html = """
            <h2>Comments</h2>
            <ul>
        """
        for c in self.comments:
            html += "<li>{}</li>".format(c)
        html += "</ul>"
        return html

    def get_html_footnotes(self):
        """HTML for footnotes."""
        return """
            <div class="endmatter">
            <h2>End matter</h2>
            <h3>Footnotes</h3>
            <ul>
                <li>Normal working day coverage: the closer this is to 100%,
                (a) the better wards and clinics run, and (b) the less
                disruption there is to anyone who isn’t employed full-time
                by the NHS, such as LTFT and academic trainees.</li>

                <li>Prospective cover: if you want to take leave when you’re
                on call and your employer will always book a locum for this,
                then you’re not working prospective cover. Normally, you’d
                have to swap shifts, meaning that you and your colleagues
                cover for each other’s absences due to annual leave; this is
                prospective cover, and is the norm.</li>

                <li>Time off rules: the interpretation used here of “24h off
                every 7 days” is that one can work 7 days in a row and the 24h
                off must occur immediately afterwards (not that there must be
                a 24h period off <i>within</i> every 7-day period). In this
                software, this is implemented by the “flexibility=2” option
                (see source). Support for this interpretation is given by the
                RCP London’s (2006, <i>Designing safer rotas for junior
                doctors</i>) explanation that a rota involving 7 night shifts
                (each 21:00–09:30) is legitimate (meaning New Deal and EWTD
                compliant) — albeit heavily discouraged on safety grounds.</li>

                <li>Banding supplements are as follows.
                Band 3: 100%.
                Band 2A: 80%.
                Band 2B: 50%.
                Band 1A: 50%.
                Band 1B: 40%.
                Band 1C: 20%.</li>

                <li>In Aug 2015, an “average” SHO salary is about £31,838
                (the midpoint of the first three scale points of the specialty
                registrar (CT/ST) grade, since core training is CT1–3), and an
                “average” SpR salary is about £37,822
                (the midpoint of the next three scale points, since higher
                specialist training is usually ST4–6) [BMA_6].
                On-call supplements don’t cost the Trust employer’s pension
                contributions, so the real cost is relatively close to the
                headline cost (I’m not sure if there are other employer
                contributions of relevance).</li>

                <li>Band 3 includes all posts that are non-compliant with the
                New Deal, and should not be in use.</li>

                <li>Band 2 includes all posts that are compliant with the New
                Deal, but where the hours of work per week are >48 and <56.
                They are not compliant with the EWTD (which stipulates a
                48h/week limit). They should no longer be in use [OHFP].</li>
            </ul>

            <h3>References</h3>
            <ul>
                <li>[BMA_1] BMA: <a href="{bma_1}">Are you being paid
                correctly? BMA pay guidance for junior doctors</a></li>

                <li>[BMA_2] BMA: <a href="{bma_2}">Junior doctors rest
                requirements under New Deal and Working Time
                Regulations</a></li>

                <li>[BMA_3] <a href="{bma}">BMA</a> → Home → Practical support
                at work → Contracts → Junior contracts → Rotas and working
                patterns → <a href="{bma_3}">Rota glossary</a></li>

                <li>[BMA_4] BMA: <a href="{bma_4}">EWTD opt-out</a></li>

                <li>[BMA_5] BMA (2008): <a href="{bma_5}">The final countdown:
                the rush to reband training posts explained</a></li>

                <li>[BMA_6] <a href="{bma}">BMA</a> → Home → Practical support
                at work → Pays, fees &amp; allowances → Pay scales
                → <a href="{bma_6}">Junior doctors England</a></li>

                <li>[BMA_7] <a href="{bma}">BMA</a> → Home → Practical support
                at work → Contracts → Junior contracts → Rotas and working
                patterns → <a href="{bma_7}">Riddell formula
                calculator</a></li>

                <li>[NHS_1] <a href="nhs_1">NHS Employers (2002–2013)
                Terms and Conditions of Service: NHS Medical and Dental Staff
                (England)</a></li>

                <li>[RCP_1] RCP London (2006) <a href="{rcp_1}">Designing safer
                rotas for junior doctors</a></li>

                <li>[OHFP] Oxford Handbook for the Foundation Programme, 2014,
                fourth edition. Oxford University Press.</li>
            </ul>

            <h3>Abbreviations</h3>
            <ul>
                <li>BMA = British Medical Association</li>
                <li>EWTD = European Working Time Directive</li>
                <li>LTFT = less than full time</li>
                <li>ND = New Deal</li>
                <li>WTR = Working Time Regulations</li>
            </ul>

            <p><i>Rota validation software by Rudolf Cardinal. Version
            {VERSION}. No guarantees as to accuracy; check with your Medical
            Staffing department as well. Software is open-source and available
            at <a href="{url}">{url}</a>.</p>

            </div>
        """.format(
            url="https://github.com/RudolfCardinal/rota",
            VERSION=VERSION,
            bma="http://bma.org.uk/",
            bma_1="http://bma.org.uk/-/media/files/pdfs/practical%20advice%20at%20work/your%20rights/pay%20fees%20allowances/are%20you%20being%20paid%20correctly.pdf",  # noqa
            bma_2="http://bma.org.uk/-/media/files/pdfs/practical%20advice%20at%20work/contracts/jdrestreqs.pdf",  # noqa
            bma_3="http://bma.org.uk/practical-support-at-work/contracts/juniors-contracts/rotas-and-working-patterns/rota-glossary",  # noqa
            bma_4="http://bma.org.uk/practical-support-at-work/ewtd/ewtd-juniors/ewtd-opt-out",  # noqa
            bma_5="https://bma.org.uk/-/media/files/pdfs/practical%20advice%20at%20work/your%20rights/pay%20fees%20allowances/finalcountdown_ewtd.pdf",  # noqa
            bma_6="http://bma.org.uk/practical-support-at-work/pay-fees-allowances/pay-scales/juniors-england",  # noqa
            bma_7="http://bma.org.uk/practical-support-at-work/contracts/juniors-contracts/rotas-and-working-patterns/rota-template-and-riddell-calculator",  # noqa
            rcp_1="https://www.rcplondon.ac.uk/sites/default/files/documents/designing_safer_rotasweb.pdf",  # noqa
            nhs_1="http://www.nhsemployers.org/~/media/Employers/Documents/Pay%20and%20reward/Terms_and_Conditions_of_Service_NHS_Medical_and_Dental_Staff_300813_bt.pdf",  # noqa
        )

    def get_html_rota_settings(self):
        """HTML for rota-wide settings."""
        return """
            <h2>Rota-wide settings</h2>
            <table>
                <tr>
                    <th>Setting</th>
                    <th>Value</th>
                </tr>
                <tr>
                    <td>Prospective cover</td>
                    <td>{prospective_cover}</td>
                </tr>
                <tr>
                    <td>Monitoring shows ≥4h work after 7pm on ≥50%
                        occasions</td>
                    <td>{work4h_after7pm_halformore}</td>
                </tr>
                <tr>
                    <td>Hours per leave week (... not covered by others)</td>
                    <td>{hours_per_leave_week}</td>
                </tr>
            </table>
        """.format(
            prospective_cover=yesno(self.prospective_cover),
            work4h_after7pm_halformore=yesno(self.work4h_after7pm_halformore),
            hours_per_leave_week=self.hours_per_leave_week,
        )

    def get_html_shifts(self):
        """HTML for shift report."""
        logger.info("- shifts")
        html = "<h2>Shifts</h2><table>"
        html += Shift.get_shift_header_tr()
        for s in self.shifts:
            html += s.get_shift_tr()
        html += "</table>"
        return html

    def get_html_rota_pattern(self, daynums=False):
        """HTML for main rota pattern."""
        logger.info("- rota")
        wholespan = Interval.dayspan(self.start_date, self.end_date)
        n_weeks = wholespan.duration_in('w')
        html = """
            <h2>Rota</h2>
            <p>Weeks in rota: {n_weeks}</p>
            <table>
        """.format(
            n_weeks=n_weeks,
        )
        # Header row
        dateheader = "Date"
        if daynums:
            dateheader += " <i>(Day number)</i>"
        html += "<tr><th>{}</th>".format(dateheader)
        for d in self.doctors:
            html += "<th>{}</th>".format(d.name)
        html += "</tr>\n"
        # Other rows
        for dayoffset in range(self.n_days):
            date = self.start_date + datetime.timedelta(days=dayoffset)
            datetext = date.strftime("%a %d %b %Y")
            if daynums:
                datetext += " <i>({})</i>".format(dayoffset + 1)
            wk = is_weekend(date)
            bh = is_bank_holiday(date)
            if bh:
                datetext += " <b>(Bank Hol)</b>"
            html += "<tr{}>".format(' class="weekend"' if wk or bh else "")
            html += "<td>{}</td>".format(datetext)
            for doctor in self.doctors:
                html += (
                    doctor.get_shift_for_day(dayoffset).get_day_td(date=date)
                )
            html += "</tr>\n"
        html += "</table>"
        return html

    def get_html_prototypes(self):
        """HTML for doctor prototypes."""
        if not self.prototypes:
            return ""
        html = "<h2>Doctor prototypes</h2>"
        for p in self.prototypes:
            d = p.copy(rotate=self.prototype_rotation_to_monday_start)
            html += """
                <h3>{name}</h3>
                <table>
                    <tr>
                        <th>Week</th>
                        <th>Mon</th>
                        <th>Tue</th>
                        <th>Wed</th>
                        <th>Thu</th>
                        <th>Fri</th>
                        <th class="weekend">Sat</th>
                        <th class="weekend">Sun</th>
                    </tr>
            """.format(
                name=d.name
            )
            for i in range(len(d.daypattern)):
                daynum = i % 7  # start on a Monday
                weeknum = i // 7 + 1  # // is integer division
                shift = d.daypattern[i]
                if daynum == 0:
                    html += "<tr><td>{}</td>".format(weeknum)
                html += shift.get_day_td(daynum=daynum)
                if daynum == 7:
                    html += "</tr>\n"
            html += "</table>\n"
        return html

    def get_html_role_coverage(self):
        """HTML for role coverage analysis."""
        logger.info("- role coverage")
        html = """
            <h2>Role coverage</h2>
            <table>
                <tr>
                    <th>Role</th>
                    <th>% time covered by at least 1 person</th>
                    <th>% time covered by &gt;1 person (but NB handover usually
                        ~2–4%)</th>
                </tr>
        """
        roles = self.get_all_roles()

        resultlist = []
        args = ((self, role) for role in roles)

        # Non-parallel version for debugging:
        # for result in map(pfunc_role_coverage, args):
        #     resultlist.append(result)

        # Parallel version (fewer error details if it fails):
        # http://stackoverflow.com/questions/6785226
        with concurrent.futures.ProcessPoolExecutor() as executor:
            for result in executor.map(pfunc_role_coverage, args):
                resultlist.append(result)

        for i in range(len(roles)):
            role = roles[i]
            result = resultlist[i]

            (at_least, more_than_one) = result
            html += """
                <tr>
                    <td>{role}</td>
                    <td>{at_least}</td>
                    <td>{more_than_one}</td>
                </tr>
            """.format(
                role=role,
                at_least=number_to_dp(100 * at_least, DP),
                more_than_one=number_to_dp(100 * more_than_one, DP),
            )
        html += "</table>"
        return html

    def get_html_nwd_coverage(self):
        """HTML for normal working day coverage."""
        logger.info("- NWD coverage")
        html = """
            <h2>Normal working day shift coverage</h2>
            <table>
                <tr>
                    <th>Doctor</th>
                    <th>% normal working days spent working a normal shift</th>
                </tr>
        """
        for d in self.doctors:
            nwd = 100 * d.get_nwd_coverage(self.nwd_shifts, self.start_date,
                                           self.end_date)
            html += """
                <tr>
                    <td>{name}</td>
                    <td>{nwd}</td>
                </tr>
            """.format(
                name=d.name,
                nwd=number_to_dp(nwd, DP),
            )
        html += "</table>"
        return html

    def get_html_bandings(self):
        """HTML for doctor/banding analysis."""
        logger.info("- bandings")
        html = """
            <h2>Bandings</h2>
            <table>
                <tr>
                    <th>Doctor</th>
                    <th>OOH person-defined roles</th>
                    <th>Hours/week (actual, ignoring leave/BHs)</th>
                    <th>Hours/week (including prospective cover)</th>
                    <th>Leave weeks/year</th>
                    <th>Fraction weekends worked</th>
                    <th>Fraction hours OOH</th>
                    <th>Max continuous hours</th>
                    <th>Min rest between duty</th>
                    <th>Max continuous duty days</th>
                    <th>24h off every 7d</th>
                    <th>48h off every 14d</th>
                    <th>48h+62h off every 21d</th>
                    <th>48h+62h off every 28d</th>
                    <th>Banding</th>
                    <th>Supplement</th>
                    <th>Working</th>
                </tr>
        """

        banding_info_list = []
        args = ((self, doc) for doc in self.doctors)

        # Non-parallel version for debugging:
        # for result in map(pfunc_doctor_banding, args):
        #     banding_info_list.append(result)

        # Parallel version (fewer error details if it fails):
        # http://stackoverflow.com/questions/6785226
        with concurrent.futures.ProcessPoolExecutor() as executor:
            for result in executor.map(pfunc_doctor_banding, args):
                banding_info_list.append(result)

        for i in range(len(self.doctors)):
            d = self.doctors[i]
            logger.debug("  ... banding calcs part 2 for {}".format(d.name))
            # banding_info = d.banding_info(self)
            banding_info = banding_info_list[i]

            (
                hours_per_week,
                hours_per_week_inc_pc,
                fraction_weekends_worked,
                fraction_hours_ooh,
                max_continuous_hours,
                max_continuous_hours_interval,
                min_rest_between_duties_h,
                min_rest_between_duties_interval,
                max_continuous_duty_days,
                max_continuous_duty_interval,
                off_24h_every_7d,
                off_24h_every_7d_ffi,
                off_48h_every_14d,
                off_48h_every_14d_ffi,
                one_48h_and_one_62h_off_every_21d,
                one_48h_and_one_62h_off_every_21d_ffi,
                one_48h_and_one_62h_off_every_28d,
                one_48h_and_one_62h_off_every_28d_ffi,
            ) = banding_info
            (banding, decisions) = d.get_banding(self, banding_info)
            html += """
                <tr>
                    <td>{name}</td>
                    <td>{ooh_roles}</td>
                    <td>{hours_per_week}</td>
                    <td>{hours_per_week_inc_pc}</td>
                    <td>{leave_weeks_per_year}</td>
                    <td>{fraction_weekends_worked}</td>
                    <td>{fraction_hours_ooh}</td>
                    <td>{max_continuous_hours}</td>
                    <td>{min_rest_between_duties_h}</td>
                    <td>{max_continuous_duty_days}</td>
                    <td>{off_24h_every_7d}</td>
                    <td>{off_48h_every_14d}</td>
                    <td>{one_48h_and_one_62h_off_every_21d}</td>
                    <td>{one_48h_and_one_62h_off_every_28d}</td>
                    <td>{banding}</td>
                    <td>{supplement}</td>
                    <td>{working}</td>
                </tr>
            """.format(
                name=d.name,
                ooh_roles=", ".join(d.ooh_roles),
                hours_per_week=number_to_dp(hours_per_week, DP),
                hours_per_week_inc_pc=number_to_dp(hours_per_week_inc_pc, DP),
                leave_weeks_per_year=d.leave_weeks_per_year,
                fraction_weekends_worked=number_to_dp(fraction_weekends_worked,
                                                      DP),
                fraction_hours_ooh=number_to_dp(fraction_hours_ooh, DP),
                max_continuous_hours=number_to_dp(max_continuous_hours, DP),
                min_rest_between_duties_h=number_to_dp(
                    min_rest_between_duties_h, DP),
                max_continuous_duty_days=max_continuous_duty_days,
                off_24h_every_7d=yesno(off_24h_every_7d),
                off_48h_every_14d=yesno(off_48h_every_14d),
                one_48h_and_one_62h_off_every_21d=yesno(
                    one_48h_and_one_62h_off_every_21d),
                one_48h_and_one_62h_off_every_28d=yesno(
                    one_48h_and_one_62h_off_every_28d),
                banding=banding,
                supplement=SUPPLEMENT[banding],
                working=" → ".join(decisions),
            )
        html += "</table>"
        return html


# =============================================================================
# Testing
# =============================================================================

def test():
    # In the UK in 2015, clocks:
    # - go forward one hour at 01:00 on 29 Mar 2015
    # - go backward one hour at 02:00 on 25 Oct 2015
    before_tz_change = datetime.datetime(2015, 3, 1)
    after_tz_change = before_tz_change + datetime.timedelta(days=120)
    print(before_tz_change)
    print(after_tz_change)
    # ... both are midnight, because these are not timezone-aware times
    # ... so that's OK. We'll ignore clock changes entirely.


# =============================================================================
# Specific rotas
# =============================================================================

# -----------------------------------------------------------------------------
# Examples/tests
# -----------------------------------------------------------------------------

def test_prospective_cover():
    # Shifts
    nwd = Shift(  # 4h am, 4h pm
        "Normal_working_day", "nwd", datetime.time(9), 8, nwd_only=True,
        resident=True, rgb=COLOURS.NWD)
    long = Shift(  # 6h am, 7h pm
        "Long", "L", datetime.time(6), 13,
        roles=["On call"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=COLOURS.LATE_F)
    night = Shift(  # 6h am, 6h pm
        "Pseudonight", "N", datetime.time(6), 12,
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=COLOURS.NIGHT_FA)
    off = Shift(
        "Off", "OFF", datetime.time(9, 15), 7.75, work=False,
        rgb=COLOURS.OFF)
    shifts = [nwd, long, night, off]

    # Doctors
    n_docs = 8
    base_doc = Doctor("Prototype doctor", [
        long, nwd, nwd, long, nwd, nwd, nwd,
        off, off, night, night, night, night, night,
        night, night, off, off, off, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, long, long, long,
        nwd, nwd, long, nwd, nwd, nwd, nwd,
        nwd, long, nwd, nwd, nwd, nwd, nwd,
    ], leave_weeks_per_year=6.5)
    doctors = []
    for i in range(n_docs):
        doctors.append(base_doc.copy("doc_" + str(i + 1), rotate=i*7))

    # Rota
    start_monday = ARBITRARY_MONDAY_NEAR_BH
    # start_monday = ARBITRARY_MONDAY_FAR_FROM_BH
    return Rota(
        "Test prospective cover calculations", shifts, doctors,
        start_date=start_monday,
        end_date=start_monday + datetime.timedelta(days=7 * n_docs - 1),
        hours_per_leave_week=40,
        nwd_shifts=[nwd],
        prototypes=[base_doc],
        comments=[
            "Tests prospective cover calculations.",
            "By Riddell formula, should be 47.86 h/week inc. PC.",
            "Rota hours mimic the example in [BMA_7] (which shows 47.9 h/week "
            "or 47.86 with rounding removed).",
            "<b>Should give identical hours/week for all doctors, even when "
            "some are working on bank holidays.</b> The test interval "
            "contains a bank holiday to check this.",
        ],
    )

# -----------------------------------------------------------------------------
# CPFT actual
# -----------------------------------------------------------------------------

def cpft_actual_aug2015_south():
    # Shifts
    # ... colours to match the original
    nwd = Shift(
        "Normal_working_day", "nwd", datetime.time(9), 8, nwd_only=True,
        resident=True, rgb=(217, 150, 148))
    late1 = Shift(
        "Late_1", "L1", datetime.time(9), 12.5,
        roles=["Late", "Late_1"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=(142, 180, 227))
    late2 = Shift(
        "Late_2", "L2", datetime.time(9), 12.5,
        roles=["Late", "Late_2"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=(0, 176, 240))
    late3 = Shift(
        "Late_3", "L3", datetime.time(9), 12.5,
        roles=["Late", "Late_3"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=(23, 55, 94))
    night1 = Shift(
        "Night_1", "N1", datetime.time(21, 15), 12,
        roles=["Night", "Night_1"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=(146, 208, 80))
    night2 = Shift(
        "Night_2", "N2", datetime.time(21, 15), 12,
        roles=["Night", "Night_2"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=(79, 98, 40))
    off = Shift(
        "Off", "OFF", datetime.time(9, 15), 7.75, work=False)
    shifts = [nwd, late1, late2, late3, night1, night2, off]

    # Doctors
    # ... so confusing to edit them manually! Just a list

    basedoc = Doctor("BASEDOC", [
        night1, night1, off, nwd, nwd,  # Mon/Tue missing
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, late1, late1, late1,
        off, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, off, night2, night2, night2,
        off, off, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        late2, late2, late2, late2, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, late3, late3, late3,
        off, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        late1, late1, late1, late1, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        late3, late3, late3, late3, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        night2, night2, night2, night2, off, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, late2, late2, late2,
        off, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        night1, night1
    ], leave_weeks_per_year=0)
    doctors = [basedoc]  # so we can operate with 1-indexing for a while
    for i in range(1, 26 + 1):
        if i in [1, 2, 8, 9, 10, 11, 12, 13, 14, 21, 25, 26]:
            name = str(i) + "_SpR"
            ooh_roles = ["Section 12 approved"]
            leave_weeks_per_year = 6
        else:
            name = str(i) + "_SHO"
            ooh_roles = []
            leave_weeks_per_year = 5
        name
        d = basedoc.copy(name=name, rotate=(i - 1) * 7)
        d.ooh_roles = ooh_roles
        d.leave_weeks_per_year = leave_weeks_per_year
        doctors.append(d)
    # Now the deviations from a simple pattern
    NWD = [nwd]
    LATE1_MTWT = [late1, late1, late1, late1]
    LATE1_FSS = [late1, late1, late1, off]
    LATE2_MTWT = [late2, late2, late2, late2]
    LATE2_FSS = [late2, late2, late2, off]
    LATE3_MTWT = [late3, late3, late3, late3]
    LATE3_FSS = [late3, late3, late3, off]
    NIGHT1_MTWT = [night1, night1, night1, night1, off]
    NIGHT1_FSS = [off, night1, night1, night1, off, off]
    NIGHT2_MTWT = [night2, night2, night2, night2, off]
    NIGHT2_FSS = [off, night2, night2, night2, off, off]
    doctors[1].set_at(135, NIGHT1_FSS)
    doctors[2].set_at(142, NIGHT1_FSS)
    doctors[3].set_at(121, NIGHT1_FSS)
    doctors[4].set_at(58, NWD * 6)  # remove NIGHT2_FSS
    doctors[4].set_at(128, NIGHT1_FSS)
    doctors[4].set_at(177, NIGHT1_FSS)  # two NIGHT1_FSS shifts
    doctors[5].set_at(107, NIGHT1_FSS)
    doctors[6].set_at(90, NWD * 4)  # remove LATE2_MTWT
    doctors[6].set_at(97, LATE2_MTWT)
    doctors[6].set_at(114, NIGHT1_FSS)
    doctors[7].set_at(79, NWD * 6)  # remove NIGHT2_FSS
    doctors[7].set_at(90, LATE2_MTWT)
    doctors[7].set_at(97, NWD * 4)  # remove LATE2_MTWT
    doctors[7].set_at(100, NIGHT1_FSS)
    doctors[8].set_at(156, NIGHT1_FSS)
    doctors[9].set_at(79, NIGHT2_FSS)
    doctors[9].set_at(93, NWD * 6)  # remove NIGHT2_FSS
    doctors[9].set_at(163, NIGHT1_FSS)
    doctors[10].set_at(45, NWD * 4)  # remove LATE2_FSS
    doctors[10].set_at(52, LATE2_FSS)
    doctors[10].set_at(170, NIGHT1_FSS)
    doctors[11].set_at(45, LATE2_FSS)
    doctors[11].set_at(52, NWD * 4)  # remove LATE2_FSS
    # no NIGHT1_FSS for doctor 11
    doctors[11].set_at(58, NIGHT2_FSS)  # but one of these
    doctors[12].set_at(2, NIGHT1_FSS)
    doctors[13].set_at(9, NIGHT1_FSS)
    doctors[14].set_at(16, NIGHT1_FSS)
    doctors[15].set_at(23, NIGHT1_FSS)
    doctors[15].set_at(41, NWD * 4)  # remove LATE3_MTWT
    doctors[15].set_at(48, LATE3_MTWT)
    doctors[15].set_at(143, LATE2_FSS)
    doctors[16].set_at(30, NIGHT1_FSS)
    doctors[16].set_at(41, LATE3_MTWT)
    doctors[16].set_at(48, NWD * 4)  # remove LATE3_MTWT
    doctors[17].set_at(37, NIGHT1_FSS)
    doctors[18].set_at(44, NIGHT1_FSS)
    doctors[18].set_at(10, NWD * 4)  # remove LATE3_FSS
    doctors[18].set_at(17, LATE3_FSS)
    doctors[18].set_at(10, LATE3_FSS)
    doctors[18].set_at(17, NWD * 4)  # remove LATE3_FSS
    doctors[19].set_at(51, NIGHT1_FSS)
    doctors[20].set_at(58, NIGHT1_FSS)
    doctors[21].set_at(65, NIGHT1_FSS)
    doctors[21].set_at(139, NWD * 6)  # remove NIGHT1_MTWT
    doctors[21].set_at(150, LATE2_FSS)
    doctors[21].set_at(157, NWD * 4)  # remove LATE1_FSS
    doctors[21].set_at(160, NIGHT1_MTWT)
    doctors[22].set_at(72, NIGHT1_FSS)
    doctors[23].set_at(79, NIGHT1_FSS)
    doctors[23].set_at(97, NWD * 4)  # remove LATE3_MTWT, *then*...
    doctors[23].set_at(93, NIGHT2_FSS)
    doctors[23].set_at(104, LATE3_MTWT)
    doctors[24].set_at(86, NIGHT1_FSS)
    doctors[23].set_at(97, LATE3_MTWT)
    doctors[23].set_at(104, NWD * 4)  # remove LATE3_MTWT
    doctors[24].set_at(143, NWD * 4)  # remove LATE2_FSS, *then*...
    doctors[24].set_at(139, NIGHT1_MTWT)
    doctors[24].set_at(160, NWD * 4)  # remove NIGHT1_MTWT
    doctors[25].set_at(93, NIGHT1_FSS)
    doctors[25].set_at(150, NWD * 4)  # remove LATE2_FSS
    doctors[25].set_at(157, LATE1_FSS)
    doctors[26].set_at(149, NIGHT1_FSS)
    # Remove the basedoc (only present for temporary 1-based indexing!)
    doctors = doctors[1:]

    # Rota
    return Rota(
        "CPFT Aug 2015 South", shifts, doctors,
        start_date=datetime.date(2015, 8, 5),  # start on a Wednesday
        nwd_shifts=[nwd],
        comments=[
            "<b>Author:</b> CPFT Medical Staffing, Aug 2015.",
            "<b>Banding:</b> Band 1B (40%).",
            "<b>Supplement cost:</b> 14 × 0.4 = 5.6 SHOs + 12 × 0.4 = 4.8 "
            "SpRs, or approx. £360,000 pa.",
        ],
    )


def cpft_actual_aug2015_north():
    # Shifts
    # ... colours to match the original
    nwd = Shift(
        "Normal_working_day", "nwd", datetime.time(9), 8, nwd_only=True,
        resident=True, rgb=(217, 150, 148))
    late1 = Shift(
        "Late_1", "L1", datetime.time(9), 12.5,
        roles=["Late", "Late_1"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=(142, 180, 227))
    late2 = Shift(
        "Late_2", "L2", datetime.time(9), 12.5,
        roles=["Late", "Late_2"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=(0, 176, 240))
    late3 = Shift(
        "Late_3", "L3", datetime.time(9), 12.5,
        roles=["Late", "Late_3"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=(23, 55, 94))
    night1 = Shift(
        "Night_1", "N1", datetime.time(21, 15), 12,
        roles=["Night", "Night_1"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=(146, 208, 80))
    night2 = Shift(
        "Night_2", "N2", datetime.time(21, 15), 12,
        roles=["Night", "Night_2"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=(79, 98, 40))
    off = Shift(
        "Off", "OFF", datetime.time(9, 15), 7.75, work=False)
    shifts = [nwd, late1, late2, late3, night1, night2, off]

    # Doctors
    # ... so confusing to edit them manually! Just a list

    basedoc = Doctor("BASEDOC", [
        night1, night1, night1, night1, off, nwd, nwd,
        nwd, nwd, nwd, nwd, late3, late3, late3,
        off, nwd, nwd, nwd, nwd, nwd, nwd,
        late1, late1, late1, late1, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, off, night2, night2, night2,
        off, off, off, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, late1, late1, late1,
        off, nwd, nwd, nwd, nwd, nwd, nwd,
        night2, night2, night2, night2, off, nwd, nwd,
        nwd, nwd, nwd, nwd, late2, late2, late2,
        off, nwd, nwd, nwd, nwd, nwd, nwd,
        late3, late3, late3, late3, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, off, night1, night1, night1,
        off, off, off, nwd, nwd, nwd, nwd,
        late2, late2, late2, late2, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
    ], leave_weeks_per_year=0)
    basedoc.rotate(-2)  # start on a Wednesday instead of a Monday
    doctors = [basedoc]  # so we can operate with 1-indexing for a while
    for i in range(1, 18 + 1):
        if i in [2, 3, 9, 10, 15, 17]:
            name = str(i) + "_SpR"
            ooh_roles = ["Section 12 approved"]
            leave_weeks_per_year = 6
        else:
            name = str(i) + "_SHO"
            ooh_roles = []
            leave_weeks_per_year = 5
        name
        d = basedoc.copy(name=name, rotate=(i - 1) * 7)
        d.ooh_roles = ooh_roles
        d.leave_weeks_per_year = leave_weeks_per_year
        doctors.append(d)
    # Now the deviations from a simple pattern
    NWD = [nwd]
    LATE1_MTWT = [late1, late1, late1, late1]
    LATE1_FSS = [late1, late1, late1, off]
    LATE2_MTWT = [late2, late2, late2, late2]
    LATE2_FSS = [late2, late2, late2, off]
    LATE3_MTWT = [late3, late3, late3, late3]
    LATE3_FSS = [late3, late3, late3, off]
    NIGHT1_MTWT = [night1, night1, night1, night1, off]
    NIGHT1_FSS = [off, night1, night1, night1, off, off]
    NIGHT2_MTWT = [night2, night2, night2, night2, off]
    NIGHT2_FSS = [off, night2, night2, night2, off, off]
    doctors[4].set_at(52, LATE2_FSS)
    doctors[16].set_at(52, NWD * 4)  # remove LATE2_FSS
    doctors[4].set_at(73, NWD * 4)  # remove LATE1_FSS
    doctors[16].set_at(73, LATE1_FSS)
    # Remove the basedoc (only present for temporary 1-based indexing!)
    doctors = doctors[1:]

    # Rota
    return Rota(
        "CPFT Aug 2015 North", shifts, doctors,
        start_date=datetime.date(2015, 8, 5),  # start on a Wednesday
        nwd_shifts=[nwd],
        comments=[
            "<b>Author:</b> CPFT Medical Staffing, Aug 2015.",
            "<b>Banding:</b> Band 2B (50%). "
            "<b><i>Band 2 rotas are not EWTD compliant.</i></b>",
            "<b>Supplement cost:</b> 12 × 0.5 = 6 SHOs + 6 × 0.5 = 3 SpRs, or "
            "approx. £305,000.",
        ],
    )


# -----------------------------------------------------------------------------
# CPFT draft
# -----------------------------------------------------------------------------

def cpft_draft_1_south():
    # Shifts
    nwd = Shift(
        "Normal_working_day", "nwd", datetime.time(9), 8, nwd_only=True,
        resident=True, rgb=COLOURS.NWD)
    sho_late1 = Shift(
        "SHO_Late_1_F", "F", datetime.time(9), 12.5,
        roles=["Fulbourn SHO"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=COLOURS.LATE_F)
    sho_late2 = Shift(
        "SHO_Late_2_A", "A", datetime.time(9), 12.5,
        roles=["Addenbrooke’s SHO"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=COLOURS.LATE_A)
    sho_night = Shift(
        "SHO_Night", "N", datetime.time(21), 12.25,
        roles=["Fulbourn SHO", "Addenbrooke’s SHO"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=COLOURS.NIGHT_FA)
    spr_late = Shift(
        "SpR_Late", "L", datetime.time(9), 12.5,
        roles=["Section 12 approved"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=COLOURS.SPR_LATE)
    spr_night = Shift(
        "SpR_Night", "N", datetime.time(21), 12.25,
        roles=["Section 12 approved"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=COLOURS.SPR_NIGHT)
    off = Shift(
        "Off", "OFF", datetime.time(9, 15), 7.75, work=False,
        rgb=COLOURS.OFF)
    shifts = [nwd, sho_late1, sho_late2, sho_night, spr_late, spr_night, off]

    # Doctors
    n_sho = 14
    base_sho = Doctor("Prototype SHO", [
        sho_night, sho_night, sho_night, sho_night, off, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        sho_late1, sho_late1, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, sho_late1, sho_late1, nwd, nwd, nwd,
        nwd, nwd, nwd, off, sho_late1, sho_late1, sho_late1,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        sho_late2, sho_late2, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, sho_late2, sho_late2, nwd, nwd, nwd,
        nwd, nwd, nwd, off, sho_late2, sho_late2, sho_late2,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, off, sho_night, sho_night, sho_night,
        off, off, nwd, nwd, nwd, nwd, nwd,
    ], leave_weeks_per_year=5)
    base_sho.pad(n_sho * 7, [nwd])
    base_sho.rotate(-2)  # from a Monday start to a Wednesday start
    shos = []
    for i in range(n_sho):
        shos.append(base_sho.copy("SHO_" + str(i + 1), rotate=i*7))

    n_spr = 12
    # Although I like weeks of nights, the RCP strongly discourages them.
    # So this is a 4-and-3. We can follow the classic RCP (2006) Table 2,
    # then remove some 'off' shifts, as it seems we can.
    base_spr = Doctor("Prototype SpR", [
        # Mon ... Sun
        spr_night, spr_night, spr_night, spr_night, off, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, spr_late, spr_late, spr_late,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, spr_late, spr_late, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        spr_late, spr_late, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, spr_night, spr_night, spr_night,
        off, nwd, nwd, nwd, nwd, nwd, nwd,
    ], leave_weeks_per_year=6)
    base_spr.pad(n_spr * 7, [nwd])
    base_spr.rotate(-2)  # from a Monday start to a Wednesday start
    sprs = []
    for i in range(n_spr):
        sprs.append(base_spr.copy("SpR_" + str(i + 1), rotate=i*7))

    doctors = shos + sprs

    # Rota
    return Rota(
        "CPFT draft 1, South", shifts, doctors,
        start_date=datetime.date(2015, 8, 5),  # Wednesday
        nwd_shifts=[nwd],
        prototypes=[base_sho, base_spr],
        comments=[
            "<b>Synopsis:</b> an improvement to the South rota by splitting "
            "SHOs and SpRs; however, it doesn’t help the North problem. "
            "See draft 2 instead.",
            "<b>Author:</b> Rudolf Cardinal, Aug 2015.",
            "<b>Banding:</b> estimated at 1B (40%).",
            "<b>Supplement cost:</b> 14 × 0.4 = 5.6 SHOs + 12 × 0.4 = 4.8 "
            "SpRs, or approx. £360,000 pa.",
        ],
    )


def cpft_draft_2_combined():
    # Shifts
    nwd = Shift(
        "Normal_working_day", "nwd", datetime.time(9), 8, nwd_only=True,
        resident=True, rgb=COLOURS.NWD)
    s_sho_late1 = Shift(
        "S_SHO_Late_1_F", "F", datetime.time(9), 12.5,
        roles=["Fulbourn SHO"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=COLOURS.LATE_F)
    s_sho_late2 = Shift(
        "S_SHO_Late_2_A", "A", datetime.time(9), 12.5,
        roles=["Addenbrooke’s SHO"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=COLOURS.LATE_A)
    s_sho_night = Shift(
        "S_SHO_Night", "N", datetime.time(21), 12.25,
        roles=["Fulbourn SHO", "Addenbrooke’s SHO"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=COLOURS.NIGHT_FA)
    n_sho_late1 = Shift(
        "N_SHO_Late_1_CC", "C", datetime.time(9), 12.5,
        roles=["Cavell Centre SHO"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=COLOURS.LATE_C)
    n_sho_late2 = Shift(
        "N_SHO_Late_2_PCH", "P", datetime.time(9), 12.5,
        roles=["PCH SHO"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=COLOURS.LATE_P)
    n_sho_night = Shift(
        "N_SHO_Night", "N", datetime.time(21), 12.25,
        roles=["Cavell Centre SHO", "PCH SHO"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=COLOURS.NIGHT_CP)
    spr_late = Shift(
        "SpR_Late", "L", datetime.time(9), 12.5,
        roles=["Section 12 approved"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=COLOURS.SPR_LATE)
    spr_night = Shift(
        "SpR_Night", "N", datetime.time(21), 12.25,
        roles=["Section 12 approved"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=COLOURS.SPR_NIGHT)
    off = Shift(
        "Off", "OFF", datetime.time(9, 15), 7.75, work=False)
    shifts = [nwd,
              s_sho_late1, s_sho_late2, s_sho_night,
              n_sho_late1, n_sho_late2, n_sho_night,
              spr_late, spr_night,
              off]

    # Doctors
    n_south_sho = 14
    south_base_sho = Doctor("South prototype SHO", [
        s_sho_night, s_sho_night, s_sho_night, s_sho_night, off, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        s_sho_late1, s_sho_late1, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, s_sho_late1, s_sho_late1, nwd, nwd, nwd,
        nwd, nwd, nwd, off, s_sho_late1, s_sho_late1, s_sho_late1,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        s_sho_late2, s_sho_late2, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, s_sho_late2, s_sho_late2, nwd, nwd, nwd,
        nwd, nwd, nwd, off, s_sho_late2, s_sho_late2, s_sho_late2,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, off, s_sho_night, s_sho_night, s_sho_night,
        off, off, nwd, nwd, nwd, nwd, nwd,
    ], leave_weeks_per_year=5)
    south_base_sho.pad(n_south_sho * 7, [nwd])
    south_base_sho.rotate(-2)  # from a Monday start to a Wednesday start
    south_shos = []
    for i in range(n_south_sho):
        south_shos.append(
            south_base_sho.copy("S" + str(i + 1), rotate=i*7))

    # Doctors
    n_north_sho = 12
    north_base_sho = Doctor("North prototype SHO", [
        # The same pattern as the South fails (i.e. yields 2B) because there
        # are slightly fewer doctors.
        n_sho_night, n_sho_night, n_sho_night, n_sho_night, off, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        n_sho_late1, n_sho_late1, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, n_sho_late1, n_sho_late1, off, nwd, nwd,  # INSERT OFF
        nwd, nwd, nwd, off, n_sho_late1, n_sho_late1, n_sho_late1,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        n_sho_late2, n_sho_late2, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, n_sho_late2, n_sho_late2, off, nwd, nwd,  # INSERT OFF
        nwd, nwd, nwd, off, n_sho_late2, n_sho_late2, n_sho_late2,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, off, n_sho_night, n_sho_night, n_sho_night,
        off, off, nwd, nwd, nwd, nwd, nwd,
    ], leave_weeks_per_year=5)
    north_base_sho.pad(n_north_sho * 7, [nwd])
    north_base_sho.rotate(-2)  # from a Monday start to a Wednesday start
    north_shos = []
    for i in range(n_north_sho):
        north_shos.append(
            north_base_sho.copy("N" + str(i + 1), rotate=i*7))

    n_spr = 18  # 18 normal; 12 works at 40%; 11 doesn't quite
    #
    base_spr = Doctor("Prototype SpR", [
        # Mon ... Sun
        spr_night, spr_night, spr_night, spr_night, off, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, spr_late, spr_late, spr_late,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, spr_late, spr_late, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        spr_late, spr_late, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, spr_night, spr_night, spr_night,
        off, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
    ], leave_weeks_per_year=6)
    base_spr.pad(n_spr * 7, [nwd])
    base_spr.rotate(-2)  # from a Monday start to a Wednesday start
    sprs = []
    for i in range(n_spr):
        sprs.append(base_spr.copy("R" + str(i + 1), rotate=i*7))

    doctors = north_shos + sprs + south_shos

    return Rota(
        "CPFT draft 2, combined North/South", shifts, doctors,
        start_date=datetime.date(2015, 8, 5),  # Wednesday
        nwd_shifts=[nwd],
        prototypes=[south_base_sho, north_base_sho, base_spr],
        comments=[
            "<b>Synopsis:</b> combine SpRs across North/South.",
            "<b>Author:</b> Rudolf Cardinal, Aug 2015.",
            "<b>Banding:</b> estimated at 1B (40%) for South SHOs and all "
            "SpRs, and 1A (50%) for North SHOs.",
            "<b>Supplement cost:</b> 14 × 0.4 + 12 × 0.5 = 11.6 SHOs "
            "+ 18 × 0.4 = 7.2 SpRs, or approx. £642,000 pa. "
            "(About £22,000 cheaper than the existing, the saving being from "
            "6 SpRs moving from 50% to 40%.)",
            """
            <b>Rationale:</b>
            <ul>
                <li>The previous rota (to Aug 2015) failed in part because
                there was too much work for the overnight SpRs for a
                non-resident on-call rota. So full shifts ± resident rotas
                are the minimum.</li>

                <li>The Aug 2015 rota gives less than 100% Tier 2/
                section 12-approved doctor coverage (with two SHOs on as the
                two night doctors in each region reasonably often, particularly
                in the North, and two SpRs sometimes).
                It confuses switchboards by failing to distinguish well
                between SHOs (CT1–3, Tier 1) and SpRs (ST4–6, Tier 2).
                These are obvious reasons to split the SpR rota from the SHO
                rota.</li>

                <li>It’s easy to improve the South (Cambridge) rota while
                keeping it at 40% banding. However, it’s impossible to split
                SHOs/SpRs in the North, because a 1:6 rota for SpRs isn’t
                achievable in a legitimate banding. Moreover, with the closure
                of the North s136 suite, there is a training gap for SpRs (in
                terms of emergency MHA work). The obvious solution to this is
                to combine the North/South SpR rotas.</li>
            </ul>
            """,
            "<b>Doctor coding:</b> N = North SHO, S = South SHO, "
            "R = registrar.",
            """
            <b>Advantages:</b>
            <ul>
                <li><b>Cheaper,</b> by about £22k/y.</li>

                <li><b>Removes Band 2,</b> which we thought had been disallowed
                since 2009 at the end of the phased requirement for EWTD
                compliance.</li>

                <li><b>100% Tier 2/section 12 coverage.</b></li>

                <li><b>Better normal working day coverage, on average.</b>
                Specifically:
                <ul>
                    <li>North SHOs improve from 54–61% to 60–64%.</li>
                    <li>South SHOs worsen from 71–76% to 68–72%.</li>
                    <li>SpRs improve from 56–60% (North) and 71–74% (South)
                    to 83–91%.</li>
                </ul>
                This has corresponding effects for doctors who work less than
                full time for the NHS (LTFT trainees, academics); the better
                the NWD coverage, the less their other job/activity is
                detrimentally impacted.</li>

                <li><b>Resumption of s136 MHA experience for North SpRs,</b>
                with a corresponding 33% reduction in s136 work for South
                SpRs.</li>

                <li><b>Less confusing for switchboards,</b>
                with a clear SHO/SpR split.</li>

                <li>Swaps are easier for SpRs (more SpRs to swap with).</li>

                <li>Can cope with a reduction in SpRs (if some leave for
                consultant posts), down to a minimum of 12 whilst maintaining
                40% banding.</li>

                <li>Handover time increased by 50%.</li>
            </ul>
            """,
            """
            <b>Disadvantages:</b>
            <ul>
                <li><b>More travel for SpRs.</b></li>

                <li><b>One SpR at night rather than two.</b> This shouldn’t be
                a problem for Mental Health Act Assessments, since there’s only
                one AMHP on in any case, so this should only really be relevant
                if both SHOs are swamped in their respective A&amp;E
                departments.</li>

                <li><b>SpRs not dedicated to a single
                location/A&amp;E.</b></li>

                <li>Fractionally (2%) longer night shifts.</li>

                <li>The concept of a ‘resident’ SpR covering two areas
                is a little unusual, though all it really means is ‘available
                for duty and expected to be working or ready to work at all
                times’. In any case, the ‘resident’ test only applies to
                true on-call rotas, not full shift rotas; this is a full
                shift rota, so we can say that the SpR should be working or
                ready to work at all times, with the expectation of normal
                breaks but not of sleep, but where they are is a matter for
                clinical judgement.</li>
            </ul>
            """,
            """
            <b>Unclear:
            <ul>
                <li>How is Hinchingbrooke covered at present?</li>
            </ul>
            </b>
            """,
        ],
    )


def cpft_draft_3_split():
    # Shifts
    nwd = Shift(
        "Normal_working_day", "nwd", datetime.time(9), 8, nwd_only=True,
        resident=True, rgb=COLOURS.NWD)
    s_sho_late1 = Shift(
        "S_SHO_Late_1_F", "F", datetime.time(9), 12.5,
        roles=["Fulbourn SHO"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=COLOURS.LATE_F)
    s_sho_late2 = Shift(
        "S_SHO_Late_2_A", "A", datetime.time(9), 12.5,
        roles=["Addenbrooke’s SHO"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=COLOURS.LATE_A)
    s_sho_night = Shift(
        "S_SHO_Night", "N", datetime.time(21), 12.25,
        roles=["Fulbourn SHO", "Addenbrooke’s SHO"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=COLOURS.NIGHT_FA)
    s_spr_night = Shift(
        "S_SpR_Night", "SN", datetime.time(21), 12.25,
        roles=["South SpR", "Section 12 approved"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=COLOURS.S_SPR_NIGHT)

    n_sho_late1 = Shift(
        "N_SHO_Late_1_CC", "C", datetime.time(9), 12.5,
        roles=["Cavell Centre SHO"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=COLOURS.LATE_C)
    n_sho_late2 = Shift(
        "N_SHO_Late_2_PCH", "P", datetime.time(9), 12.5,
        roles=["PCH SHO"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=COLOURS.LATE_P)
    n_sho_night = Shift(
        "N_SHO_Night", "N", datetime.time(21), 12.25,
        roles=["Cavell Centre SHO", "PCH SHO"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=COLOURS.NIGHT_CP)
    n_spr_night = Shift(
        "N_SpR_Night", "NN", datetime.time(21), 12.25,
        roles=["North SpR", "Section 12 approved"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=COLOURS.N_SPR_NIGHT)

    spr_late = Shift(
        "SpR_Late", "L", datetime.time(9), 12.5,
        roles=["North SpR", "South SpR", "Section 12 approved"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=COLOURS.SPR_LATE)

    off = Shift(
        "Off", "OFF", datetime.time(9, 15), 7.75, work=False)
    shifts = [nwd,
              s_sho_late1, s_sho_late2, s_sho_night, s_spr_night,
              n_sho_late1, n_sho_late2, n_sho_night, n_spr_night,
              spr_late,
              off]

    # Doctors
    n_south_sho = 14
    south_base_sho = Doctor("South prototype SHO", [
        s_sho_night, s_sho_night, s_sho_night, s_sho_night, off, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        s_sho_late1, s_sho_late1, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, s_sho_late1, s_sho_late1, nwd, nwd, nwd,
        nwd, nwd, nwd, off, s_sho_late1, s_sho_late1, s_sho_late1,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        s_sho_late2, s_sho_late2, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, s_sho_late2, s_sho_late2, nwd, nwd, nwd,
        nwd, nwd, nwd, off, s_sho_late2, s_sho_late2, s_sho_late2,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, off, s_sho_night, s_sho_night, s_sho_night,
        off, off, nwd, nwd, nwd, nwd, nwd,
    ], leave_weeks_per_year=5)
    south_base_sho.pad(n_south_sho * 7, [nwd])
    south_base_sho.rotate(-2)  # from a Monday start to a Wednesday start
    south_shos = []
    for i in range(n_south_sho):
        south_shos.append(
            south_base_sho.copy("S" + str(i + 1), rotate=i*7))

    # Doctors
    n_north_sho = 12
    north_base_sho = Doctor("North prototype SHO", [
        # The same pattern as the South fails (i.e. yields 2B) because there
        # are slightly fewer doctors.
        n_sho_night, n_sho_night, n_sho_night, n_sho_night, off, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        n_sho_late1, n_sho_late1, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, n_sho_late1, n_sho_late1, off, nwd, nwd,  # INSERT OFF
        nwd, nwd, nwd, off, n_sho_late1, n_sho_late1, n_sho_late1,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        n_sho_late2, n_sho_late2, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, n_sho_late2, n_sho_late2, off, nwd, nwd,  # INSERT OFF
        nwd, nwd, nwd, off, n_sho_late2, n_sho_late2, n_sho_late2,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, off, n_sho_night, n_sho_night, n_sho_night,
        off, off, nwd, nwd, nwd, nwd, nwd,
    ], leave_weeks_per_year=5)
    north_base_sho.pad(n_north_sho * 7, [nwd])
    north_base_sho.rotate(-2)  # from a Monday start to a Wednesday start
    north_shos = []
    for i in range(n_north_sho):
        north_shos.append(
            north_base_sho.copy("N" + str(i + 1), rotate=i*7))

    n_north_spr = 9
    north_base_spr = Doctor("North prototype SpR", [
        # Mon ... Sun

        # block with no lates:
        n_spr_night, n_spr_night, n_spr_night, n_spr_night, off, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, n_spr_night, n_spr_night, n_spr_night,
        off, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,

        # block with lates:
        n_spr_night, n_spr_night, n_spr_night, n_spr_night, off, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, spr_late, spr_late, spr_late,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, n_spr_night, n_spr_night, n_spr_night,
        off, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, spr_late, spr_late, off, nwd, nwd,  # INSERT OFF
        spr_late, spr_late, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
    ], leave_weeks_per_year=6)
    north_base_spr.pad(n_north_spr * 7, [nwd])
    north_base_spr.rotate(-2)  # from a Monday start to a Wednesday start
    north_sprs = []
    for i in range(n_north_spr):
        north_sprs.append(north_base_spr.copy("NR" + str(i + 1), rotate=i*7))

    n_south_spr = 9
    south_base_spr = Doctor("South prototype SpR", [
        # Mon ... Sun
        # block with lates:
        s_spr_night, s_spr_night, s_spr_night, s_spr_night, off, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, spr_late, spr_late, spr_late,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, s_spr_night, s_spr_night, s_spr_night,
        off, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, spr_late, spr_late, off, nwd, nwd,  # INSERT OFF
        spr_late, spr_late, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,

        # block with no lates:
        s_spr_night, s_spr_night, s_spr_night, s_spr_night, off, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, s_spr_night, s_spr_night, s_spr_night,
        off, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
    ], leave_weeks_per_year=6)
    south_base_spr.pad(n_south_spr * 7, [nwd])
    south_base_spr.rotate(-2)  # from a Monday start to a Wednesday start
    south_sprs = []
    for i in range(n_south_spr):
        south_sprs.append(south_base_spr.copy("SR" + str(i + 1), rotate=i*7))

    doctors = north_shos + north_sprs + south_sprs + south_shos

    return Rota(
        "CPFT draft 3, combined North/South SpRs in evenings, split "
        "North/South at night", shifts, doctors,
        start_date=datetime.date(2015, 8, 5),  # Wednesday
        nwd_shifts=[nwd],
        prototypes=[south_base_sho, north_base_sho,
                    north_base_spr, south_base_spr],
        comments=[
            "<b>Synopsis:</b> North and South SpRs are split at night (two 1:9"
            " rotas), but there is a combined (1:18) evening SpR rota."
            " A single SpR covers both halves of the county in the evening.",
            "<b>Author:</b> Rudolf Cardinal, Aug 2015.",
            "<b>Banding:</b> estimated at 1B (40%) for South SHOs and all "
            "SpRs, and 1A (50%) for North SHOs.",
            "The starting point is draft 2 <b>(q.v.)</b>.",
            "1:9 North SpRs and 1:9 South SpRs for both evenings and nights "
            "fails (Band 2)",
            "If a North and a South SpR are required at night, but only one "
            "SpR is required in the evening (when 4 SHOs are on across the "
            "county, this pattern works at Band 1.",
            "<b>Three SpRs who are South by day will be North at night, "
            "and considered North in this pattern.</b>",
            "Section 136 North/South inequality persists, but is not "
            "complete, since North SpRs regularly cover Cambridge s136 work "
            "for the evening shifts.",
            "Remains Band 1B with 8+8 SpRs; not with 7+7.",
        ],
    )


def cpft_draft_4_split():
    # Shifts
    nwd = Shift(
        "Normal_working_day", "nwd", datetime.time(9), 8, nwd_only=True,
        resident=True, rgb=COLOURS.NWD)
    s_sho_late1 = Shift(
        "S_SHO_Late_1_F", "F", datetime.time(9), 12.5,
        roles=["Fulbourn SHO"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=COLOURS.LATE_F)
    s_sho_late2 = Shift(
        "S_SHO_Late_2_A", "A", datetime.time(9), 12.5,
        roles=["Addenbrooke’s SHO"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=COLOURS.LATE_A)
    s_sho_night = Shift(
        "S_SHO_Night", "N", datetime.time(21), 12.25,
        roles=["Fulbourn SHO", "Addenbrooke’s SHO"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=COLOURS.NIGHT_FA)
    s_spr_late = Shift(
        "S_SpR_Late", "SL", datetime.time(9), 12.5,
        roles=["South SpR", "Section 12 approved"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=COLOURS.S_SPR_LATE)

    n_sho_late1 = Shift(
        "N_SHO_Late_1_CC", "C", datetime.time(9), 12.5,
        roles=["Cavell Centre SHO"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=COLOURS.LATE_C)
    n_sho_late2 = Shift(
        "N_SHO_Late_2_PCH", "P", datetime.time(9), 12.5,
        roles=["PCH SHO"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=COLOURS.LATE_P)
    n_sho_night = Shift(
        "N_SHO_Night", "N", datetime.time(21), 12.25,
        roles=["Cavell Centre SHO", "PCH SHO"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=COLOURS.NIGHT_CP)
    n_spr_late = Shift(
        "N_SpR_Late", "NL", datetime.time(9), 12.5,
        roles=["North SpR", "Section 12 approved"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=COLOURS.N_SPR_LATE)

    spr_night = Shift(
        "SpR_Night", "N", datetime.time(21), 12.25,
        roles=["North SpR", "South SpR", "Section 12 approved"],
        resident=True, shift_type=SHIFT_TYPES.FULL, rgb=COLOURS.SPR_NIGHT)


    off = Shift(
        "Off", "OFF", datetime.time(9, 15), 7.75, work=False)
    shifts = [nwd,
              s_sho_late1, s_sho_late2, s_sho_night, s_spr_late,
              n_sho_late1, n_sho_late2, n_sho_night, n_spr_late,
              spr_night,
              off]

    # Doctors
    n_south_sho = 14
    south_base_sho = Doctor("South prototype SHO", [
        s_sho_night, s_sho_night, s_sho_night, s_sho_night, off, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        s_sho_late1, s_sho_late1, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, s_sho_late1, s_sho_late1, nwd, nwd, nwd,
        nwd, nwd, nwd, off, s_sho_late1, s_sho_late1, s_sho_late1,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        s_sho_late2, s_sho_late2, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, s_sho_late2, s_sho_late2, nwd, nwd, nwd,
        nwd, nwd, nwd, off, s_sho_late2, s_sho_late2, s_sho_late2,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, off, s_sho_night, s_sho_night, s_sho_night,
        off, off, nwd, nwd, nwd, nwd, nwd,
    ], leave_weeks_per_year=5)
    south_base_sho.pad(n_south_sho * 7, [nwd])
    south_base_sho.rotate(-2)  # from a Monday start to a Wednesday start
    south_shos = []
    for i in range(n_south_sho):
        south_shos.append(
            south_base_sho.copy("S" + str(i + 1), rotate=i*7))

    # Doctors
    n_north_sho = 12
    north_base_sho = Doctor("North prototype SHO", [
        # The same pattern as the South fails (i.e. yields 2B) because there
        # are slightly fewer doctors.
        n_sho_night, n_sho_night, n_sho_night, n_sho_night, off, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        n_sho_late1, n_sho_late1, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, n_sho_late1, n_sho_late1, off, nwd, nwd,  # INSERT OFF
        nwd, nwd, nwd, off, n_sho_late1, n_sho_late1, n_sho_late1,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        n_sho_late2, n_sho_late2, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, n_sho_late2, n_sho_late2, off, nwd, nwd,  # INSERT OFF
        nwd, nwd, nwd, off, n_sho_late2, n_sho_late2, n_sho_late2,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, off, n_sho_night, n_sho_night, n_sho_night,
        off, off, nwd, nwd, nwd, nwd, nwd,
    ], leave_weeks_per_year=5)
    north_base_sho.pad(n_north_sho * 7, [nwd])
    north_base_sho.rotate(-2)  # from a Monday start to a Wednesday start
    north_shos = []
    for i in range(n_north_sho):
        north_shos.append(
            north_base_sho.copy("N" + str(i + 1), rotate=i*7))

    n_north_spr = 9
    north_base_spr = Doctor("North prototype SpR", [
        # Mon ... Sun

        # block with no nights:
        nwd, nwd, n_spr_late, n_spr_late, off, nwd, nwd,  # INSERT OFF
        n_spr_late, n_spr_late, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, n_spr_late, n_spr_late, n_spr_late,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,

        # block with nights:
        nwd, nwd, n_spr_late, n_spr_late, off, nwd, nwd,  # INSERT OFF
        n_spr_late, n_spr_late, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        spr_night, spr_night, spr_night, spr_night, off, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, n_spr_late, n_spr_late, n_spr_late,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, spr_night, spr_night, spr_night,
        off, nwd, nwd, nwd, nwd, nwd, nwd,
    ], leave_weeks_per_year=6)
    north_base_spr.pad(n_north_spr * 7, [nwd])
    north_base_spr.rotate(-2)  # from a Monday start to a Wednesday start
    north_sprs = []
    for i in range(n_north_spr):
        north_sprs.append(north_base_spr.copy("NR" + str(i + 1), rotate=i*7))

    n_south_spr = 9
    south_base_spr = Doctor("South prototype SpR", [
        # Mon ... Sun
        # block with nights:
        nwd, nwd, s_spr_late, s_spr_late, off, nwd, nwd,  # INSERT OFF
        s_spr_late, s_spr_late, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        spr_night, spr_night, spr_night, spr_night, off, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, s_spr_late, s_spr_late, s_spr_late,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, spr_night, spr_night, spr_night,
        off, nwd, nwd, nwd, nwd, nwd, nwd,

        # block with no nights:
        nwd, nwd, s_spr_late, s_spr_late, off, nwd, nwd,  # INSERT OFF
        s_spr_late, s_spr_late, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, s_spr_late, s_spr_late, s_spr_late,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
        nwd, nwd, nwd, nwd, nwd, nwd, nwd,
    ], leave_weeks_per_year=6)
    south_base_spr.pad(n_south_spr * 7, [nwd])
    south_base_spr.rotate(-2)  # from a Monday start to a Wednesday start
    south_sprs = []
    for i in range(n_south_spr):
        south_sprs.append(south_base_spr.copy("SR" + str(i + 1), rotate=i*7))

    doctors = north_shos + north_sprs + south_sprs + south_shos

    return Rota(
        "CPFT draft 4, split North/South SpRs in evenings, combined "
        "North/South SpRs at night", shifts, doctors,
        start_date=datetime.date(2015, 8, 5),  # Wednesday
        nwd_shifts=[nwd],
        prototypes=[south_base_sho, north_base_sho,
                    north_base_spr, south_base_spr],
        comments=[
            "<b>Synopsis:</b> North and South SpRs are split in the evening "
            "(two 1:9 rotas), but combined (1:18) at night.",
            "<b>Author:</b> Rudolf Cardinal, Aug 2015.",
            "<b>Banding:</b> estimated at 1B (40%) for South SHOs and all "
            "SpRs, and 1A (50%) for North SHOs.",
            "<b>See also drafts 2/3.</b>",
            "<b>Three SpRs who are South by day will be North at night, "
            "and considered North in this pattern.</b>",
        ],
    )


# =============================================================================
# Rota map
# =============================================================================

ROTA_GENERATORS = OrderedDict([
    ('test_prospective_cover', test_prospective_cover),
    ('cpft_actual_aug2015_south', cpft_actual_aug2015_south),
    ('cpft_actual_aug2015_north', cpft_actual_aug2015_north),
    ('cpft_draft_1_south', cpft_draft_1_south),
    ('cpft_draft_2_combined', cpft_draft_2_combined),
    ('cpft_draft_3_split', cpft_draft_3_split),
    ('cpft_draft_4_split', cpft_draft_4_split),
])


def process_rota(rotaname, daynums=False, skipanalytics=False,
                 noprototypes=False):
    if rotaname not in ROTA_GENERATORS:
        raise ValueError("Invalid rota: " + rotaname)
    fn = ROTA_GENERATORS[rotaname]
    rota = fn()
    filename = rotaname + ".html"
    rota.print_html(filename, daynums, skipanalytics, noprototypes)


# =============================================================================
# Main
# =============================================================================

def main():
    description = """
Medical rota checker. By Rudolf Cardinal (rudolf@pobox.com), Aug 2015.
Version {}
""".format(VERSION)
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "rota", nargs="*",
        help=("Rota name(s). Options are: "
              + ", ".join(ROTA_GENERATORS.keys()) + "."))
    parser.add_argument(
        "-a", "--all", action="store_true",
        help="Process all known rotas")
    parser.add_argument(
        "-d", "--daynums", action="store_true",
        help="Show day numbers as well as dates")
    parser.add_argument(
        "-s", "--skipanalytics", action="store_true",
        help="Skip analyses")
    parser.add_argument(
        "-n", "--noprototypes", action="store_true",
        help="Skip prototype patterns")
    args = parser.parse_args()

    if args.all:
        for rotaname in ROTA_GENERATORS.keys():
            process_rota(rotaname, args.daynums, args.skipanalytics,
                         args.noprototypes)
    else:
        if not args.rota:
            parser.print_help()
            sys.exit(1)
        for rotaname in args.rota:
            process_rota(rotaname, args.daynums, args.skipanalytics,
                         args.noprototypes)


# =============================================================================
# Command-line entry point
# =============================================================================

if __name__ == '__main__':
    main()
    # test()
