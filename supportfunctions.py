"""Copyright [2014] [Athindran Ramesh Kumar]

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
"""

from math import atan2, cos, floor, pi, sin, sqrt
import numpy as np

def four_digit_year(year):
  if year <= 79:
    year += 2000
  elif year > 79:
    year += 1900
  return year


def leapseconds(date):
  leaptimes = np.array([[30, 6, 2012], [31, 12, 2008], [31, 12, 2005],
                        [31, 12, 1998], [30, 6, 1997], [31, 12, 1995],
                        [30, 6, 1994], [30, 6, 1993], [30, 6, 1992],
                        [31, 12, 1990], [31, 12, 1989], [31, 12, 1987],
                        [30, 6, 1985], [30, 6, 1983], [30, 6, 1982],
                        [30, 6, 1981]])
  for tind, t in enumerate(leaptimes):
    diff = date - t
    if (diff[2] > 0 or
        (diff[2] == 0 and diff[1] > 0) or
        (diff[2] == 0 and diff[1] == 0 and diff[0] > 0)):
      return 16 - tind

def date2gps(datelist, leaps):
  year, month, day, hour, minute, second = datelist
  if month <= 2:
    year += -1
    month += 12
  yl = floor(year / 100)
  jdn = ((2 - yl + floor(yl / 4)) + floor(365.25 * year)
         + floor(30.6001 * (month + 1)) + day + 1720994.5)
  al = floor(jdn + 0.5)
  nl = (al - 7 * floor(al / 7) + 1) % 7
  secondspastmidnightgmt = nl * 86400 + hour * 3600 + minute * 60 + second
  gpstime = secondspastmidnightgmt + leaps
  return gpstime



