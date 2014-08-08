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

import numpy as np
import re
import supportfunctions as gpt

fl = re.compile('[-+]?(\d+(\.\d*)?|\.\d+)([e][-+]?\d+)?')
global cmapping
cmapping = {
    'G': 'GPS',
    'R': 'GLONASS',
    'E': 'Galileo',
    'C': 'BeiDou',
    'J': 'QZSS',
    'S': 'SBSS'
}


def parserinexnav(filename, constellation):
  f = open(filename, 'r')
  header = 0

  mappingo = np.array([
      ['--', 'a_f0', 'a_f1', 'a_f2'],
      ['IODE', 'C_rs', 'deltan', 'M_0'],
      ['C_uc', 'e', 'C_us', 'roota'],
      ['t_oe', 'C_ic', 'omega_0', 'C_is'],
      ['i_0', 'C_rc', 'omega', 'Omega_dot'],
      ['iDot', 'L2code', 'Weekno', 'L2P'],
      ['acc', 'health', 'T_GD', 'IODC'],
      ['transmittime', '--', '--', '--']])

  mappingl = (
      np.array(
          [['--', 'TauN', 'GammaN', 'tk'],
           ['X', 'VX', 'AX', 'health'], ['Y', 'VY', 'AY', 'e'],
           ['Z', 'VZ', 'AZ', '--']]))

  while (header == 0):
    linedata = f.readline()
    if('RINEX VERSION / TYPE' in linedata):
      version = float(linedata[0:8])
    if('END OF HEADER' in linedata):
      header = 1
  ephem = [{} for t in xrange(86400)]

  index = 0
  linedata = f.readline()
  while(len(linedata) != 0):
    fail = 1
    while(fail):
      if(linedata.strip() != '' and 'COMMENT' not in linedata):
        linedata = linedata.replace('D', 'e')
        sysid = linedata[0]
        if(sysid == 'R'):
          nlines = 3
          mapping = mappingl
          linedata = linedata[1:]
        elif(sysid == 'G' or sysid == 'E' or sysid == 'C' or sysid == 'J' or sysid == 'S'):
          nlines = 7
          mapping = mappingo
          linedata = linedata[1:]
        elif(filename[-1] == 'g' and (sysid == ' ' or sysid.isdigit())):
          nlines = 3
          sysid = 'R'
          mapping = mappingl
        else:
          sysid = 'G'
          nlines = 7
          mapping = mappingo
        headerdata = fl.finditer(linedata)

        for data in headerdata:
          if data.start() < 2:
            svprn = float(data.group())
          elif data.start() >= 2 and data.start() < 6:
            year = float(data.group())
          elif data.start() >= 6 and data.start() < 9:
            month = float(data.group())
            fail = 0
          elif data.start() >= 9 and data.start() < 12:
            day = float(data.group())
            fail = 0
          elif data.start() >= 12 and data.start() < 15:
            hour = float(data.group())
            fail = 0
          elif data.start() >= 15 and data.start() < 18:
            minute = float(data.group())
            fail = 0
          elif data.start() >= 18 and data.start() < 22:
            second = float(data.group())
            fail = 0
          elif data.start() >= 22 and data.start() < 41:
            value0 = float(data.group())
          elif data.start() >= 41 and data.start() < 60:
            value1 = float(data.group())
          elif data.start() >= 60:
            value2 = float(data.group())
      if(fail == 1):
        linedata = f.readline()
        if(len(linedata) == 0):
          break
      ephem[index]['constellation'] = sysid
    if(version < 3):
      year = gpt.four_digit_year(year)
    toc = gpt.date2gps([year, month, day, hour, minute, second], 0)
    ephem[index]['PRN'] = svprn
    ephem[index]['t_oc'] = toc

    ephem[index][mapping[0][1]] = value0
    ephem[index][mapping[0][2]] = value1
    ephem[index][mapping[0][3]] = value2

    if(constellation[cmapping[sysid]]['Enabled']):
      for t in xrange(nlines):
        lineread = f.readline()
        fail = 1
        while(fail):
          if(lineread.strip() != '' and ('COMMENT' not in lineread)):
            lineread = lineread.replace('D', 'e')
            block = fl.finditer(lineread)
            for m in block:
              if m.start() >= 2 and m.start() < 22:
                ephem[index][mapping[t + 1][0]] = float(m.group())
                fail = 0
              elif m.start() >= 22 and m.start() < 40 and t != 6:
                ephem[index][mapping[t + 1][1]] = float(m.group())
                fail = 0
              elif m.start() >= 40 and m.start() < 60 and t != 6:
                ephem[index][mapping[t + 1][2]] = float(m.group())
                fail = 0
              elif m.start() >= 60 and m.start() < 79 and t != 6:
                ephem[index][mapping[t + 1][3]] = float(m.group())
                fail = 0
          if(fail == 1):
            lineread = f.readline()
      index = index + 1
    else:
      for t in xrange(nlines):
        lineread = f.readline()
        while(lineread.strip() == '' or 'COMMENT' in lineread):
          lineread = f.readline()

    linedata = f.readline()

  ephem = ephem[0:index]
  f.close()
  return ephem


def parserinexsp3(filename, constellation):
  from scipy.interpolate import lagrange
  f = open(filename, 'r')
  linedata = f.readline()
  MAX_SIZE = int(linedata[32:39])
  linedata = f.readline()
  linedata = f.readline()
  nsats = int(linedata[4:6])
  ephem = {}
  for i in range(19):
    linedata = f.readline()

  ind = -1
  ephem['filename'] = filename
  ephem['Time'] = -1 * np.ones(MAX_SIZE)
  ephem['sysid'] = [['' for j in range(nsats)] for i in range(MAX_SIZE)]
  ephem['sats'] = -1 * np.ones((MAX_SIZE, nsats))
  ephem['Positions'] = -1 * np.ones((MAX_SIZE, nsats, 3))

  linedata = f.readline()
  while('EOF' not in linedata):
    if(linedata[0] == '*'):
      ind = ind + 1

      if(ind == MAX_SIZE):
        break
      year = int(linedata[3:7])
      month = int(linedata[8:10])
      day = int(linedata[11:13])
      hour = float(linedata[14:16])
      minute = float(linedata[17:19])
      second = float(linedata[20:31])
      ephem['Time'][ind], leaps = gpt.date2gps(
          [year, month, day, hour, minute, second], 0)
      for sat in range(nsats):
        linedata = f.readline()
        if('EOF' in linedata):
          return ephem
        while(linedata[0] != 'P'):
          linedata = f.readline()
          if('EOF' in linedata):
            return ephem

        ephem['sysid'][ind][sat] = linedata[1]
        ephem['sats'][ind][sat] = int(linedata[2:4])
        ephem['Positions'][ind][sat][0] = float(linedata[5:18]) * 1E3
        ephem['Positions'][ind][sat][1] = float(linedata[18:32]) * 1E3
        ephem['Positions'][ind][sat][2] = float(linedata[32:46]) * 1E3
    linedata = f.readline()
  return ephem

