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
import math
import supportfunctions as gpt

global mapping
mapping = [' ' for x in xrange(15)]
global MAX_SIZE
MAX_SIZE = 86400
global constellations
constellations = {}
global alldata
alldata = {}
global constellationmap
constellationmap = {
    'G': 'GPS',
    'R': 'GLONASS',
    'E': 'Galileo',
    'C': 'BeiDou',
    'J': 'QZSS',
    'S': 'SBAS'
}
fl = re.compile('[-+]?(\d+(\.\d*)?|\.\d+)([e][-+]?\d+)?')
integ = re.compile('[-+]?\d+')


def initializedata():
  global alldata
  global constellations
  for key in constellations:
    if constellations[key]['Enabled'] == 1:
      alldata[key] = {
          'sats': [[] for x in xrange(MAX_SIZE)],
          'numcsats': np.zeros((MAX_SIZE, 1))
      }
      alldata[key]['PRNseq'] = constellations[key]['PRN']


def initializeconstellation(key):
  global alldata
  alldata[key] = {
      'sats': [[] for x in xrange(MAX_SIZE)],
      'numcsats': np.zeros((MAX_SIZE, 1))
  }
  alldata[key]['PRNseq'] = constellations[key]['PRN']


def mapfunction(word):
  return float(word[0])


def checkforkeywords(stri):
  answer = 1
  stri = stri.strip()
  stri = stri.upper()
  keywords = [
      'COMMENT', 'MARKERNAME', 'MARKERNUMBER', 'APPROXPOSITIONXYZ',
      'ANTENNA:DELTAH/E/N', 'ENDOFHEADER']
  if not filter(lambda x: x in stri, keywords):
    answer = 0
  return answer


def floatm(word):
  word = word.strip()
  word = word.replace('D', 'e')
  try:
    out = float(word)
  except:
    out = -1
  return out

# Modify this function later


def organizedatatypes():
  global constellations
  global mapping
  for typ in mapping:
    if(typ == 'L1' or typ == 'LA'):
      for key in constellations:
        if(constellations[key]['Enabled']):
          alldata[key][
              'L1'] = -1 * np.ones((MAX_SIZE, constellations[key]['Nsats']))
          alldata[key][
              'L1quality'] = -1 * np.ones((MAX_SIZE, constellations[key]['Nsats']))
    elif(typ == 'L2' or typ == 'LC'):
      for key in constellations:
        if(constellations[key]['Enabled']):
          alldata[key][
              'L2'] = -1 * np.ones((MAX_SIZE, constellations[key]['Nsats']))
    elif(typ == 'C1' or typ == 'CA'):
      for key in constellations:
        if(constellations[key]['Enabled']):
          alldata[key][
              'C1'] = -1 * np.ones((MAX_SIZE, constellations[key]['Nsats']))
    elif(typ == 'P2' or typ == 'CC'):
      for key in constellations:
        if(constellations[key]['Enabled']):
          alldata[key][
              'P2'] = -1 * np.ones((MAX_SIZE, constellations[key]['Nsats']))
    elif(typ == 'S1' or typ == 'SA'):
      for key in constellations:
        if(constellations[key]['Enabled']):
          alldata[key][
              'S1'] = -1 * np.ones((MAX_SIZE, constellations[key]['Nsats']))
    elif(typ == 'S2' or typ == 'SC'):
      for key in constellations:
        if(constellations[key]['Enabled']):
          alldata[key][
              'S2'] = -1 * np.ones((MAX_SIZE, constellations[key]['Nsats']))
    elif(typ == 'D1' or typ == 'DA'):
      for key in constellations:
        if(constellations[key]['Enabled']):
          alldata[key][
              'D1'] = -1 * np.ones((MAX_SIZE, constellations[key]['Nsats']))
    elif(typ == 'D2' or typ == 'DC'):
      for key in constellations:
        if(constellations[key]['Enabled']):
          alldata[key][
              'D2'] = -1 * np.ones((MAX_SIZE, constellations[key]['Nsats']))
    elif(typ != ' '):
      for key in constellations:
        if(constellations[key]['Enabled']):
          alldata[key][
              typ] = -1 * np.ones((MAX_SIZE, constellations[key]['Nsats']))


def parseepoch(f):
  eof = 0
  time = -1
  sats = -1
  numsats = -1
  sepcounts = {
      'GPS': 0,
      'GLONASS': 0,
      'Galileo': 0,
      'BeiDou': 0,
      'QZSS': 0,
      'SBAS': 0
  }
  while(eof == 0):
    origdata = f.readline()
    answer = checkforkeywords(origdata)
    # Read valid line unless eof
    while((answer == 1 or origdata.strip() == '') and origdata != ''):
      origdata = f.readline()

      answer = checkforkeywords(origdata)
    if(origdata == ''):
      eof = 1
      return [-1, -1, -1, -1, 1]

    temparray = ['G', 'R', 'E', 'C', 'J', 'S']
    # Valid line read
    # REVIEW THIS CONDITION
    # RINEX version 2.xx
    if(origdata[0] != '>' and any(x in origdata for x in temparray)):
      timparams = fl.findall(origdata[0:26])
      if(len(timparams) < 4):
        continue
      elif(len(timparams) < 6):
        time = -1
      else:
        [year, month, day, hour, minute, second] = map(
            mapfunction, timparams[0:6])
        year = gpt.four_digit_year(year)
        time = gpt.date2gps([year, month, day, hour, minute, second], 0)
      try:
        numsats = int(integ.search(origdata[29:32]).group())
      except:
        continue
      nlines = int(math.ceil(numsats / 12.0))
      substring = ''
      substring = substring + origdata[32:68]
      for n in range(nlines):
        if(n > 0):
          origdata = f.readline()
          while(origdata.strip() == '' and origdata != ''):
            origdata = f.readline()
          if(origdata == ''):
            eof = 1
            break
          substring = substring + origdata[32:68]
      sats = -1 * np.ones(numsats)
      sattypes = [' ' for t in range(numsats)]
      for s in range(numsats):
        sattypes[s] = substring[3 * s]
        sats[s] = floatm(substring[3 * s + 1:3 * s + 3])
        sepcounts[constellationmap[sattypes[s]]] += 1
      # return time,numsats,sepcounts,sattypes,sats,eof;
      return time, numsats, sattypes, sats, eof
    elif(origdata[0] == '>'):
      timparams = fl.findall(origdata[1:29])
      if(len(timparams) < 4):
        continue
      elif(len(timparams) < 6):
        time = -1
      else:
        [year, month, day, hour, minute, second] = map(
            mapfunction, timparams[0:6])
        year = gpt.four_digit_year(year)
        time = gpt.date2gps([year, month, day, hour, minute, second], 0)
      try:
        numsats = int(integ.search(origdata[32:35]).group())
      except:
        continue
      return time, numsats, -1, -1, eof


def parseheader(f):
  """ Header operations """
  global alldata
  global mapping
  global constellations

  fvisit = 0
  origdata = f.readline()
  linedata = ''.join(origdata.split())
  linedata = linedata.upper()
  found_types = 1
  while(linedata.find('ENDOFHEADER') == -1 and origdata != ''):
    if(linedata.find('LEAPSECONDS') != -1):
      listdata = origdata.split()
      alldata['leapseconds'] = int(listdata[0])
      # RINEX v2.xx
    elif(linedata.find('#/TYPESOFOBSERV') != -1):
      alldata['Version'] = 2
      if(fvisit == 0):
        initializedata()
      fvisit = 1
      nObs = int(integ.search(origdata[0:7]).group())
      nlines = int(math.ceil(nObs / 9.0))
      alldata['nObs'] = nObs
      mapping = mapping[0:nObs]
      for i in range(nlines):
        if(i > 0):
          origdata = f.readline()
        nwords = min(nObs, 9)
        for k in range(nwords):
          typ = ''.join(origdata[(k + 1) * 6:(k + 1) * 6 + 6].split())
          # outtyp=organizedatatypes2(typ,2);
          # outtyp=typ[0:2];
          mapping[i * 9 + k] = typ
        nObs = nObs - 9
      found_types = 1
      # RINEX v3.xx
    elif(linedata.find('SYS/#/OBSTYPES') != -1):
      alldata['Version'] = 3
      alldata['sysid'].append(origdata[0])
      if(constellations[constellationmap[origdata[0]]]['Enabled']):
        initializeconstellation(constellationmap[origdata[0]])
        nObs = int(integ.search(origdata[1:7]).group())
        nlines = int(math.ceil(nObs / 13.0))
        alldata[constellationmap[origdata[0]]]['nObs'] = nObs
        alldata[constellationmap[origdata[0]]]['obs'] = [' '
                                                         for x in range(nObs)]
        for i in range(nlines):
          if(i > 0):
            origdata = f.readline()
          nwords = min(nObs, 13)
          for k in range(nwords):
            typ = ''.join(origdata[6 + k * 4:6 + k * 4 + 4].split())
            outtyp = typ[0:3]
            alldata[constellationmap[origdata[0]]][
                outtyp] = -1 * np.ones((MAX_SIZE, constellations[constellationmap[origdata[0]]]['Nsats']))
            alldata[constellationmap[origdata[0]]]['obs'][i * 9 + k] = outtyp
            mapping[i * 9 + k] = outtyp
          nObs = nObs - 13
        found_types = 1
    elif(linedata.find('APPROXPOSITIONXYZ') != -1):
      X = float(fl.search(origdata[0:14]).group())
      Y = float(fl.search(origdata[14:28]).group())
      Z = float(fl.search(origdata[28:42]).group())
      alldata['Position'] = np.array([X, Y, Z])
    elif(linedata.find('INTERVAL') != -1):
      alldata['Interval'] = float(fl.search(origdata[0:10]).group())
    origdata = f.readline()
    linedata = ''.join(origdata.split())
  return found_types


def parserinexobs(filename, i):
  global constellations
  global mapping
  global constellationmap
  global alldata
  global MAX_SIZE

  constellations = i
  mapping = [' ' for x in xrange(12)]
  index = filename.rfind('\\')
  surname = filename[index + 1:]
  alldata = {
      'filename': surname,
      'sysid': [],
      'leapseconds': 16,
      'GPStime': np.zeros((MAX_SIZE, 1)),
      'numsats': np.zeros((MAX_SIZE, 1))
  }
  f = open(filename, 'r')
  success = parseheader(f)
  if(alldata['Version'] == 2):
    organizedatatypes()
  if(success == 1):
    eof = 0
    ind = 0
    while(eof == 0):
      [alldata['GPStime'][ind], alldata['numsats']
       [ind], sattypes, sats, eof] = parseepoch(f)
      if(eof == 1):
        break
      if(alldata['Version'] == 2):
        nsats = alldata['numsats'][ind]
        nObs = alldata['nObs']
        # assignsatsize(sepcounts,ind);
        nlines = int(math.ceil(alldata['nObs'] / 5.0))
        for n1 in range(nsats):
          if(eof == 1):
            break
          if(constellations[constellationmap[sattypes[n1]]]['Enabled']):
            alldata[constellationmap[sattypes[n1]]]['numcsats'][ind] += 1
            alldata[constellationmap[sattypes[n1]]][
                'sats'][ind].append(int(sats[n1]))
            for n2 in range(nlines):
              origdata = f.readline()
              answer = checkforkeywords(origdata)
              # while((answer==1 or origdata.strip()=='') and origdata!=''):
              while(answer == 1 and origdata != ''):
                origdata = f.readline()
                answer = checkforkeywords(origdata)
              if(origdata == ''):
                eof = 1
                print 'Itshere'
                break
              iterlength = len(origdata)
            # Assuming a maximum of 5 observations per line
              if(sattypes[n1] == 'J'):
                sub = sats[n1] - 198
              elif(sattypes[n1] == 'S'):
                sub = alldata['SBAS']['PRNseq'].index((int(sats[n1]) + 100))
              else:
                sub = sats[n1] - 1

              if(iterlength > 3 and n2 * 5 + 0 < nObs):
                alldata[constellationmap[sattypes[n1]]][
                    mapping[n2 * 5 + 0]][ind][sub] = floatm(origdata[1:14])
                if(mapping[n2 * 5 + 0] == 'L1'):
                  alldata[constellationmap[sattypes[n1]]][
                      'L1quality'][ind][sub] = floatm(origdata[15:16])
              if(iterlength > 18 and n2 * 5 + 1 < nObs):
                alldata[constellationmap[sattypes[n1]]][
                    mapping[n2 * 5 + 1]][ind][sub] = floatm(origdata[17:30])
                if(mapping[n2 * 5 + 1] == 'L1'):
                  alldata[constellationmap[sattypes[n1]]][
                      'L1quality'][ind][sub] = floatm(origdata[31:32])
              if(iterlength > 34 and n2 * 5 + 2 < nObs):
                alldata[constellationmap[sattypes[n1]]][
                    mapping[n2 * 5 + 2]][ind][sub] = floatm(origdata[33:46])
                if(mapping[n2 * 5 + 2] == 'L1'):
                  alldata[constellationmap[sattypes[n1]]][
                      'L1quality'][ind][sub] = floatm(origdata[47:48])
              if(iterlength > 50 and n2 * 5 + 3 < nObs):
                alldata[constellationmap[sattypes[n1]]][
                    mapping[n2 * 5 + 3]][ind][sub] = floatm(origdata[49:62])
                if(mapping[n2 * 5 + 3] == 'L1'):
                  alldata[constellationmap[sattypes[n1]]][
                      'L1quality'][ind][sub] = floatm(origdata[63:64])
              if(iterlength > 65 and n2 * 5 + 4 < nObs):
                alldata[constellationmap[sattypes[n1]]][
                    mapping[n2 * 5 + 4]][ind][sub] = floatm(origdata[65:78])
                if(mapping[n2 * 5 + 4] == 'L1'):
                  alldata[constellationmap[sattypes[n1]]][
                      'L1quality'][ind][sub] = floatm(origdata[79:80])
          else:
            # Skip lines
            for n2 in range(nlines):
              origdata = f.readline()
              answer = checkforkeywords(origdata)
              # while((answer==1 or origdata.strip()=='') and origdata!=''):
              while(answer == 1 and origdata != ''):
                origdata = f.readline()
                answer = checkforkeywords(origdata)
              if(origdata == ''):
                eof = 1
                break
      else:
        nsats = alldata['numsats'][ind]
        for n1 in range(nsats):
          origdata = f.readline()
          answer = checkforkeywords(origdata)
          # while((answer==1 or origdata.strip()=='') and origdata!=''):
          while(answer == 1 and origdata != ''):
            origdata = f.readline()
            answer = checkforkeywords(origdata)
          if(origdata == ''):
            eof = 1
            break
          cid = origdata[0]
          sysid = ''.join(alldata['sysid'])
          if (
              constellations[constellationmap[cid]]['Enabled'] and (
                  sysid.find(
                      cid) != -
                  1)):
            satid = floatm(origdata[1:3])
            if(cid == 'J'):
              sub = satid - 198
            elif(cid == 'S'):
              sub = alldata['SBAS']['PRNseq'].index((satid + 100))
            else:
              sub = satid - 1
            alldata[constellationmap[cid]]['numcsats'][ind] += 1
            alldata[constellationmap[cid]]['sats'][ind].append(int(satid))
            substring = origdata[3:]
            nObs = alldata[constellationmap[cid]]['nObs']
            startindex = 1
            l = len(substring)
            data = substring[startindex:startindex + 13]
            count = 0
            while(startindex <= l - 13):
              outtyp = alldata[constellationmap[origdata[0]]]['obs'][count]
              alldata[constellationmap[origdata[0]]][
                  outtyp][ind][sub] = floatm(data)
              startindex = startindex + 16
              data = substring[startindex:startindex + 13]
              count = count + 1
      ind = ind + 1
      # Add indices if necessary
      if(ind >= MAX_SIZE):
        print('Reached limit')
        eof = 1
  else:
    print 'Parsing failed'
    return

  alldata['GPStime'] = alldata['GPStime'][0:ind]
  alldata['numsats'] = alldata['numsats'][0:ind]
  if(alldata['Version'] == 2):
    for key in constellations:
      if constellations[key]['Enabled'] == 1:
        for k in range(alldata['nObs']):
          alldata[key][mapping[k]] = alldata[key][mapping[k]][0:ind]
        alldata[key]['sats'] = alldata[key]['sats'][0:ind]
        alldata[key]['numcsats'] = alldata[key]['numcsats'][0:ind]
  elif(alldata['Version'] == 3):
    for key in constellations:
      if constellations[key]['Enabled'] == 1:
        for k in range(alldata[key]['nObs']):
          val = alldata[key]['obs'][k]
          alldata[key][val] = alldata[key][val][0:ind]
        alldata[key]['sats'] = alldata[key]['sats'][0:ind]
        alldata[key]['numcsats'] = alldata[key]['numcsats'][0:ind]
  return alldata

