"""This is a stand alone function specifically designed to get METAR observations from UWYO,
    May update later to run inline with other class objects within the package, otherwise, it works fine as is,
    so why "fix it" """

import urllib  # the lib that handles the url stuff
import re
import datetime as DT
import pandas as pd
import os
import glob as glob
import numpy as np
import sys

def main():
    """Main function where all other functions are called from"""
    try:
        input_namelist=sys.argv[1]
    except:
        print("NO NAMELIST ARGUMENT, DEFAULTING TO 'wrf_post_namelist.txt'")
        input_namelist='stn_namelist.txt'

    namelist_dictionary=parse_namelist(namelist_file=input_namelist)

    print(namelist_dictionary)

    stns=namelist_dictionary['stns']
    date_start=namelist_dictionary['time_start']
    date_end=namelist_dictionary['time_end']
    tfmt=namelist_dictionary['time_format']
    output_path=namelist_dictionary['output_path']
    unit_choice=namelist_dictionary['units']
    pfx=namelist_dictionary['pfx']

    date_start_strp=DT.datetime.strptime(date_start,tfmt)
    date_end_strp=DT.datetime.strptime(date_end,tfmt)

    total_days=(date_end_strp-date_start_strp).total_seconds()/86400.

    keys=['TIME','PMSL','ALTM','TMP','DEW','RH','DIR','SPD','VIS','Weather']

    if unit_choice == 'M':
        unit_dict={'TIME':'NA','PMSL':'hPa','ALTM':'hPa','TMP':'C','DEW':'C','RH':'%',
                   'DIR':'Degrees CC from North','SPD':'m/s','VIS':'km'}
    if unit_choice == 'A':
        unit_dict={'TIME':'NA','PMSL':'hPa','ALTM':'inHg','TMP':'F','DEW':'F','RH':'%',
                   'DIR':'Degrees CC from North','SPD':'kt','VIS':'mile'}

    print('-----------------------------------------------------')
    print('|  Collecting Data from UWYO Surface data archive   |')

    for stn in stns:
        start_collect=False
        outputfile=output_path+'/%s_%s.csv'%(pfx,stn)
        print('-----------------------------------------------------')
        print("--- Saving data to: %s"%outputfile)
        print("--- Start Date    : %s"%date_start_strp.strftime('%c'))
        print("--- End Date      : %s"%date_end_strp.strftime('%c'))
        print("--- Total Days    : %i"%total_days)
        print("--- Units in %s (M = metric, A= Imperial)"%unit_choice)
        print("--- Saving Vars   : "+ ','.join(['%s'%i for i in keys]))

        ## GET HEADER INFORMATION ##
        stn_search='<H3>'
        loc_search='<H5>'
        head_search='STN'
        start_collect=False
        unit_search='DD/HHMM'
        start_search="=="
        date_search='<I>'
        ###

        data_dict={}
        for i in keys:
            data_dict[i]=[] ## initialize lists for appending!

        ### NOW DO THE DATA THING!  ##

        for i in range(int(total_days)+1):
            current_date_strf=(date_start_strp+DT.timedelta(hours=24*i)).strftime('%Y%m%d')
            target_url='http://weather.uwyo.edu/cgi-bin/wyowx.fcgi?TYPE=sflist&DATE=%s&HOUR=24&UNITS=%s&STATION=%s'%(current_date_strf,unit_choice,stn)
            data = urllib.request.urlopen(target_url) # it's a file like object and works just like a file

            for line in data: # files are iterable
                if start_collect == True:
                    if 'PRE' in line.decode("utf-8"):
                        ## END OF LINE:
                        start_collect=False
                        break
                    line_data=line.decode("utf-8")
                    str_idx=0 ## start at ZERO index.
                    for idx, i in enumerate(length_id):
                        if header[idx] in keys:
                            ## GET THIS DATA! ##
                            if header[idx] == 'TIME':#
                                if int(line_data[str_idx:str_idx+i+1].strip().split('/')[0]) == int(day1):
                                    month,year=month1,year1
                                else:
                                    month,year=month2,year2
                                data_dict[header[idx]].append(DT.datetime.strptime('%s-%s-%s'%(year,month,line_data[str_idx:str_idx+i+1].strip()),'%Y-%b-%d/%H%M'))
                                ## SPECIAL CASE, TIME VALUES! ##
                            else:
                                ## IT's NOT TIME
                                data_dict[header[idx]].append(line_data[str_idx:str_idx+i+1].strip())

                        str_idx+=(i+1) ## move down the line! ##

                    for k in no_key:
                        data_dict[k].append('NaN')

                else:

                    if loc_search in line.decode("utf-8"):
                        lat,lon=parse_latlon(line.decode("utf-8"))

                    if date_search in line.decode("utf-8"):
                        date_str=line.decode("utf-8").split()
                        day1,month1,year1=date_str[1:4]
                        day2,month2,year2=date_str[6:9]
                        year2=year2[:4]

                    if head_search in line.decode("utf-8"):
                        header=line.decode("utf-8").split()
                        no_key=[]
                        for k in keys:
                            if k not in header:
                                no_key.append(k)

                    if start_search in line.decode("utf-8"):
                        length_id=[len(i) for i in line.decode("utf-8").split()]
                        c=1
                        while(len(length_id) > len(header)):
                            ## NEED TO ACCOUNT FOR THE FACT THAT CLOUDS CAN HAVE MULTIPLE LEVELS. ##
                            if c == 1:
                                idx=header.index('CLOUDS')
                            else:
                                idx=header.index('CLOUDS%i'%(c-1))
                            header.insert(idx+1, 'CLOUDS%i'%c)
                            c=c+1
                        start_collect=True

                    if unit_search in line.decode('utf-8'):
                        units=line.decode("utf-8").split()

                    if stn_search in line.decode("utf-8"):
                        metaname=line.decode("utf-8").split('>')[1].split('<')[0]
                        if i == 0:
                            print("--- Station Name : %s"%metaname)
                            print('-----------------------------------------------------')
                        unitname='|'.join([k+':'+unit_dict[k] for k in unit_dict])
            print("Data Saved from: %s"%target_url)
            data.close()
            
        dataframe=pd.DataFrame.from_dict(data_dict)
        dataframe['latitude']=lat
        dataframe['longitude']=lon
        ## OUTPUT TO CSV FILE ##
        dataframe.sort_values(by=['TIME'], inplace=True)
        with open(outputfile, 'w') as f:
            f.write('%s\n'%metaname)
            f.write('%s\n'%unitname)
        dataframe.to_csv(outputfile,index=False,mode='a')

        print("Finished! Data written to %s"%outputfile)
    print("FINISHED ALL FILES!")

def parse_latlon(llstr):
    latlon=re.findall(r"[-+]?\d*\.\d+|\d+", llstr)
    latlon=[float(i) for i in latlon]
    if 'W' in llstr:
        ## WESTERN HEMISPHERE: LON = LON*-1
        latlon[2]=latlon[2]*-1
    return latlon[1:3]


def parse_namelist(namelist_file='stn_namelist.txt'):
    import os
    cwd=os.getcwd()

    lookforkeys=['stns','time_start','time_end','time_format','units','output_path','pfx']
    list_keys=['stns']

    try:
        namelistfile=glob.glob(cwd+'/'+namelist_file)[0]
    except:
        print("NAMELIST FILE: %s does not exist, check your paths are try again!")
        sys.exit("EXITING...")


    ## SET DEFAULTS HERE ##
    namelist_dictionary={}

    namelist_dictionary['stns']=['KLEB']
    namelist_dictionary['time_format']='%Y%m%d'
    namelist_dictionary['pfx']=''
    namelist_dictionary['units']='M'
    namelist_dictionary['time_start']='20191001'
    namelist_dictionary['time_end']='20191005'
    namelist_dictionary['output_path']=cwd

    ## OKAY, DEFAULTS SET! ##
    with open(namelistfile) as fp:
        for line in fp:
            var_key=line.split('=')[0].strip(' \t\n\r')
            if var_key in lookforkeys and '#' not in var_key: ## if not a comment!
                values=line.split('=')[1].split('#')[0].strip(' \t\n\r').split(',')
                values=[x for x in values if x != '']
                if len(values) > 1:
                    if var_key in list_keys:
                        namelist_dictionary[var_key]=values
                    else:
                        print("WARNING, THIS IS NOT A LIST VARIABLE")
                        print("SETTING KEY VALUE TO FIRST OPTION")
                        namelist_dictionary[var_key]=values[0]
                else:
                    if var_key in list_keys:
                        namelist_dictionary[var_key]=values
                    else:
                        namelist_dictionary[var_key]=values[0]

    return namelist_dictionary

if __name__ == '__main__':
    main()
