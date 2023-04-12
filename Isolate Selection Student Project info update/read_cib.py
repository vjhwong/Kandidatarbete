# Author(s): Yasmine Sundelin Tjärnström
# Date: Feb/March 2023

# Methods for reading, interpreting and preparing CIB for isolate selection

import re

def cut_ranges(data,fast, a,ranges,abx_abbr,parameters):
    
    #Cut BMD data to match ASTar reportable ranges
    #To avoid false on-scale etc

    abx_abbr = {v: k for k, v in abx_abbr.items()}
    try:
        kit=parameters['Kit Software Version']
    except KeyError:
        print('invalid kit software version')

    try:
        abx=abx_abbr[a]

        if fast=='Fastidious':
            try:
                abx_temp=abx+'_fast'
                range=ranges[kit][abx_temp].split(' - ')
            except:
                range=ranges[kit][abx].split(' - ')
        else:
            range=ranges[kit][abx].split(' - ')
    except:       
        return data
    
    if data[0]==0:
        return data

    if (float(data[2]) < float(range[0])) | (float(data[2]) == float(range[0])):       #if outside lower range, change to lower range
        if float(range[0])>1:
            new=range[0].split('.')[0]
        data[2] = new
        if data[1]=='=':
            data[1]='<='
            data[3]='off-scale'


    if float(data[2]) > float(range[1]):      #if outside higher range, change to higher range
        if float(range[1])>1:
            new=range[1].split('.')[0]

        data[2] = new
        if data[1]=='=':
            data[1]='>'
            data[3]='off-scale'

    return data

def D_test(iso,data):
    
    #Find true positive D-test

    CLI=iso.loc['Clindamycin'].split(' ')[0]
    ERY=iso.loc['Erythromycin'].split(' ')[0]
    D=iso.loc['D-test'].split(' ')[0]

    if D=='R':
        if (CLI=='S') & (ERY=='R'):
            data[3]='POS'
            return data
        else:
            data[3]='NEG'
            return data

    
    elif D=='S':
        data[3]='NEG'
        return data

    else:
        data[3]='-'   
        return data

def extract_data(d,j,a,ranges,abx_abbr,fast,parameters):

    #put data from CIB to variables
    tempdata=d.iloc[j][a]

    # extract and separate SIR, sign and value
    if any(char.isdigit() for char in tempdata):
    
        temp=tempdata
        if 'Missing BP' in temp:
            [MBP,temp] = temp.split('Missing BP')
            temp='Missing_BP'+temp
        [SIR, rest] = temp.split(' ')[:-1]
        temp=re.split('(\d+)', rest)[:-1]
        SIGN=temp[0]
        digits=[i for i in range(len(temp)) if temp[i].isdigit()]
        VALUE = ''.join(map(str, temp[digits[0]:]))

        #get on-scale/off-scale information
        if a=='Gentamicin':
            SCALE='on-scale'
        elif ('<' in SIGN) | ('>' in SIGN):
            SCALE='off-scale'
        else:
            SCALE='on-scale'
        
        if a=='D-test':
            [SIR,SIGN,VALUE,SCALE]=D_test(d.iloc[j],[SIR,SIGN,VALUE,SCALE])
        else:
            [SIR, SIGN,VALUE,SCALE]=cut_ranges([SIR,SIGN,VALUE,SCALE],fast,a,ranges,abx_abbr,parameters)


    else:
        [SIR,SIGN,VALUE,SCALE] = [0,0,0,0]
        

    return [SIR, SIGN,VALUE,SCALE]

def comp_data(extracted_data):

    #only relevant if both US and EU dataset
    #Compare and decide which values are valid (i.e. in panel), if both: keep both

    #both the same
    #this also includes if both has no values
    if extracted_data['US']==extracted_data['EU']:
        final_data = {'US+EU': extracted_data['US']}
        
    #only EU valid
    elif extracted_data['US'][0] == 0:
        final_data = {'EU': extracted_data['EU']}

    #only US valid
    elif extracted_data['EU'][0] == 0:
        final_data ={'US': extracted_data['US']}
    
    #both valid but diff 
    else:
        final_data = {'US': extracted_data['US'], 'EU': extracted_data['EU']}   

    return final_data

def get_data(data, j, a,ranges,abx_abbr,fast,parameters):

    #get and compare extracted data

    final_data={}    

    for t,d in data.items():

        try:
            temp=extract_data(d,j,a,ranges,abx_abbr,fast,parameters)
        
        except: #if current a not in dataset
            temp=[0,0,0,0]
        
    
        final_data[t]=temp

    if len(final_data)>1:
        final_data=comp_data(final_data)
    
    return final_data

def rank_system(res: dict, point_system: dict): 

   #Find info about SIR and on/offscale and give point (predefined)
   #Add point to isolate, which will be included in final dataset and used to sort it
    
    points=0
    res_temp=dict(list(res.items())[3:])
    for key1, value1 in res_temp.items():
        for key2, value2 in value1.items():
            for v in [value2[0],value2[3]]:

                try:
                    points+=point_system[v]
                except KeyError:
                    points+=0

    res['Q-rank']=points

    return res