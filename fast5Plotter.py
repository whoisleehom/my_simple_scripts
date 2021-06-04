import os
import math
import numpy as np
import h5py as h5
import matplotlib as mpl
import matplotlib.pyplot as plt

__author__ = 'Yuan GAO'
__organization__ = 'AGIS'
__doc__ = ('A script for drawing current intensity from ONT fast5 file.')

def drawBatch(work_dir:str,batch_size:int,window_size:int,name_head:str):
    
    raw_data_list = []
    file_list = os.listdir(work_dir)
    
    for file in file_list:#read fast5 file and exract certain info
        if file[-5:] == 'fast5':
            raw_data = list(h5.File(file,'r')['/Raw/Reads'].values())[0]
            read_id  = raw_data.attrs['read_id'].decode('utf_8')
            read_len = raw_data.attrs['duration']
            raw_signal = np.asarray(raw_data[('Signal')])
            raw_data_list.append((read_id,read_len,raw_signal))
    
    raw_data_list=sorted(raw_data_list,key=lambda a:a[1])
    print
    batch_list = []
    for i in range(0,len(raw_data_list),batch_size):
        batch_list.append(raw_data_list[i:i+batch_size])
    
    fig_num = 0
    
    for j in range(len(batch_list)):
        fig = plt.figure(num=j,figsize=(10,15),dpi=300)
        fig_num += 1
        max_len = 0
        for read in batch_list[j]:
            max_len = read[1] if max_len < read[1] else max_len
        x_np= np.arange(max_len)
        #plt.title('raw data length: {}'.format(max_len))

        for i in range(len(batch_list[j])):
            pic = fig.add_subplot(batch_size,1,i+1)
            pic.set_title('read ID: {0}; raw data length: {1} + {2}'.format(batch_list[j][i][0],batch_list[j][i][1],max_len-int(batch_list[j][i][1])),fontsize=8)
            pic.set_ylabel('current intensity',fontsize=8)
            raw_np = batch_list[j][i][2]
            add_value = raw_np[-10:].mean()
            new_np = np.pad(raw_np,(0,max_len-int(batch_list[j][i][1])),'constant',constant_values=(0,add_value))
            #left_size = int(batch[i][1])%window_size
            #if  left_size !=0:
                #fill_num = window_size-left_size
                #fill_value = batch[i][2][-left_size:].mean()
                #raw_np = np.pad(raw_np,(0,fill_num),'constant',constant_values=(0,fill_value))
            #mean_np = raw_np.reshape(-1,window_size).mean(axis = -1)
            pic.plot(x_np,new_np,color='red',alpha=0.75)
            pic.tick_params(
                axis='x',          # changes apply to the x-axis
                which='both',      # both major and minor ticks are affected
                bottom=False,      # ticks along the bottom edge are off
                top=False,         # ticks along the top edge are off
                labelbottom=False) # labels along the bottom edge are off
        
        plt.savefig('{0}_batch{1}.jpg'.format(name_head,fig_num),dpi=400,bbox_inches='tight')

drawBatch('./',10,20,'RNA')