import os
import math
import numpy as np
import h5py as h5
import argparse
import matplotlib
import matplotlib.pyplot as plt
import multiprocessing as mp
matplotlib.use('Agg')

__author__ = 'Yuan GAO'
__organization__ = 'AGIS'
__doc__ = ('A script for drawing current intensity from ONT fast5 file.')


def getFileBatch(work_dir:str,batch_size:int,name_head:str):
    files = os.scandir(work_dir)
    file_count = 0
    file_info = []
    for file in files:#read fast5 file and exract certain info
        if file.name[-5:] == 'fast5':
            if file_count < 1000:
                read_len = list(h5.File(file.path,'r')['/Raw/Reads'].values())[0].attrs['duration']
                file_info.append((file.path,read_len))
                file_count += 1
            else:
                break

    file_info=sorted(file_info,key=lambda a:a[1])
    batch_list = []
    batch_count = 1
    for i in range(0,len(file_info),batch_size):
        batch_list.append((file_info[i:i+batch_size],batch_count,name_head))
        batch_count += 1
    return batch_list


def drawBatch(batch:list):
    batch_size = len(batch[0])
    batch_count = batch[1]
    name_head = batch[2]

    raw_data_list = []
    for file in batch[0]:#read fast5 file and exract certain info
        raw_data = list(h5.File(file[0],'r')['/Raw/Reads'].values())[0]
        read_id  = raw_data.attrs['read_id'].decode('utf_8')
        read_len = raw_data.attrs['duration']
        raw_signal = np.asarray(raw_data[('Signal')])
        raw_data_list.append((read_id,read_len,raw_signal))
    fig = plt.figure(num=batch_count,figsize=(10,15),dpi=300)
    max_len = 0
    for read in raw_data_list:
        max_len = read[1] if max_len < read[1] else max_len
    x_np= np.arange(max_len)
    #plt.title('raw data length: {}'.format(max_len))

    for j in range(len(raw_data_list)):
        pic = fig.add_subplot(batch_size,1,j+1)
        pic.set_title('read ID: {0}; raw data length: {1} + {2}'.format(raw_data_list[j][0],raw_data_list[j][1],max_len-raw_data_list[j][1]),fontsize=math.floor(7*10/batch_size +1))
        #pic.set_ylabel('current intensity',fontsize=8)
        raw_np = raw_data_list[j][2]
        add_value = raw_np[-10:].mean()
        new_np = np.pad(raw_np,(0,max_len-int(raw_data_list[j][1])),'constant',constant_values=(0,add_value))
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
        
        plt.savefig('{0}_batch{1}.jpg'.format(name_head,batch_count),dpi=400,bbox_inches='tight')


def getArgs():
    parser = argparse.ArgumentParser(description=__doc__,prog='uniqKmerGenerator.py')
    parser.add_argument('-d',action='store',dest='fast5_dir',required=True,type=str,
        help='the path to the folder of your fast5 files')
    parser.add_argument('-b',action='store',dest='batch_size',required=True,type=int,
        help='the number of fast5 files you want to plot in each picture')
    parser.add_argument('-o',action='store',dest='file_header',required=True,type=str,
        help='head of your jpg files')
    parser.add_argument('-t',action='store',dest='threads_num',type=int,default=5,
        help='threads you want to use. default = 5')
    return parser.parse_args()

if __name__ == '__main__':
    args = getArgs()
    batch_list = getFileBatch(args.fast5_dir,args.batch_size,args.file_header)
    real_threads = args.threads_num
    with mp.Pool(real_threads) as p:
        p.map(drawBatch,batch_list)