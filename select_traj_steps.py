## for giving trajectory/ies chose every defined number of frames
import numpy as np
import os

def find_lines_in_model(my_file, out): #find first HEADER and calculate number of positions for one model
    find = 0
    while 1:
        line = my_file.readline()
        if "HEADER" in line:
            pos_in_line = len(line)
            find += 1
            if find ==2: break
        out.write(line)
    pos_for_model = my_file.tell()-pos_in_line
    return pos_for_model
    
def write_model(m_file, out_f):
    find = 0
    while 1:
        line = m_file.readline()
        if "HEADER" in line:
            find += 1
            if find == 2: break
        out_f.write(line)


N = 10 #number of skipping models    
traj = open("test_A.pdb")
out_traj = open("output_test.pdb","w")
file_size = os.fstat(traj.fileno()).st_size
lines_for_model  = find_lines_in_model(traj, out_traj)
print lines_for_model
traj.seek(0)
while 1:
    traj.seek((N-2)*lines_for_model, 1)
    print traj.tell(), file_size
    if (traj.tell() + (N-2)*lines_for_model) >= file_size: 
        out_traj.close()
        break
    else: write_model(traj, out_traj)

 

