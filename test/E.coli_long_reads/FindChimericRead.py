# 想要根据chimeric read的特征，将每条read的overlap分为比对到前半段和后半段的，
# 判别这两部分是否有交叉。
# 比对到前半段的，即overlap在该read的起点处hangingout很短，计在read_F中，
# 后半段，即overlap在该read的终点处hangingout很短，计到read_R中。
# 记录overlap远端到这一端的距离

def DCountOverlap(file):
    readLength = {}
    Hangingout = 50
    # suspious_reads = []
    with open(file, 'r') as blasrout:
        read_olp_count = dict()
        for line in blasrout:
            record = line.strip().split()
            (A_Orientation, B_Orientation, score)=list(map(int,record[2:5]))
            A_name,B_name = map(str,[record[0].split('.')[1],record[1].split('.')[1]])
            if A_name==B_name:
                continue
            (B_start, B_end, B_length, A_start, A_end, A_length) = list(map(int,record[6:12]))
            readLength[A_name] = A_length
            readLength[B_name] = B_length
            if B_Orientation==0 and B_length-B_end<Hangingout and A_start<Hangingout:
                if A_name+'_F' in read_olp_count:
                    read_olp_count[A_name+'_F'] += [A_end]
                else:
                    read_olp_count[A_name+'_F'] = [A_end]
                if B_name+'_R' in read_olp_count:
                    read_olp_count[B_name+'_R'] += [B_length-B_start]
                else:
                    read_olp_count[B_name+'_R'] = [B_length-B_start]
            elif B_Orientation==0 and A_length-A_end<Hangingout and B_start<Hangingout:
                if B_name+'_F' in read_olp_count:
                    read_olp_count[B_name+'_F'] += [B_end]
                else:
                    read_olp_count[B_name+'_F'] = [B_end]
                if A_name+'_R' in read_olp_count:
                    read_olp_count[A_name+'_R'] += [A_length-A_start]
                else:
                    read_olp_count[A_name+'_R'] = [A_length-A_start]
            elif B_Orientation==1 and B_length-B_end<Hangingout and A_start<Hangingout:
                if A_name+'_F' in read_olp_count:
                    read_olp_count[A_name+'_F'] += [A_end]
                else:
                    read_olp_count[A_name+'_F'] = [A_end]
                if B_name+'_F' in read_olp_count:
                    read_olp_count[B_name+'_F'] += [B_length-B_start]
                else:
                    read_olp_count[B_name+'_F'] = [B_length-B_start]
            elif B_Orientation==1 and A_length-A_end<Hangingout and B_start<Hangingout:
                if B_name+'_R' in read_olp_count:
                    read_olp_count[B_name+'_R'] += [B_end]
                else:
                    read_olp_count[B_name+'_R'] = [B_end]
                if A_name+'_R' in read_olp_count:
                    read_olp_count[A_name+'_R'] += [A_length-A_start]
                else:
                    read_olp_count[A_name+'_R'] = [A_length-A_start]
            # else:
                # print("Hangingout!",record)
                # suspious_reads+=[A_name,B_name]
    return read_olp_count,readLength

def FindChimeria(file_name,N):
    Read_olp_count,readLength = DCountOverlap(file_name)
    candidate_chimera = []
    for i in range(1,N+1):
        x = str(i)
        if x+'_F' in Read_olp_count and x+'_R' in Read_olp_count:
            LF = max(Read_olp_count[x+'_F'])+max(Read_olp_count[x+'_R'])
            # print(x,LF,readLength[x], LF>readLength[x])
            if LF<readLength[x]:
                candidate_chimera.append(x)
        else:
            candidate_chimera.append(x)
    return candidate_chimera