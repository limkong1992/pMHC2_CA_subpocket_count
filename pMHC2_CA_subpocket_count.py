import numpy as np
import pandas as pd
from biopandas.pdb import PandasPdb

# The dist_cutoff is set to 2.5A
dist_cutoff = 2.5 
# The number of subpockets of the 9 subpockets that count 
match_count_cutoff = 0
# The reference peptide
ref_pep_file = 'ref_core.pdb'
# The peptides aligned to the reference peptide
peptide_file_list = ['p1.pdb','p2.pdb','p3.pdb']
# The output csv file
out_csv_file = 'result_25A.csv'

# calculate the eucli distance
def eucliDist(A,B):
    return np.sqrt(sum(np.power((A-B),2)))

# transform the amino acid name from three letter to one letter form
def three_to_one_seq(list3):
    one_letter = {'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', 'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',    \
'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A', 'GLY':'G', 'PRO':'P', 'CYS':'C'}
    three_letter = {'V':'VAL', 'I':'ILE', 'L':'LEU', 'E':'GLU', 'Q':'GLN', 'D':'ASP', 'N':'ASN', 'H':'HIS', 'W':'TRP', 'F':'PHE', 'Y':'TYR',    \
'R':'ARG', 'K':'LYS', 'S':'SER', 'T':'THR', 'M':'MET', 'A':'ALA', 'G':'GLY', 'P':'PRO', 'C':'CYS'}
    list_1 = list(range(len(list3)))
    for i in range(len(list3)):
        list_1[i] = one_letter[list3[i]]

    string_list = ''.join(list_1)
    return string_list

def three_to_one_reverse_seq(list3):
    one_letter = {'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', 'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',    \
'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A', 'GLY':'G', 'PRO':'P', 'CYS':'C'}
    three_letter = {'V':'VAL', 'I':'ILE', 'L':'LEU', 'E':'GLU', 'Q':'GLN', 'D':'ASP', 'N':'ASN', 'H':'HIS', 'W':'TRP', 'F':'PHE', 'Y':'TYR',    \
'R':'ARG', 'K':'LYS', 'S':'SER', 'T':'THR', 'M':'MET', 'A':'ALA', 'G':'GLY', 'P':'PRO', 'C':'CYS'}
    list_1 = list(range(len(list3)))
    for i in range(len(list3)):
        list_1[i] = one_letter[list3[i]]
    list_1.reverse()
    string_list = ''.join(list_1)
    return string_list


def find_best_matche_position(peptide_file):
    # read the pdb files
    p1 = PandasPdb().read_pdb(peptide_file)
    data1 = p1.df['ATOM']
    
    # screen the CA coordinates
    p1_ca = data1.loc[data1['atom_name']=='CA'][['residue_name','x_coord', 'y_coord', 'z_coord','b_factor']]
    p1_ca.index = range(len(p1_ca))
    
    p1_seq = list(p1_ca['residue_name'])
    p1_seq_one = three_to_one_seq(p1_seq)
    p1_seq_one_r = three_to_one_reverse_seq(p1_seq)
    #save the b_factor
    p1_plddt = list(p1_ca['b_factor'])
    p1_plddt_r = list(reversed(p1_plddt))

    for pref_ca_id in range(len(pref_ca)):
        pref_ca_co = np.array(pref_ca.loc[pref_ca_id][1:4])
        for p1_ca_id in range(len(p1_ca)):
            p1_ca_co = np.array(p1_ca.loc[p1_ca_id][1:4])
            ca_dist_i = eucliDist(pref_ca_co,p1_ca_co)
            if ca_dist_i<=dist_cutoff:
                matched_df.loc[(p1_ca_id, pref_ca_id)]=1
            else:
                matched_df.loc[(p1_ca_id, pref_ca_id)]=0
            ca_dist_df.loc[(p1_ca_id, pref_ca_id)]= round(ca_dist_i,1)  
    #print(ca_dist_df)
    print(matched_df)
    
    # select the forward matched position
    #best_matched = []
    matched_count = match_count_cutoff
    s1matchpid = 0
    result_align_list = []
    
    for j in range(len(p1_ca)-8):
        i = 0
        match_j = j
        align_list = []
        while j < len(p1_ca) and i < len(pref_ca):
            align_list.append(matched_df.loc[(j,i)])
            i = i+1
            j = j+1
        if sum(align_list) > matched_count:
            matched_count = sum(align_list)
            s1matchpid = match_j+1
            print('Forward: peptide no %d match peptide-ref no %d best. ' %(match_j+1, 1))
            print(align_list)
            result_align_list = align_list
    fr_l = p1_seq_one[:s1matchpid-1]
    core_seq = p1_seq_one[s1matchpid-1:s1matchpid+8]
    fr_r = p1_seq_one[s1matchpid+8:]
    
    fr_l_plddt = p1_plddt[:s1matchpid-1]
    core_plddt = p1_plddt[s1matchpid-1:s1matchpid+8]
    fr_r_plddt = p1_plddt[s1matchpid+8:]
    core_plddt_mean = round(sum(core_plddt)/9,2)
    
    # select the reversed matched position
    #best_matched_r = []
    matched_count_r = match_count_cutoff
    s9matchpid = 0
    result_align_list_r = []

    for j in range(len(p1_ca)-8):
        i = 8
        match_j = j
        align_list = []
        while j < len(p1_ca) and i >-1:
            align_list.append(matched_df.loc[(j,i)])
            i = i-1
            j = j+1
        if sum(align_list) > matched_count_r:
            matched_count_r = sum(align_list)
            s9matchpid = match_j+1
            print('Reverse: peptide no %d match peptide-ref no %d best. ' %(match_j+1, 9))
            print(align_list)
            result_align_list_r = align_list
    print(p1_seq_one_r)
    print(s9matchpid)
    print(p1_seq_one_r)
    if s9matchpid == 1:
        fr_l_r = p1_seq_one_r[:-s9matchpid-8]
        core_seq_r = p1_seq_one_r[-9:]
        fr_r_r = []
        fr_l_plddt_r = p1_plddt_r[:-s9matchpid-8]
        core_plddt_r = p1_plddt_r[-9:]
        fr_r_plddt_r = []
    else:
        fr_l_r = p1_seq_one_r[:-s9matchpid-8]
        core_seq_r = p1_seq_one_r[-s9matchpid-8:-s9matchpid+1]
        fr_r_r = p1_seq_one_r[-s9matchpid+1:]

        fr_l_plddt_r = p1_plddt_r[:-s9matchpid-8]
        core_plddt_r = p1_plddt_r[-s9matchpid-8:-s9matchpid+1]
        fr_r_plddt_r = p1_plddt_r[-s9matchpid+1:]


    core_plddt_mean_r = round(sum(core_plddt_r)/9,2)

    if matched_count_r > matched_count:
        matched_count_result = matched_count_r
        matched_direction = 'Reverse'
        core_plddt_mean_result = core_plddt_mean_r 
    else:
        matched_count_result = matched_count
        matched_direction = 'Forward'
        core_plddt_mean_result = core_plddt_mean
        

    return [p1_seq_one, matched_direction, matched_count_result, core_plddt_mean_result,\
            matched_count, s1matchpid, fr_l,core_seq, fr_r, fr_l_plddt, core_plddt, fr_r_plddt, result_align_list, core_plddt_mean,\
            matched_count_r,s9matchpid,fr_l_r,core_seq_r,fr_r_r,fr_l_plddt_r,core_plddt_r,fr_r_plddt_r,result_align_list_r, core_plddt_mean_r]
    

pref = PandasPdb().read_pdb(ref_pep_file)
dataref = pref.df['ATOM']
pref_ca = dataref.loc[dataref['atom_name']=='CA'][['residue_name','x_coord', 'y_coord', 'z_coord']]
pref_ca.index = range(len(pref_ca))
# caculate the CA-CA distance
ca_dist_df =  pd.DataFrame(columns = range(len(pref_ca)))
matched_df =  pd.DataFrame(columns = range(len(pref_ca)))

mcl = []
for peptide_file in peptide_file_list:
    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"+peptide_file)
    mc = find_best_matche_position(peptide_file)
    mcl.append(mc)
#print(mcl)
result_df = pd.DataFrame(mcl, columns = ['Sequence','Match_Direction','Match_Count_Result','Core_Plddt_Mean_Result',\
        'Match_Count_Forward','Site_1_AA_Id_Forward','FR_Left_Forward','Core_Forward','FR_Right_Forward','FR_Left_Forward_Plddt','Core_Forward_Plddt','FR_Right_Forward_Plddt','AA_In_Site_Forward','Core_Forward_Plddt_Mean',\
        'Match_Count_Reverse','Site_9_AA_Id_Reverse','FR_Left_Reverse','Core_Reverse','FR_Right_Reverse','FR_Left_Reverse_Plddt','Core_Reverse_Plddt','FR_Right_Reverse_Plddt','AA_In_Site_Reverse','Core_Reverse_Plddt_Mean'])
print(result_df)
result_df.to_csv(out_csv_file,index=False, header=True, sep=',')
