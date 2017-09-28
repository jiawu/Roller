import pdb
inf_methods = ['Dionesus']
organisms = [("yeast", "Yeast"), ("ecoli","Ecoli")]
                             
fobj = open("job_params_min_max_tdwindow.txt", "w")

#params = ["1_2_3_4_5","6_7_8_9_10","11_12_13_14_15","16_17_18_19_20_21"] 
params = ["0_0","0_1","0_2","0_3","0_4", "1_1","1_2","1_3", "1_4", "2_2", "2_3", "2_4","3_3", "3_4","4_4"]
window_size = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]

p = ["_".join([params[x],str(window_size[y])]) for x in range(len(params)) for y in range(len(window_size))]
# remove invalid windows
final_params = []
for param in p:
    min_l, max_l, td_window = param.split('_')
    lag_space = 21-int(td_window)
    # 21-20 = 1
    if lag_space >= int(max_l):
        final_params.append(param)

#final_p = [ "_".join([x,y,z]) for x,y,z in zip(final_params[0::3], final_params[1::3], final_params[2::3])]
nth = 75
final_p = ["_".join(final_params[x:x+nth]) for x in range(0, len(final_params), nth)]

for inf_method in inf_methods:
    for organism in organisms:
        for i in range(1,21):
            for param in final_p:
                #long_str = "/projects/p20519/roller_output/large_networks/{}/{}-insilico_size100_{} /projects/p20519/roller_output/large_networks/{}/{}-insilico_size100_{}_ /home/jjw036/Roller/data/gnw_insilico/network_data/{}/{}-{}_timeseries.tsv min_lag^max_lag pair_{}\n".format(inf_method, organism[0], i, inf_method, organism[0], i, organism[1], organism[1], i, param)
                long_str = "/projects/p20519/roller_output/gnw/{}/{}-insilico_size10_{} /projects/p20519/roller_output/gnw/{}/{}-insilico_size10_{}_ /home/jjw036/Roller/data/gnw_insilico/network_data/{}/{}-{}_timeseries.tsv min_lag^max_lag^td_window triplet_{}\n".format(inf_method, organism[0], i, inf_method, organism[0], i, organism[1], organism[1], i, param)
                fobj.write(long_str)

fobj.close()
