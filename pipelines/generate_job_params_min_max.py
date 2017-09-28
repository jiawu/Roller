inf_methods = ['Dionesus']
organisms = [("yeast", "Yeast"), ("ecoli","Ecoli")]
                             
fobj = open("job_params_min_max.txt", "w")

#params = ["1_2_3_4_5","6_7_8_9_10","11_12_13_14_15","16_17_18_19_20_21"] 
params = ["0_0_0_1_0_2_0_3_0_4", "1_1_1_2_1_3_1_4", "2_2_2_3_2_4","3_3_3_4","4_4"]

for inf_method in inf_methods:
    for organism in organisms:
        for i in range(1,21):
            for param in params:
                #long_str = "/projects/p20519/roller_output/large_networks/{}/{}-insilico_size100_{} /projects/p20519/roller_output/large_networks/{}/{}-insilico_size100_{}_ /home/jjw036/Roller/data/gnw_insilico/network_data/{}/{}-{}_timeseries.tsv min_lag^max_lag pair_{}\n".format(inf_method, organism[0], i, inf_method, organism[0], i, organism[1], organism[1], i, param)
                long_str = "/projects/p20519/roller_output/gnw/{}/{}-insilico_size10_{} /projects/p20519/roller_output/gnw/{}/{}-insilico_size10_{}_ /home/jjw036/Roller/data/gnw_insilico/network_data/{}/{}-{}_timeseries.tsv min_lag^max_lag pair_{}\n".format(inf_method, organism[0], i, inf_method, organism[0], i, organism[1], organism[1], i, param)
                fobj.write(long_str)

fobj.close()
