inf_methods = ['Dionesus', 'Lasso', 'RandomForest']
organisms = [("Yeast100", "Yeast100"), ("Ecoli100","Ecoli100")]

fobj = open("job_params_large_networks.txt", "w")

params = ["1_2_3_4_5","6_7_8_9_10","11_12_13_14_15","16_17_18_19_20_21"] 

for inf_method in inf_methods:
    for organism in organisms:
        for i in range(1,21):
            for param in params:
                long_str = "/projects/p20519/roller_output/large_networks/{}/{}-insilico_size100_{} /projects/p20519/roller_output/large_networks/{}/{}-insilico_size100_{}_ /home/jjw036/Roller/data/gnw_insilico/network_data/{}/{}-{}_timeseries.tsv td_window num_{}\n".format(inf_method, organism[0], i, inf_method, organism[0], i, organism[1], organism[1], i, param)
                #long_str = "/projects/p20519/roller_output/gnw/{}/{}-insilico_size10_{} /projects/p20519/roller_output/gnw/{}/{}-insilico_size10_{}_ /home/jjw036/Roller/data/gnw_insilico/network_data/{}/{}-{}_timeseries.tsv td_window num_{}\n".format(inf_method, organism[0], i, inf_method, organism[0], i, organism[1], organism[1], i, param)
                fobj.write(long_str)

fobj.close()
