inf_methods = ['Lasso', 'RandomForest']
organisms = [("yeast", "Yeast"), ("ecoli","Ecoli")]

fobj = open("job_params_small_window_sizes.txt", "w")

params = ["1","2","3","4","5"] 

for inf_method in inf_methods:
    for organism in organisms:
        for i in range(1,21):
            for param in params:
                #long_str = "/projects/p20519/roller_output/large_networks/{}/{}-insilico_size100_{} /projects/p20519/roller_output/large_networks/{}/{}-insilico_size100_{}_ /home/jjw036/Roller/data/gnw_insilico/network_data/{}/{}-{}_timeseries.tsv td_window num_{}\n".format(inf_method, organism[0], i, inf_method, organism[0], i, organism[1], organism[1], i, param)
                long_str = "/projects/p20519/roller_output/gnw/{}/{}-insilico_size10_{} /projects/p20519/roller_output/gnw/{}/{}-insilico_size10_{}_ /home/jjw036/Roller/data/gnw_insilico/network_data/{}/{}-{}_timeseries.tsv td_window num_{}\n".format(inf_method, organism[0], i, inf_method, organism[0], i, organism[1], organism[1], i, param)
                fobj.write(long_str)

fobj.close()
