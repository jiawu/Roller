inf_methods = ['RandomForest', 'Dionesus','Lasso']
organisms = [("yeast", "Yeast", "Yeast10"), ("ecoli","Ecoli", "Ecoli10")]

fobj = open("job_params_high_sampling_3.txt", "w")

intervals = [10]
window_sizes = [30,40, 50, 60, 70]
min_maxes=['1_1', '1_2', '1_3', '1_4','1_5','1_6', '1_7', '1_8', '2_2', '2_3', '2_4', '2_5', '2_6', '2_7', '2_8', '3_3', '3_4', '3_5', '3_6', '3_7', '3_8']


for inf_method in inf_methods:
    for organism in organisms:
        for i in range(1,6):
            for interval in intervals:
                for window_size in window_sizes:
                    param_str = "{}_{}".format(window_size, min_maxes[0])
                    for min_max in min_maxes[1:]:
                        param_str = "{}_{}_{}".format(param_str, window_size, min_max)
                    max_window = int(1000/interval+1)
                    params = [param_str]
                    print(interval, params)
                    for param in params:
                        long_str = "/projects/p20519/roller_output/high_sampling/{}/{}-insilico_size10_{} /projects/p20519/roller_output/high_sampling/{}/{}-insilico_size10_{}_ /home/jjw036/Roller/data/gnw_insilico/high_sampling/{}/{}-{}_{}_timeseries.tsv td_window^min_lag^max_lag triplet_{}\n".format(inf_method, organism[0], i, inf_method, organism[0], i, organism[2], organism[1], i, interval, param)
                        fobj.write(long_str)

fobj.close()
