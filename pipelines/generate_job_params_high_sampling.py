inf_methods = ['Lasso']
organisms = [("yeast", "Yeast", "Yeast10"), ("ecoli","Ecoli", "Ecoli10")]

fobj = open("job_params_high_sampling.txt", "w")

intervals = [10,30,50,100,200,333,500]

for inf_method in inf_methods:
    for organism in organisms:
        for i in range(1,21):
            for interval in intervals:
                max_window = int(1000/interval+1)
                param_str = "{}_{}".format(max_window, int(max_window*(2/3)))
                params = [param_str]
                print(interval, params)
                for param in params:
                    long_str = "/projects/p20519/roller_output/high_sampling/{}/{}-insilico_size10_{} /projects/p20519/roller_output/high_sampling/{}/{}-insilico_size10_{}_ /home/jjw036/Roller/data/gnw_insilico/high_sampling/{}/{}-{}_{}_timeseries.tsv td_window num_{}\n".format(inf_method, organism[0], i, inf_method, organism[0], i, organism[2], organism[1], i, interval, param)
                    fobj.write(long_str)

fobj.close()
