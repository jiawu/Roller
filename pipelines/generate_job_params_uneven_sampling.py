inf_methods = ['Dionesus', 'Lasso', 'RandomForest']
organisms = [("yeast", "Yeast", "Yeast10"), ("ecoli","Ecoli", "Ecoli10")]

fobj = open("job_params_even_sampling.txt", "w")

intervals = [10,30,50,100,200,333,500]
max_window = 7
test_window = 4

for inf_method in inf_methods:
    for organism in organisms:
        for i in range(1,21):
            param_str = "{}_{}".format(max_window, test_window)
            params = [param_str]
            print(params)
            for param in params:
                long_str = "/projects/p20519/roller_output/high_sampling/{}/{}-insilico_size10_{} /projects/p20519/roller_output/high_sampling/{}/{}-insilico_size10_{}_ /home/jjw036/Roller/data/gnw_insilico/high_sampling/{}/{}-{}_even_timeseries.tsv td_window num_{}\n".format(inf_method, organism[0], i, inf_method, organism[0], i, organism[2], organism[1], i,  param)
                fobj.write(long_str)

fobj.close()
