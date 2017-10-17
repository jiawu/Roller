inf_methods = ['RandomForest','Dionesus','Lasso']

fobj = open("job_params_invitro_sampling2.txt", "w")

intervals = [2,3,4]

for inf_method in inf_methods:
    for interval in intervals:
        max_window = int(round(14/interval))
        param_str = "{}_{}".format(max_window, int(round(max_window*(0.5))))
        params = [param_str]
        print(interval, params)
        for param in params:
            long_str = "/projects/p20519/roller_output/gardner_out/{}/ /projects/p20519/roller_output/gardner_out/{}/ /home/jjw036/Roller/data/invitro/gardner_{}_timeseries.tsv td_window num_{}\n".format(inf_method, inf_method, interval, param)
            fobj.write(long_str)

fobj.close()
